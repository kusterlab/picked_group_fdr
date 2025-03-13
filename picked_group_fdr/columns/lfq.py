from __future__ import annotations

import collections
import itertools
from typing import List, Dict, Tuple
import logging

import numpy as np
from scipy.sparse.linalg import lsqr
from scipy.sparse import csr_matrix
import bottleneck as bn
from job_pool import JobPool

from .. import helpers
from .base import ProteinGroupColumns
from .sum_and_ibaq import _get_intensities, get_silac_channels
from .peptide_count import _unique_peptide_counts_per_experiment

# imports for typing
from .. import results
from .. import precursor_quant


logger = logging.getLogger(__name__)


class LFQIntensityColumns(ProteinGroupColumns):
    minPeptideRatiosLFQ: int
    stabilizeLargeRatiosLFQ: bool

    def __init__(
        self,
        minPeptideRatiosLFQ: int,
        stabilizeLargeRatiosLFQ: bool,
        numThreads: int = 1,
        protein_group_fdr_threshold: float = 0.01,
    ) -> None:
        self.minPeptideRatiosLFQ = minPeptideRatiosLFQ
        self.stabilizeLargeRatiosLFQ = stabilizeLargeRatiosLFQ
        self.numThreads = numThreads

        # only used for reporting number of protein groups at the given threshold
        self.protein_group_fdr_threshold = protein_group_fdr_threshold

    def is_valid(self, protein_group_results: results.ProteinGroupResults):
        return (
            len(protein_group_results.experiments) > 1
            and protein_group_results.num_tmt_channels <= 0
        )

    def append_headers(
        self,
        protein_group_results: results.ProteinGroupResults,
    ) -> None:
        silac_channels = get_silac_channels(protein_group_results.num_silac_channels)
        for experiment in protein_group_results.experiments:
            if protein_group_results.num_silac_channels > 0:
                for silacChannel in silac_channels:
                    protein_group_results.append_header(
                        "LFQ Intensity " + silacChannel + " " + experiment
                    )
            else:
                protein_group_results.append_header("LFQ Intensity " + experiment)

    def append_columns(
        self,
        protein_group_results: results.ProteinGroupResults,
        post_err_prob_cutoff: float,
    ) -> None:
        experiment_to_idx_map = protein_group_results.get_experiment_to_idx_map()

        logger.info("Doing quantification: MaxLFQ intensity")
        silac_channels = get_silac_channels(protein_group_results.num_silac_channels)
        num_silac_channels = len(silac_channels)

        if self.numThreads > 1:
            processingPool = JobPool(
                processes=self.numThreads,
                maxtasksperchild=10,
                max_jobs_queued=self.numThreads*3,
                write_progress_to_logger=True,
                total_jobs=len(protein_group_results),
            )

        allIntensities = list()
        for i, pgr in enumerate(protein_group_results):
            args = [
                pgr.precursorQuants,
                experiment_to_idx_map,
                post_err_prob_cutoff,
                self.minPeptideRatiosLFQ,
                self.stabilizeLargeRatiosLFQ,
                num_silac_channels,
            ]
            if self.numThreads > 1:
                processingPool.applyAsync(_getLFQIntensities, args)
            else:
                if i % 100 == 0:
                    logger.info(f"Processing protein {i}/{len(protein_group_results)}")
                allIntensities.append(_getLFQIntensities(*args))

        if self.numThreads > 1:
            allIntensities = processingPool.checkPool(printProgressEvery=100)

        proteinGroupCounts = np.zeros(
            len(experiment_to_idx_map) * max(1, num_silac_channels), dtype="int"
        )
        for i, (pgr, intensities) in enumerate(
            zip(protein_group_results, allIntensities)
        ):
            pgr.extend(intensities)

            if pgr.qValue < self.protein_group_fdr_threshold:
                proteinGroupCounts += np.array(
                    [1 if intensity > 0 else 0 for intensity in intensities]
                )

        logger.info(
            f"#Protein groups quantified ({self.protein_group_fdr_threshold*100:g}% protein group-level FDR, LFQ):"
        )
        for experiment, numProteinGroups in zip(
            experiment_to_idx_map.keys(),
            helpers.chunks(proteinGroupCounts, max(1, num_silac_channels)),
        ):
            if num_silac_channels > 0:
                for silacIdx, silacChannel in enumerate(silac_channels):
                    logger.info(
                        f"    {experiment} {silacChannel}: {numProteinGroups[silacIdx]}"
                    )
            else:
                logger.info(f"    {experiment}: {numProteinGroups[0]}")


def _getLFQIntensities(
    precursor_list: List[precursor_quant.PrecursorQuant],
    experiment_to_idx_map: Dict[str, int],
    postErrProbCutoff: float,
    minPeptideRatiosLFQ: int = 2,
    stabilizeLargeRatiosLFQ: bool = False,
    numSilacChannels: int = 0,
) -> List[float]:
    # in case of SILAC, we output LFQ intensites for each of the channels
    # only, not for the experiment (=sum over channels) itself, which is
    # what MQ does as well
    numExperiments = len(experiment_to_idx_map) * max(1, numSilacChannels)

    peptideIntensities, totalIntensity = _getPeptideIntensities(
        precursor_list,
        experiment_to_idx_map,
        postErrProbCutoff,
        numSilacChannels,
        numExperiments,
    )
    if len(peptideIntensities) == 0:
        return [0.0] * numExperiments

    logMedianPeptideRatios = _getLogMedianPeptideRatios(
        peptideIntensities, minPeptideRatiosLFQ
    )
    if len(logMedianPeptideRatios) == 0:
        return [0.0] * numExperiments

    if stabilizeLargeRatiosLFQ:
        logMedianPeptideRatios = _applyLargeRatioStabilization(
            logMedianPeptideRatios,
            precursor_list,
            experiment_to_idx_map,
            postErrProbCutoff,
            numSilacChannels,
        )

    matrix, vector = _buildLinearSystem(logMedianPeptideRatios, numExperiments)
    if len(vector) == 0:
        return [0.0] * numExperiments

    intensities = _solveLinearSystem(matrix, vector)
    intensities = _scaleEqualSum(intensities, totalIntensity)

    return intensities


def _getPeptideIntensities(
    precursor_list: List[precursor_quant.PrecursorQuant],
    experimentToIdxMap: Dict[str, int],
    postErrProbCutoff: float,
    numSilacChannels: int,
    numExperiments: int,
) -> Tuple[Dict[Tuple[str, int], List[float]], float]:
    """
    Collects all precursor intensities per experiment
    """

    def filterMissingAndUnidentified(p: precursor_quant.PrecursorQuant):
        return p.intensity > 0.0 and (
            helpers.is_mbr(p.post_err_prob) or p.post_err_prob <= postErrProbCutoff
        )

    precursor_list = filter(filterMissingAndUnidentified, precursor_list)

    # for each (peptide, charge, experiment, fraction) tuple, sort the
    # lowest PEP (= most confident PSM) on top
    def orderByPEP(p: precursor_quant.PrecursorQuant):
        return (
            p.peptide,
            p.charge,
            p.experiment,
            p.fraction,
            -1 * p.intensity,
            p.post_err_prob,
        )

    precursor_list = sorted(precursor_list, key=orderByPEP)

    peptideIntensities = collections.defaultdict(lambda: [0.0] * numExperiments)
    totalIntensity = 0.0
    prevExpFrac = (None, None)
    prevPrecursor = (None, None)
    for precursor in precursor_list:
        expIdx = experimentToIdxMap[precursor.experiment]
        # for each (peptide, charge, experiment, fraction) tuple, only
        # use the intensity of the PSM with the lowest PEP (= most
        # confident PSM)
        currPrecursor = (precursor.peptide, precursor.charge)
        currExpFrac = (precursor.experiment, precursor.fraction)
        if prevExpFrac != currExpFrac or prevPrecursor != currPrecursor:
            if numSilacChannels > 0:
                for silacIdx, silacIntensity in enumerate(precursor.silac_intensities):
                    silacExpIdx = expIdx * numSilacChannels + silacIdx
                    peptideIntensities[currPrecursor][silacExpIdx] += silacIntensity
                    totalIntensity += silacIntensity
            else:
                peptideIntensities[currPrecursor][expIdx] += precursor.intensity
                totalIntensity += precursor.intensity
            prevExpFrac = currExpFrac
            prevPrecursor = currPrecursor
    return peptideIntensities, totalIntensity


def _applyLargeRatioStabilization(
    logMedianPeptideRatios: Dict[Tuple[int, int], float],
    precursor_list: List[precursor_quant.PrecursorQuant],
    experimentToIdxMap: Dict[str, int],
    postErrProbCutoff: float,
    numSilacChannels: int,
) -> Dict[Tuple[int, int], float]:
    summedIntensities = _get_intensities(
        precursor_list, experimentToIdxMap, postErrProbCutoff, numSilacChannels
    )
    if numSilacChannels > 0:
        # removes the column intensities per experiment with the SILAC
        # channels summed up
        del summedIntensities[:: (numSilacChannels + 1)]

    peptideCounts = _unique_peptide_counts_per_experiment(
        precursor_list, experimentToIdxMap, postErrProbCutoff
    )
    peptideCounts = np.repeat(peptideCounts, numSilacChannels)
    for (i, j), (si1, si2), (pc1, pc2) in zip(
        itertools.combinations(range(len(summedIntensities)), 2),
        itertools.combinations(summedIntensities, 2),
        itertools.combinations(peptideCounts, 2),
    ):
        if pc1 == 0 or pc2 == 0 or (i, j) not in logMedianPeptideRatios:
            continue

        peptideCountRatio = _getMaxRatio(pc1, pc2)
        if peptideCountRatio > 5:
            logMedianPeptideRatios[(i, j)] = np.log(si1 / si2)
        elif peptideCountRatio > 2.5:
            w = (peptideCountRatio - 2.5) / 2.5
            logMedianPeptideRatios[(i, j)] = (
                w * np.log(si1 / si2) + (1 - w) * logMedianPeptideRatios[(i, j)]
            )

    return logMedianPeptideRatios


def _getLogMedianPeptideRatios(
    peptideIntensities: Dict[Tuple[str, int], List[float]], minPeptideRatiosLFQ: int
) -> Dict[Tuple[int, int], float]:
    """
    :param peptideIntensities: rows are peptides, columns are experiments
    :param minPeptideRatiosLFQ: minimum valid ratios needed to perform LFQ
    returns: dictionary of (sample_i, sample_j) -> log(median(ratios))
    """
    intensityMatrix = np.array(list(peptideIntensities.values()))

    nonzeros_per_column = np.count_nonzero(intensityMatrix, axis=0)
    valid_columns = np.argwhere(nonzeros_per_column >= minPeptideRatiosLFQ).flatten()

    valid_counts_matrix = (intensityMatrix[:, valid_columns] > 0).astype(int)
    valid_vals = np.dot(valid_counts_matrix.T, valid_counts_matrix)

    # replace zeroes by NaNs to automatically filter ratios with one missing value
    intensityMatrix[intensityMatrix == 0] = np.nan
    cols = [(idx, i, intensityMatrix[:, i]) for idx, i in enumerate(valid_columns)]
    peptideRatios, experimentPairs = list(), list()
    # this for loop is the reason for quadratic runtime increase for more samples
    # (30 seconds per protein for 2400 samples and 600 peptides)
    for (idx_i, i, col_i), (idx_j, j, col_j) in itertools.combinations(cols, 2):
        if valid_vals[idx_i, idx_j] < minPeptideRatiosLFQ:
            continue

        ratios = col_i / col_j
        # vectorization of computing medians is 25% faster than computing
        # them in a loop here but requires more memory
        peptideRatios.append(bn.nanmedian(ratios))
        experimentPairs.append((i, j))

    return dict(zip(experimentPairs, np.log(peptideRatios)))


def _getMaxRatio(pc1: float, pc2: float) -> float:
    assert pc1 > 0 and pc2 > 0
    if pc1 < pc2:
        return pc2 / pc1
    else:
        return pc1 / pc2


def _buildLinearSystem(
    logMedianPeptideRatios: Dict[Tuple[int, int], float], numExperiments: int
) -> Tuple[np.array, np.array]:
    numRatios = len(logMedianPeptideRatios)
    rows, cols, vals = list(), list(), list()
    seenColumns = set()
    for idx, (i, j) in enumerate(logMedianPeptideRatios.keys()):
        rows.append(idx)
        cols.append(i)
        vals.append(1.0)
        seenColumns.add(i)

        rows.append(idx)
        cols.append(j)
        vals.append(-1.0)
        seenColumns.add(j)

    # anchor the (log) average and experiments without ratios to 0,
    # otherwise the matrix is ill-conditioned and least squares sometimes
    # does not converge
    num_non_zero_columns = len(seenColumns)
    zero_columns = [x for x in range(numExperiments) if x not in seenColumns]
    num_zero_columns = len(zero_columns)

    rows.extend([numRatios] * num_non_zero_columns)
    cols.extend(list(seenColumns))
    vals.extend([1.0] * num_non_zero_columns)

    rows.extend(range(numRatios + 1, numRatios + 1 + num_zero_columns))
    cols.extend(zero_columns)
    vals.extend([1.0] * num_zero_columns)

    matrix = csr_matrix((vals, (rows, cols)), shape=(max(rows) + 1, max(cols) + 1))
    vector = np.array(
        list(logMedianPeptideRatios.values()) + [0.0] * (num_zero_columns + 1)
    )

    # np.savetxt("matrix.csv", matrix, delimiter="\t")
    # np.savetxt("vector.csv", vector, delimiter="\t")

    return matrix, vector


def _solveLinearSystem(matrix: csr_matrix, vector: np.array) -> np.array:
    intensities = np.exp(lsqr(matrix, vector)[0])

    zero_columns = np.array(matrix.astype(bool).sum(axis=0))[0] == 1
    intensities[zero_columns] = 0.0
    return intensities


def _scaleEqualSum(intensities: np.array, totalIntensity: float) -> List[float]:
    """
    scale intensities such that the total sum stays the same as it was
    before

    :param intensities: array of intensities per sample
    :param totalIntensity: sum of intensities prior to method application
    :returns: array of intensities scaled such that sum equals
              totalIntensity
    """
    intensity_sum = np.sum(intensities)
    if intensity_sum > 0:
        return ((totalIntensity / intensity_sum) * intensities).tolist()
    else:
        return intensities.tolist()
