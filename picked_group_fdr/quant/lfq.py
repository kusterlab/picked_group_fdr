import collections
import itertools
import warnings
from typing import List, Dict, Tuple
import logging

import numpy as np
from scipy.sparse.linalg import lsqr
from scipy.sparse import csr_matrix
import bottleneck as bn

from .. import helpers
from .base import ProteinGroupColumns
from .sum_and_ibaq import SummedIntensityAndIbaqColumns
from .peptide_count import UniquePeptideCountColumns
from ..utils import multiprocessing_pool as pool

# imports for typing
from ..results import ProteinGroupResults
from .precursor_quant import PrecursorQuant


logger = logging.getLogger(__name__)


class LFQIntensityColumns(ProteinGroupColumns):
    silacChannels: List[str]
    minPeptideRatiosLFQ: int
    stabilizeLargeRatiosLFQ: bool

    def __init__(self, silacChannels: List[str], minPeptideRatiosLFQ: int,
                 stabilizeLargeRatiosLFQ: bool, numThreads: int = 1) -> None:
        self.silacChannels = silacChannels
        self.minPeptideRatiosLFQ = minPeptideRatiosLFQ
        self.stabilizeLargeRatiosLFQ = stabilizeLargeRatiosLFQ
        self.numThreads = numThreads

    def append_headers(self, proteinGroupResults: ProteinGroupResults,
                       experiments: int) -> None:
        if len(experiments) <= 1:
            return

        numSilacChannels = len(self.silacChannels)
        for experiment in experiments:
            if numSilacChannels > 0:
                for silacChannel in self.silacChannels:
                    proteinGroupResults.append_header(
                        'LFQ Intensity ' + silacChannel + ' ' + experiment)
            else:
                proteinGroupResults.append_header(
                    'LFQ Intensity ' + experiment)

    def append_columns(self, proteinGroupResults: ProteinGroupResults,
                       experimentToIdxMap: Dict[str, int],
                       postErrProbCutoff: float) -> None:
        if len(experimentToIdxMap) <= 1:
            return

        logger.info("Doing quantification: MaxLFQ intensity")
        numSilacChannels = len(self.silacChannels)
        
        if self.numThreads > 1:
            processingPool = pool.JobPool(processes=self.numThreads, maxtasksperchild=10)
        
        allIntensities = list()
        for i, pgr in enumerate(proteinGroupResults):
            args = [pgr.precursorQuants, experimentToIdxMap, postErrProbCutoff, 
                    self.minPeptideRatiosLFQ, self.stabilizeLargeRatiosLFQ, numSilacChannels]
            if self.numThreads > 1:
                processingPool.applyAsync(_getLFQIntensities, args)
            else:
                if i % 100 == 0:
                    logger.info(f"Processing protein {i}/{len(proteinGroupResults)}")
                allIntensities.append(_getLFQIntensities(*args))

        if self.numThreads > 1:
            allIntensities = processingPool.checkPool(printProgressEvery=100)

        proteinGroupCounts = np.zeros(
            len(experimentToIdxMap) * max(1, numSilacChannels), dtype='int')
        for i, (pgr, intensities) in enumerate(zip(proteinGroupResults, allIntensities)):
            pgr.extend(intensities)

            if pgr.qValue < 0.01:
                proteinGroupCounts += np.array(
                    [1 if intensity > 0 else 0 for intensity in intensities])

        logger.info("#Protein groups quantified (1% protein group-level FDR, LFQ):")
        for experiment, numProteinGroups in zip(
                experimentToIdxMap.keys(),
                helpers.chunks(proteinGroupCounts, max(1, numSilacChannels))):
            if numSilacChannels > 0:
                for silacIdx, silacChannel in enumerate(self.silacChannels):
                    logger.info(f"    {experiment} {silacChannel}: {numProteinGroups[silacIdx]}")
            else:
                logger.info(f"    {experiment}: {numProteinGroups[0]}")


def _getLFQIntensities(peptideIntensityList: List[PrecursorQuant],
                       experimentToIdxMap: Dict[str, int],
                       postErrProbCutoff: float,
                       minPeptideRatiosLFQ: int = 2,
                       stabilizeLargeRatiosLFQ: bool = False,
                       numSilacChannels: int = 0) -> List[float]:
    # in case of SILAC, we output LFQ intensites for each of the channels
    # only, not for the experiment (=sum over channels) itself, which is
    # what MQ does as well
    numExperiments = len(experimentToIdxMap) * max(1, numSilacChannels)

    peptideIntensities, totalIntensity = _getPeptideIntensities(
        peptideIntensityList, experimentToIdxMap, postErrProbCutoff,
        numSilacChannels, numExperiments)
    if len(peptideIntensities) == 0:
        return [0.0] * numExperiments

    logMedianPeptideRatios = _getLogMedianPeptideRatios(
        peptideIntensities, minPeptideRatiosLFQ)
    if len(logMedianPeptideRatios) == 0:
        return [0.0] * numExperiments

    if stabilizeLargeRatiosLFQ:
        logMedianPeptideRatios = _applyLargeRatioStabilization(
            logMedianPeptideRatios, peptideIntensityList,
            experimentToIdxMap, postErrProbCutoff, numSilacChannels)

    matrix, vector = _buildLinearSystem(
        logMedianPeptideRatios, numExperiments)
    if len(vector) == 0:
        return [0.0] * numExperiments

    intensities = _solveLinearSystem(matrix, vector)
    intensities = _scaleEqualSum(intensities, totalIntensity)

    return intensities


def _getPeptideIntensities(
        peptideIntensityList: List[PrecursorQuant],
        experimentToIdxMap: Dict[str, int],
        postErrProbCutoff: float,
        numSilacChannels: int,
        numExperiments: int) \
        -> Tuple[Dict[Tuple[str, int], List[float]], float]:
    """
    Collects all precursor intensities per experiment
    """
    def filterMissingAndUnidentified(p): 
        return (p.intensity > 0.0 and
            (helpers.isMbr(p.postErrProb) or p.postErrProb <= postErrProbCutoff))
    peptideIntensityList = filter(filterMissingAndUnidentified, peptideIntensityList)
    
    # for each (peptide, charge, experiment, fraction) tuple, sort the
    # lowest PEP (= most confident PSM) on top
    def orderByPEP(p): return (p.peptide, p.charge, p.experiment,
                               p.fraction, -1 * p.intensity, p.postErrProb)
    peptideIntensityList = sorted(peptideIntensityList, key=orderByPEP)
    
    peptideIntensities = collections.defaultdict(
        lambda: [0.0] * numExperiments)
    totalIntensity = 0.0
    prevExpFrac = (None, None)
    prevPrecursor = (None, None)
    for precursor in peptideIntensityList:
        expIdx = experimentToIdxMap[precursor.experiment]
        # for each (peptide, charge, experiment, fraction) tuple, only
        # use the intensity of the PSM with the lowest PEP (= most
        # confident PSM)
        currPrecursor = (precursor.peptide, precursor.charge)
        currExpFrac = (precursor.experiment, precursor.fraction)
        if prevExpFrac != currExpFrac or prevPrecursor != currPrecursor:
            if numSilacChannels > 0:
                for silacIdx, silacIntensity in enumerate(
                        precursor.silacIntensities):
                    silacExpIdx = expIdx * numSilacChannels + silacIdx
                    peptideIntensities[currPrecursor][silacExpIdx] += (
                            silacIntensity)
                    totalIntensity += silacIntensity
            else:
                peptideIntensities[currPrecursor][expIdx] += (
                        precursor.intensity)
                totalIntensity += precursor.intensity
            prevExpFrac = currExpFrac
            prevPrecursor = currPrecursor
    return peptideIntensities, totalIntensity


def _applyLargeRatioStabilization(
        logMedianPeptideRatios: Dict[Tuple[int, int], float],
        peptideIntensityList: List[PrecursorQuant],
        experimentToIdxMap: Dict[str, int],
        postErrProbCutoff: float,
        numSilacChannels: int) -> Dict[Tuple[int, int], float]:
    summedIntensities = SummedIntensityAndIbaqColumns.getIntensities(
        peptideIntensityList, experimentToIdxMap, postErrProbCutoff,
        numSilacChannels)
    if numSilacChannels > 0:
        # removes the column intensities per experiment with the SILAC
        # channels summed up
        del summedIntensities[::(numSilacChannels + 1)]

    peptideCounts = \
        UniquePeptideCountColumns.uniquePeptideCountsPerExperiment(
            peptideIntensityList, experimentToIdxMap, postErrProbCutoff)
    peptideCounts = np.repeat(peptideCounts, numSilacChannels)
    for (i, j), (si1, si2), (pc1, pc2) in zip(
            itertools.combinations(range(len(summedIntensities)), 2),
            itertools.combinations(summedIntensities, 2),
            itertools.combinations(peptideCounts, 2)):
        if pc1 == 0 or pc2 == 0 or (i, j) not in logMedianPeptideRatios:
            continue

        peptideCountRatio = _getMaxRatio(pc1, pc2)
        if peptideCountRatio > 5:
            logMedianPeptideRatios[(i, j)] = np.log(si1 / si2)
        elif peptideCountRatio > 2.5:
            w = (peptideCountRatio - 2.5) / 2.5
            logMedianPeptideRatios[(i, j)] = (
                    w * np.log(si1 / si2) +
                    (1 - w) * logMedianPeptideRatios[(i, j)])

    return logMedianPeptideRatios


def _getLogMedianPeptideRatios(
        peptideIntensities: Dict[Tuple[str, int], List[float]],
        minPeptideRatiosLFQ: int) -> Dict[Tuple[int, int], float]:
    """
    :param peptideIntensities: rows are peptides, columns are experiments
    :param minPeptideRatiosLFQ: minimum valid ratios needed to perform LFQ
    returns: dictionary of (sample_i, sample_j) -> log(median(ratios))
    """        
    intensityMatrix = np.array(list(peptideIntensities.values()))
    
    nonzeros_per_column = np.count_nonzero(intensityMatrix, axis=0)
    valid_columns = np.argwhere(
        nonzeros_per_column >= minPeptideRatiosLFQ).flatten()
    
    valid_counts_matrix = (intensityMatrix[:, valid_columns] > 0).astype(int)
    valid_vals = np.dot(valid_counts_matrix.T, valid_counts_matrix)
    
    # replace zeroes by NaNs to automatically filter ratios with one missing value
    intensityMatrix[intensityMatrix == 0] = np.nan
    columns = [(idx, i, intensityMatrix[:, i]) for idx, i in enumerate(valid_columns)]
    peptideRatios, experimentPairs = list(), list()
    for ((idx_i, i, col_i), (idx_j, j, col_j)) in itertools.combinations(
            columns, 2):
        if valid_vals[idx_i, idx_j] < minPeptideRatiosLFQ:
            continue
        
        ratios = col_i / col_j
        # vectorization of computing medians is 25% faster than computing 
        # them in a loop here but requires more memory
        peptideRatios.append(bn.nanmedian(ratios))
        experimentPairs.append((i,j))
    
    return dict(zip(experimentPairs, np.log(peptideRatios)))


def _getMaxRatio(pc1: float, pc2: float) -> float:
    assert pc1 > 0 and pc2 > 0
    if pc1 < pc2:
        return pc2 / pc1
    else:
        return pc1 / pc2


def _buildLinearSystem(
        logMedianPeptideRatios: Dict[Tuple[int, int], float],
        numExperiments: int) -> Tuple[np.array, np.array]:
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
    
    rows.extend(range(numRatios+1, numRatios+1 + num_zero_columns))
    cols.extend(zero_columns)
    vals.extend([1.0] * num_zero_columns)
    
    matrix = csr_matrix((vals, (rows, cols)), shape=(max(rows)+1, max(cols)+1))
    vector = np.array(list(logMedianPeptideRatios.values()) + [0.0]*(num_zero_columns+1))
    
    # np.savetxt("matrix.csv", matrix, delimiter="\t")
    # np.savetxt("vector.csv", vector, delimiter="\t")

    return matrix, vector


def _solveLinearSystem(matrix: csr_matrix,
                       vector: np.array) -> np.array:
    intensities = np.exp(lsqr(matrix, vector)[0])
    
    zero_columns = np.array(matrix.astype(bool).sum(axis=0))[0] == 1
    intensities[zero_columns] = 0.0
    return intensities


def _scaleEqualSum(intensities: np.array,
                   totalIntensity: float) -> List[float]:
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

