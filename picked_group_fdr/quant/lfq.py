import collections
import itertools
import warnings
from typing import List, Dict, Tuple
import logging

import numpy as np
from scipy.sparse.linalg import lsqr
from scipy.sparse import csr_matrix

from .. import helpers
from .base import ProteinGroupColumns
from .sum_and_ibaq import SummedIntensityAndIbaqColumns
from .peptide_count import UniquePeptideCountColumns
# imports for typing
from ..results import ProteinGroupResults
from .precursor_quant import PrecursorQuant


logger = logging.getLogger(__name__)


class LFQIntensityColumns(ProteinGroupColumns):
    silacChannels: List[str]
    minPeptideRatiosLFQ: int
    stabilizeLargeRatiosLFQ: bool

    def __init__(self, silacChannels: List[str], minPeptideRatiosLFQ: int,
                 stabilizeLargeRatiosLFQ: bool) -> None:
        self.silacChannels = silacChannels
        self.minPeptideRatiosLFQ = minPeptideRatiosLFQ
        self.stabilizeLargeRatiosLFQ = stabilizeLargeRatiosLFQ

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
        proteinGroupCounts = np.zeros(
            len(experimentToIdxMap) * max(1, numSilacChannels), dtype='int')
        for i, pgr in enumerate(proteinGroupResults):
            if i % 100 == 0:
                logger.info(f"Processing protein {i}/{len(proteinGroupResults)}")
            intensities = self._getLFQIntensities(
                pgr.precursorQuants, experimentToIdxMap, postErrProbCutoff)
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

    def _getLFQIntensities(self, peptideIntensityList: List[PrecursorQuant],
                           experimentToIdxMap: Dict[str, int],
                           postErrProbCutoff: float) -> List[float]:
        numSilacChannels = len(self.silacChannels)

        # in case of SILAC, we output LFQ intensites for each of the channels
        # only, not for the experiment (=sum over channels) itself, which is
        # what MQ does as well
        numExperiments = len(experimentToIdxMap) * max(1, numSilacChannels)

        peptideIntensities, totalIntensity = self._getPeptideIntensities(
            peptideIntensityList, experimentToIdxMap, postErrProbCutoff,
            numSilacChannels, numExperiments)
        if len(peptideIntensities) == 0:
            return [0.0] * numExperiments

        logMedianPeptideRatios = self._getLogMedianPeptideRatios(
            peptideIntensities, self.minPeptideRatiosLFQ)
        if len(logMedianPeptideRatios) == 0:
            return [0.0] * numExperiments

        if self.stabilizeLargeRatiosLFQ:
            logMedianPeptideRatios = self._applyLargeRatioStabilization(
                logMedianPeptideRatios, peptideIntensityList,
                experimentToIdxMap, postErrProbCutoff, numSilacChannels)

        matrix, vector = self._buildLinearSystem(
            logMedianPeptideRatios, numExperiments)
        if len(vector) == 0:
            return [0.0] * numExperiments

        intensities = self._solveLinearSystem(matrix, vector)
        intensities = self._scaleEqualSum(intensities, totalIntensity)

        return intensities

    def _getPeptideIntensities(
            self, peptideIntensityList: List[PrecursorQuant],
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
            self,
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

            peptideCountRatio = self._getMaxRatio(pc1, pc2)
            if peptideCountRatio > 5:
                logMedianPeptideRatios[(i, j)] = np.log(si1 / si2)
            elif peptideCountRatio > 2.5:
                w = (peptideCountRatio - 2.5) / 2.5
                logMedianPeptideRatios[(i, j)] = (
                        w * np.log(si1 / si2) +
                        (1 - w) * logMedianPeptideRatios[(i, j)])

        return logMedianPeptideRatios

    def _getLogMedianPeptideRatios(
            self,
            peptideIntensities: Dict[Tuple[str, int], List[float]],
            minPeptideRatiosLFQ: int) -> Dict[Tuple[int, int], float]:
        """
        :param peptideIntensities: rows are peptides, columns are experiments
        :param minPeptideRatiosLFQ: minimum valid ratios needed to perform LFQ
        """
        peptideRatios, experimentPairs = list(), list()
        maxRatiosLength = 0
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            intensityMatrix = np.array(
                list(peptideIntensities.values()))
            nonzeros_per_column = np.count_nonzero(intensityMatrix, axis=0)
            valid_columns = np.argwhere(
                nonzeros_per_column >= minPeptideRatiosLFQ).flatten()
            columns = [(i, intensityMatrix[:, i]) for i in valid_columns]
            for ((i, col_i), (j, col_j)) in itertools.combinations(
                    columns, 2):
                ratios = col_i / col_j
                # remove ratios where one of the intensities was 0
                ratios = ratios[(np.isfinite(ratios)) & (ratios != 0)]
                if len(ratios) >= minPeptideRatiosLFQ:
                    peptideRatios.append(ratios)
                    experimentPairs.append((i,j))
                    maxRatiosLength = max(len(ratios), maxRatiosLength)

        # ratio lists are not always of equal size, pad those with NaNs
        peptideRatiosArray = self._createJaggedArray(peptideRatios, maxRatiosLength)
        logPeptideRatios = np.log(np.nanmedian(peptideRatiosArray, axis=1))

        peptideRatiosDict = {p: r for p, r in zip(experimentPairs, logPeptideRatios)}
        return peptideRatiosDict

    def _createJaggedArray(self, peptideRatios: List[np.array], maxRatiosLength: int):
        peptideRatiosArray = np.empty((len(peptideRatios), maxRatiosLength))
        peptideRatiosArray[:] = np.nan
        for idx, row in enumerate(peptideRatios):
            peptideRatiosArray[idx, :len(row)] = row
        return peptideRatiosArray

    def _getMaxRatio(self, pc1: float, pc2: float) -> float:
        assert pc1 > 0 and pc2 > 0
        if pc1 < pc2:
            return pc2 / pc1
        else:
            return pc1 / pc2

    def _buildLinearSystem(
            self, logMedianPeptideRatios: Dict[Tuple[int, int], float],
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

    def _solveLinearSystem(self, matrix: csr_matrix,
                           vector: np.array) -> np.array:
        intensities = np.exp(lsqr(matrix, vector)[0])
        
        zero_columns = np.array(matrix.astype(bool).sum(axis=0))[0] == 1
        intensities[zero_columns] = 0.0
        return intensities

    def _scaleEqualSum(self, intensities: np.array,
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

