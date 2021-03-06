import pytest
import numpy as np

from picked_group_fdr.quant.sum_and_ibaq import SummedIntensityAndIbaqColumns as pq
from picked_group_fdr.quantification import retainOnlyIdentifiedPrecursors


class TestGetIntensities:
  def test_getIntensities_3files(self, peptideIntensityListThreeFiles, experimentToIdxMap):
    np.testing.assert_almost_equal(pq.getIntensities(peptideIntensityListThreeFiles, experimentToIdxMap, 0.1, numSilacChannels=0), [25.0, 10.0, 5.0])

  def test_getIntensities_5files(self, peptideIntensityListFiveFiles, experimentToIdxMapFiveFiles):
    np.testing.assert_almost_equal(pq.getIntensities(peptideIntensityListFiveFiles, experimentToIdxMapFiveFiles, 0.1, numSilacChannels=0), [25.0, 10.0, 5.0, 0.0, 7.0])

  def test_getIntensities_missingValuesZeroes(self, peptideIntensityListMissingValuesZeroes, experimentToIdxMap):  
    np.testing.assert_almost_equal(pq.getIntensities(peptideIntensityListMissingValuesZeroes, experimentToIdxMap, 0.1, numSilacChannels=0), [25.0, 0.0, 0.0])

  def test_getIntensities_missingValuesNans(self, peptideIntensityListMissingValuesNans, experimentToIdxMap):
    np.testing.assert_almost_equal(pq.getIntensities(peptideIntensityListMissingValuesNans, experimentToIdxMap, 0.1, numSilacChannels=0), [25.0, 0.0, 0.0])

  def test_getIntensities_2peptides(self, peptideIntensityListMultiplePeptides, experimentToIdxMap):
    np.testing.assert_almost_equal(pq.getIntensities(peptideIntensityListMultiplePeptides, experimentToIdxMap, 0.1, numSilacChannels=0), [40.0, 20.0, 10.0])

  def test_getIntensities_duplicatePeptides(self, peptideIntensityListDuplicatePeptides, experimentToIdxMap):
    np.testing.assert_almost_equal(pq.getIntensities(peptideIntensityListDuplicatePeptides, experimentToIdxMap, 0.1, numSilacChannels=0), [35.0, 10.0, 5.0])

  def test_getIntensities_2fractions(self, peptideIntensityListMultipleFractions, experimentToIdxMap):
    np.testing.assert_almost_equal(pq.getIntensities(peptideIntensityListMultipleFractions, experimentToIdxMap, 0.1, numSilacChannels=0), [35.0, 10.0, 5.0])

  def test_getIntensities_silac(self, peptideIntensityListSilac, experimentToIdxMap):
    np.testing.assert_almost_equal(pq.getIntensities(peptideIntensityListSilac, experimentToIdxMap, 0.1, numSilacChannels = 3), [35.0, 13.0, 17.0, 5.0, 10.0, 3.0, 4.0, 3.0, 5.0, 1.5, 2.5, 1.0])
  
  def test_getIntensities_matchBetweenRuns(self, peptideIntensityListMatchBetweenRuns, experimentToIdxMap):
    postErrProbCutoff = 0.1
    peptideIntensityList = retainOnlyIdentifiedPrecursors(peptideIntensityListMatchBetweenRuns, postErrProbCutoff)
    np.testing.assert_almost_equal(pq.getIntensities(peptideIntensityList, experimentToIdxMap, postErrProbCutoff, numSilacChannels=0), [25.0, 10.0, 5.0])

  def test_getIntensities_matchBetweenRunsWithUnidentified(self, peptideIntensityListMatchBetweenRunsWithUnidentified, experimentToIdxMap):
    postErrProbCutoff = 0.1
    peptideIntensityList = retainOnlyIdentifiedPrecursors(peptideIntensityListMatchBetweenRunsWithUnidentified, postErrProbCutoff)
    np.testing.assert_almost_equal(pq.getIntensities(peptideIntensityList, experimentToIdxMap, postErrProbCutoff, numSilacChannels=0), [0.0, 0.0, 0.0])

