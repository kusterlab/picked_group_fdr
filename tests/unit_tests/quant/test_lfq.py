import pytest
import numpy as np

from picked_group_fdr.quant.lfq import LFQIntensityColumns


class TestGetLFQIntensities:
  def test_getLFQIntensities_3files(self, pq, peptideIntensityListThreeFiles, experimentToIdxMap):
    np.testing.assert_almost_equal(pq._getLFQIntensities(peptideIntensityListThreeFiles, experimentToIdxMap, 0.1), [25.0, 10.0, 5.0])

  def test_getLFQIntensities_5files(self, pq, peptideIntensityListFiveFiles, experimentToIdxMapFiveFiles):
    np.testing.assert_almost_equal(pq._getLFQIntensities(peptideIntensityListFiveFiles, experimentToIdxMapFiveFiles, 0.1), [25.0, 10.0, 5.0, 0.0, 7.0])

  def test_getLFQIntensities_missingValuesZeroes(self, pq, peptideIntensityListMissingValuesZeroes, experimentToIdxMap):
    np.testing.assert_almost_equal(pq._getLFQIntensities(peptideIntensityListMissingValuesZeroes, experimentToIdxMap, 0.1), [0.0, 0.0, 0.0])

  def test_getLFQIntensities_missingValuesNans(self, pq, peptideIntensityListMissingValuesNans, experimentToIdxMap):
    np.testing.assert_almost_equal(pq._getLFQIntensities(peptideIntensityListMissingValuesNans, experimentToIdxMap, 0.1), [0.0, 0.0, 0.0])

  def test_getLFQIntensities_missingValuesNans_minPeptideRatiosLFQ2(self, pq2, peptideIntensityListMissingValuesNans, experimentToIdxMap):
    np.testing.assert_almost_equal(pq2._getLFQIntensities(peptideIntensityListMissingValuesNans, experimentToIdxMap, 0.1), [0.0, 0.0, 0.0])

  def test_getLFQIntensities_2peptides(self, pq2, peptideIntensityListMultiplePeptides, experimentToIdxMap):
    np.testing.assert_almost_equal(pq2._getLFQIntensities(peptideIntensityListMultiplePeptides, experimentToIdxMap, 0.1), [40.0, 20.0, 10.0])

  def test_getLFQIntensities_2peptidesJagged(self, pq2, peptideIntensityListMultiplePeptidesJagged, experimentToIdxMap):
    np.testing.assert_almost_equal(pq2._getLFQIntensities(peptideIntensityListMultiplePeptidesJagged, experimentToIdxMap, 0.1), [50.0, 20.0, 0.0])

  def test_getLFQIntensities_tooFewPeptideRatios(self, pq3, peptideIntensityListMultiplePeptides, experimentToIdxMap):
    np.testing.assert_almost_equal(pq3._getLFQIntensities(peptideIntensityListMultiplePeptides, experimentToIdxMap, 0.1), [0.0, 0.0, 0.0])

  def test_getLFQIntensities_duplicatePeptides(self, pq, peptideIntensityListDuplicatePeptides, experimentToIdxMap):
    np.testing.assert_almost_equal(pq._getLFQIntensities(peptideIntensityListDuplicatePeptides, experimentToIdxMap, 0.1), [25.0, 10.0, 5.0])
  
  def test_getLFQIntensities_multipleCharges(self, pq2, peptideIntensityListMultipleCharges, experimentToIdxMap):
    """MaxLFQ considers different charge states of the same peptide as two different peptides for the sake of minPeptideRatiosLFQ"""
    np.testing.assert_almost_equal(pq2._getLFQIntensities(peptideIntensityListMultipleCharges, experimentToIdxMap, 0.1), [40.0, 20.0, 10.0])
  
  def test_getLFQIntensities_2fractions(self, pq, peptideIntensityListMultipleFractions, experimentToIdxMap):
    np.testing.assert_almost_equal(pq._getLFQIntensities(peptideIntensityListMultipleFractions, experimentToIdxMap, 0.1), [35.0, 10.0, 5.0])
  
  def test_getLFQIntensities_silac(self, peptideIntensityListSilac, experimentToIdxMap):
    pq = LFQIntensityColumns(silacChannels=['l', 'm', 'h'], minPeptideRatiosLFQ=1, stabilizeLargeRatiosLFQ=False)
    np.testing.assert_almost_equal(pq._getLFQIntensities(peptideIntensityListSilac, experimentToIdxMap, 0.1), [13.0, 17.0, 5.0, 3.0, 4.0, 3.0, 1.5, 2.5, 1.0])


class TestStabilizeLargeRatios:
  def test_getMaxRatio(self, pq):
    assert pq._getMaxRatio(5, 2) == 5 / 2
    
  def test_getMaxRatio_reverse(self, pq):
    assert pq._getMaxRatio(2, 5) == 5 / 2
  
  def test_getMaxRatio_both_zero(self, pq):
    with pytest.raises(AssertionError):
      pq._getMaxRatio(0, 0)
  
  def test_getLFQIntensities_3files(self, pq, peptideIntensityListThreeFiles, experimentToIdxMap):
    np.testing.assert_almost_equal(pq._getLFQIntensities(peptideIntensityListThreeFiles, experimentToIdxMap, 0.1), [25.0, 10.0, 5.0])


class TestGetPeptideIntensities:
  def test_getPeptideIntensities(self, pq, peptideIntensityListThreeFiles, experimentToIdxMap):
    postErrProbCutoff = 0.1
    numSilacChannels = 0
    numExperiments = 3
    np.testing.assert_almost_equal(pq._getPeptideIntensities(peptideIntensityListThreeFiles, experimentToIdxMap, postErrProbCutoff, numSilacChannels, numExperiments)[0][('_APEPTIDE_', 2)], [25.0, 10.0, 5.0])


class TestCalculateRatios:
  def test_getLogMedianPeptideRatios(self, pq, peptideIntensities):
    assert pq._getLogMedianPeptideRatios(peptideIntensities, 1) == {(0,1) : np.log(2.5), (0, 2): np.log(5.0), (1, 2): np.log(2.0)}


@pytest.fixture
def pq():
  return LFQIntensityColumns(silacChannels=[], minPeptideRatiosLFQ=1, stabilizeLargeRatiosLFQ=False)
  

@pytest.fixture
def pq2():
  return LFQIntensityColumns(silacChannels=[], minPeptideRatiosLFQ=2, stabilizeLargeRatiosLFQ=False)

@pytest.fixture
def pq3():
  return LFQIntensityColumns(silacChannels=[], minPeptideRatiosLFQ=3, stabilizeLargeRatiosLFQ=False)

@pytest.fixture
def peptideIntensities():
  return {('_APEPTIDE_', 2): [25.0, 10.0, 5.0]}
