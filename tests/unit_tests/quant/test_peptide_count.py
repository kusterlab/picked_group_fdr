import pytest
import numpy as np

from picked_group_fdr.quant.peptide_count import UniquePeptideCountColumns as pq


class TestUniquePeptideCountsPerExperiment:
  def test_uniquePeptideCountsPerExperiment_3files(self, peptideIntensityListThreeFiles, experimentToIdxMap):
    np.testing.assert_equal(pq.uniquePeptideCountsPerExperiment(peptideIntensityListThreeFiles, experimentToIdxMap, 0.1), [1, 1, 1])

  def test_uniquePeptideCountsPerExperiment_missingValuesZeroes(self, peptideIntensityListMissingValuesZeroes, experimentToIdxMap):
    np.testing.assert_equal(pq.uniquePeptideCountsPerExperiment(peptideIntensityListMissingValuesZeroes, experimentToIdxMap, 0.1), [1, 1, 1])

  def test_uniquePeptideCountsPerExperiment_missingValuesNans(self, peptideIntensityListMissingValuesNans, experimentToIdxMap):
    np.testing.assert_equal(pq.uniquePeptideCountsPerExperiment(peptideIntensityListMissingValuesNans, experimentToIdxMap, 0.1), [1, 0, 0])

  def test_uniquePeptideCountsPerExperiment_duplicatePeptides(self, peptideIntensityListDuplicatePeptides, experimentToIdxMap):
    np.testing.assert_equal(pq.uniquePeptideCountsPerExperiment(peptideIntensityListDuplicatePeptides, experimentToIdxMap, 0.1), [1, 1, 1])

  def test_uniquePeptideCountsPerExperiment_2peptides(self, peptideIntensityListMultiplePeptides, experimentToIdxMap):
    np.testing.assert_equal(pq.uniquePeptideCountsPerExperiment(peptideIntensityListMultiplePeptides, experimentToIdxMap, 0.1), [2, 2, 2])

