import numpy as np

from picked_group_fdr.columns.peptide_count import _unique_peptide_counts_per_experiment


class TestUniquePeptideCountsPerExperiment:
    def test_uniquePeptideCountsPerExperiment_3files(
        self, peptideIntensityListThreeFiles, experimentToIdxMap
    ):
        np.testing.assert_equal(
            _unique_peptide_counts_per_experiment(
                peptideIntensityListThreeFiles, experimentToIdxMap, 0.1
            ),
            [1, 1, 1],
        )

    def test_uniquePeptideCountsPerExperiment_missingValuesZeroes(
        self, peptideIntensityListMissingValuesZeroes, experimentToIdxMap
    ):
        np.testing.assert_equal(
            _unique_peptide_counts_per_experiment(
                peptideIntensityListMissingValuesZeroes, experimentToIdxMap, 0.1
            ),
            [1, 1, 1],
        )

    def test_uniquePeptideCountsPerExperiment_missingValuesNans(
        self, peptideIntensityListMissingValuesNans, experimentToIdxMap
    ):
        np.testing.assert_equal(
            _unique_peptide_counts_per_experiment(
                peptideIntensityListMissingValuesNans, experimentToIdxMap, 0.1
            ),
            [1, 0, 0],
        )

    def test_uniquePeptideCountsPerExperiment_duplicatePeptides(
        self, peptideIntensityListDuplicatePeptides, experimentToIdxMap
    ):
        np.testing.assert_equal(
            _unique_peptide_counts_per_experiment(
                peptideIntensityListDuplicatePeptides, experimentToIdxMap, 0.1
            ),
            [1, 1, 1],
        )

    def test_uniquePeptideCountsPerExperiment_2peptides(
        self, peptideIntensityListMultiplePeptides, experimentToIdxMap
    ):
        np.testing.assert_equal(
            _unique_peptide_counts_per_experiment(
                peptideIntensityListMultiplePeptides, experimentToIdxMap, 0.1
            ),
            [2, 2, 2],
        )
