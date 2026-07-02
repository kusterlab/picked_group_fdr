import pytest
import numpy as np

from picked_group_fdr.columns.lfq import (
    _get_lfq_intensities,
    _get_max_ratio,
    _get_peptide_intensities,
    _get_log_median_peptide_ratios,
    _apply_large_ratio_stabilization,
)


class TestGetLFQIntensities:
    def test_get_lfq_intensities_3files(
        self, peptideIntensityListThreeFiles, experimentToIdxMap
    ):
        np.testing.assert_almost_equal(
            _get_lfq_intensities(
                peptideIntensityListThreeFiles, experimentToIdxMap, 0.1, 1
            ),
            [25.0, 10.0, 5.0],
        )

    def test_get_lfq_intensities_5files(
        self, peptideIntensityListFiveFiles, experimentToIdxMapFiveFiles
    ):
        np.testing.assert_almost_equal(
            _get_lfq_intensities(
                peptideIntensityListFiveFiles, experimentToIdxMapFiveFiles, 0.1, 1
            ),
            [25.0, 10.0, 5.0, 0.0, 7.0],
        )

    def test_get_lfq_intensities_missingValuesZeroes_fallBackSingleSampleProtein(
        self, peptideIntensityListMissingValuesZeroes, experimentToIdxMap
    ):
        np.testing.assert_almost_equal(
            _get_lfq_intensities(
                peptideIntensityListMissingValuesZeroes, experimentToIdxMap, 0.1, 1
            ),
            [25.0, 0.0, 0.0],
        )

    def test_get_lfq_intensities_missingValuesNans_fallBackSingleSampleProtein(
        self, peptideIntensityListMissingValuesNans, experimentToIdxMap
    ):
        np.testing.assert_almost_equal(
            _get_lfq_intensities(
                peptideIntensityListMissingValuesNans, experimentToIdxMap, 0.1, 1
            ),
            [25.0, 0.0, 0.0],
        )

    def test_get_lfq_intensities_missingValuesNans_minPeptideRatiosLFQ2(
        self, peptideIntensityListMissingValuesNans, experimentToIdxMap
    ):
        np.testing.assert_almost_equal(
            _get_lfq_intensities(
                peptideIntensityListMissingValuesNans, experimentToIdxMap, 0.1, 2
            ),
            [0.0, 0.0, 0.0],
        )

    def test_get_lfq_intensities_2peptides(
        self, peptideIntensityListMultiplePeptides, experimentToIdxMap
    ):
        np.testing.assert_almost_equal(
            _get_lfq_intensities(
                peptideIntensityListMultiplePeptides, experimentToIdxMap, 0.1, 2
            ),
            [40.0, 20.0, 10.0],
        )

    def test_get_lfq_intensities_2peptidesJagged(
        self, peptideIntensityListMultiplePeptidesJagged, experimentToIdxMap
    ):
        np.testing.assert_almost_equal(
            _get_lfq_intensities(
                peptideIntensityListMultiplePeptidesJagged, experimentToIdxMap, 0.1, 2
            ),
            [50.0, 20.0, 0.0],
        )

    def test_get_lfq_intensities_tooFewPeptideRatios(
        self, peptideIntensityListMultiplePeptides, experimentToIdxMap
    ):
        np.testing.assert_almost_equal(
            _get_lfq_intensities(
                peptideIntensityListMultiplePeptides, experimentToIdxMap, 0.1, 3
            ),
            [0.0, 0.0, 0.0],
        )

    def test_get_lfq_intensities_duplicatePeptides(
        self, peptideIntensityListDuplicatePeptides, experimentToIdxMap
    ):
        np.testing.assert_almost_equal(
            _get_lfq_intensities(
                peptideIntensityListDuplicatePeptides, experimentToIdxMap, 0.1, 1
            ),
            [25.0, 10.0, 5.0],
        )

    def test_get_lfq_intensities_multipleCharges(
        self, peptideIntensityListMultipleCharges, experimentToIdxMap
    ):
        """MaxLFQ considers different charge states of the same peptide as two different peptides for the sake of minPeptideRatiosLFQ"""
        np.testing.assert_almost_equal(
            _get_lfq_intensities(
                peptideIntensityListMultipleCharges, experimentToIdxMap, 0.1, 2
            ),
            [40.0, 20.0, 10.0],
        )

    def test_get_lfq_intensities_2fractions(
        self, peptideIntensityListMultipleFractions, experimentToIdxMap
    ):
        np.testing.assert_almost_equal(
            _get_lfq_intensities(
                peptideIntensityListMultipleFractions, experimentToIdxMap, 0.1, 1
            ),
            [35.0, 10.0, 5.0],
        )

    def test_get_lfq_intensities_silac(
        self, peptideIntensityListSilac, experimentToIdxMap
    ):
        np.testing.assert_almost_equal(
            _get_lfq_intensities(
                precursor_list=peptideIntensityListSilac,
                experiment_to_idx_map=experimentToIdxMap,
                post_err_prob_cutoff=0.1,
                min_peptide_ratios_lfq=1,
                stabilize_large_ratios_lfq=False,
                fast_lfq_graph=None,
                num_silac_channels=3,
            ),
            [13.0, 17.0, 5.0, 3.0, 4.0, 3.0, 1.5, 2.5, 1.0],
        )


class TestStabilizeLargeRatios:
    def test_get_max_ratio(self):
        assert _get_max_ratio(5, 2) == 5 / 2

    def test_get_max_ratio_reverse(self):
        assert _get_max_ratio(2, 5) == 5 / 2

    def test_get_max_ratio_both_zero(self):
        with pytest.raises(AssertionError):
            _get_max_ratio(0, 0)

    def test_apply_large_ratio_stabilization_largeRatio_replacesWithIntensityRatio(
        self, peptideIntensityListLargeRatioStabilizationHighRatio, experimentToIdxMap
    ):
        """Test that when peptide_count_ratio > 5, the log ratio is replaced with log(si1/si2)"""
        log_median_peptide_ratios = {(0, 1): np.log(2.0)}  # original ratio

        result = _apply_large_ratio_stabilization(
            log_median_peptide_ratios=log_median_peptide_ratios,
            precursor_list=peptideIntensityListLargeRatioStabilizationHighRatio,
            experiment_to_idx_map=experimentToIdxMap,
            post_err_prob_cutoff=0.01,
            num_silac_channels=0,
        )
        assert isinstance(result, dict)

        # file1 has 1 peptide, file2 has 10 peptides, ratio = 10 > 5
        # intensities are equal (both 100), so log ratio should be log(100/100) = 0
        expected = np.log(100.0 / 100.0)
        np.testing.assert_almost_equal(result[(0, 1)], expected)

    def test_apply_large_ratio_stabilization_mediumRatio_usesWeightedAverage(
        self, peptideIntensityListLargeRatioStabilizationMediumRatio, experimentToIdxMap
    ):
        """Test that when 2.5 < peptide_count_ratio <= 5, weighted average is used"""
        original_ratio = np.log(2.0)
        log_median_peptide_ratios = {(0, 1): original_ratio}

        result = _apply_large_ratio_stabilization(
            log_median_peptide_ratios=log_median_peptide_ratios,
            precursor_list=peptideIntensityListLargeRatioStabilizationMediumRatio,
            experiment_to_idx_map=experimentToIdxMap,
            post_err_prob_cutoff=0.01,
            num_silac_channels=0,
        )

        # file1 has 2 peptides, file2 has 8 peptides, ratio = 4.0
        # w = (4.0 - 2.5) / 2.5 = 0.6
        # intensities are equal (both 100), so intensity_ratio = 0
        # result = 0.6 * 0 + 0.4 * log(2.0) = 0.4 * log(2.0)
        w = (4.0 - 2.5) / 2.5
        expected = w * np.log(100.0 / 100.0) + (1 - w) * original_ratio
        np.testing.assert_almost_equal(result[(0, 1)], expected)

    def test_apply_large_ratio_stabilization_smallRatio_unchanged(
        self, peptideIntensityListLargeRatioStabilizationSmallRatio, experimentToIdxMap
    ):
        """Test that when peptide_count_ratio <= 2.5, the ratio remains unchanged"""
        original_ratio = np.log(3.0)
        log_median_peptide_ratios = {(0, 1): original_ratio}

        result = _apply_large_ratio_stabilization(
            log_median_peptide_ratios=log_median_peptide_ratios,
            precursor_list=peptideIntensityListLargeRatioStabilizationSmallRatio,
            experiment_to_idx_map=experimentToIdxMap,
            post_err_prob_cutoff=0.01,
            num_silac_channels=0,
        )

        # file1 has 4 peptides, file2 has 8 peptides, ratio = 2.0 < 2.5
        # should remain unchanged
        np.testing.assert_almost_equal(result[(0, 1)], original_ratio)

    def test_apply_large_ratio_stabilization_zeroPeptideCount_skipped(
        self, peptideIntensityListLargeRatioStabilizationZeroCount, experimentToIdxMap
    ):
        """Test that pairs with zero peptide counts are skipped"""
        log_median_peptide_ratios = {(0, 1): np.log(2.0), (0, 2): np.log(3.0)}

        result = _apply_large_ratio_stabilization(
            log_median_peptide_ratios=log_median_peptide_ratios,
            precursor_list=peptideIntensityListLargeRatioStabilizationZeroCount,
            experiment_to_idx_map=experimentToIdxMap,
            post_err_prob_cutoff=0.01,
            num_silac_channels=0,
        )

        # file2 has 0 peptides, so (0, 1) pair should be skipped and unchanged
        np.testing.assert_almost_equal(result[(0, 1)], np.log(2.0))
        # (0, 2) pair has both peptides, ratio = 1.0, should not change (ratio <= 2.5)
        np.testing.assert_almost_equal(result[(0, 2)], np.log(3.0))

    def test_apply_large_ratio_stabilization_missingPair_skipped(
        self, peptideIntensityListThreeFiles, experimentToIdxMap
    ):
        """Test that pairs not in log_median_peptide_ratios dict are skipped"""
        log_median_peptide_ratios = {(0, 1): np.log(2.0)}  # only (0,1) present

        result = _apply_large_ratio_stabilization(
            log_median_peptide_ratios=log_median_peptide_ratios,
            precursor_list=peptideIntensityListThreeFiles,
            experiment_to_idx_map=experimentToIdxMap,
            post_err_prob_cutoff=0.01,
            num_silac_channels=0,
        )

        # Only (0, 1) should be present in result
        assert len(result) == 1
        assert (0, 1) in result
        np.testing.assert_almost_equal(result[(0, 1)], np.log(2.0))


class TestGetPeptideIntensities:
    def test_get_peptide_intensities(
        self, peptideIntensityListThreeFiles, experimentToIdxMap
    ):
        postErrProbCutoff = 0.1
        numSilacChannels = 0
        numExperiments = 3
        np.testing.assert_almost_equal(
            _get_peptide_intensities(
                peptideIntensityListThreeFiles,
                experimentToIdxMap,
                postErrProbCutoff,
                numSilacChannels,
                numExperiments,
            )[0][("APEPTIDE", 2)],
            [25.0, 10.0, 5.0],
        )


class TestCalculateRatios:
    def test_get_log_median_peptide_ratios(self, peptideIntensities):
        assert _get_log_median_peptide_ratios(peptideIntensities, 1) == {
            (0, 1): np.log(2.5),
            (0, 2): np.log(5.0),
            (1, 2): np.log(2.0),
        }


@pytest.fixture
def peptideIntensities():
    return {("APEPTIDE", 2): [25.0, 10.0, 5.0]}
