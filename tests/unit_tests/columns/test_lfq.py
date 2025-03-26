import pytest
import numpy as np

from picked_group_fdr.columns.lfq import (
    _getLFQIntensities,
    _getMaxRatio,
    _getPeptideIntensities,
    _getLogMedianPeptideRatios,
)


class TestGetLFQIntensities:
    def test_getLFQIntensities_3files(
        self, peptideIntensityListThreeFiles, experimentToIdxMap
    ):
        np.testing.assert_almost_equal(
            _getLFQIntensities(
                peptideIntensityListThreeFiles, experimentToIdxMap, 0.1, 1
            ),
            [25.0, 10.0, 5.0],
        )

    def test_getLFQIntensities_5files(
        self, peptideIntensityListFiveFiles, experimentToIdxMapFiveFiles
    ):
        np.testing.assert_almost_equal(
            _getLFQIntensities(
                peptideIntensityListFiveFiles, experimentToIdxMapFiveFiles, 0.1, 1
            ),
            [25.0, 10.0, 5.0, 0.0, 7.0],
        )

    def test_getLFQIntensities_missingValuesZeroes(
        self, peptideIntensityListMissingValuesZeroes, experimentToIdxMap
    ):
        np.testing.assert_almost_equal(
            _getLFQIntensities(
                peptideIntensityListMissingValuesZeroes, experimentToIdxMap, 0.1, 1
            ),
            [0.0, 0.0, 0.0],
        )

    def test_getLFQIntensities_missingValuesNans(
        self, peptideIntensityListMissingValuesNans, experimentToIdxMap
    ):
        np.testing.assert_almost_equal(
            _getLFQIntensities(
                peptideIntensityListMissingValuesNans, experimentToIdxMap, 0.1, 1
            ),
            [0.0, 0.0, 0.0],
        )

    def test_getLFQIntensities_missingValuesNans_minPeptideRatiosLFQ2(
        self, peptideIntensityListMissingValuesNans, experimentToIdxMap
    ):
        np.testing.assert_almost_equal(
            _getLFQIntensities(
                peptideIntensityListMissingValuesNans, experimentToIdxMap, 0.1, 2
            ),
            [0.0, 0.0, 0.0],
        )

    def test_getLFQIntensities_2peptides(
        self, peptideIntensityListMultiplePeptides, experimentToIdxMap
    ):
        np.testing.assert_almost_equal(
            _getLFQIntensities(
                peptideIntensityListMultiplePeptides, experimentToIdxMap, 0.1, 2
            ),
            [40.0, 20.0, 10.0],
        )

    def test_getLFQIntensities_2peptidesJagged(
        self, peptideIntensityListMultiplePeptidesJagged, experimentToIdxMap
    ):
        np.testing.assert_almost_equal(
            _getLFQIntensities(
                peptideIntensityListMultiplePeptidesJagged, experimentToIdxMap, 0.1, 2
            ),
            [50.0, 20.0, 0.0],
        )

    def test_getLFQIntensities_tooFewPeptideRatios(
        self, peptideIntensityListMultiplePeptides, experimentToIdxMap
    ):
        np.testing.assert_almost_equal(
            _getLFQIntensities(
                peptideIntensityListMultiplePeptides, experimentToIdxMap, 0.1, 3
            ),
            [0.0, 0.0, 0.0],
        )

    def test_getLFQIntensities_duplicatePeptides(
        self, peptideIntensityListDuplicatePeptides, experimentToIdxMap
    ):
        np.testing.assert_almost_equal(
            _getLFQIntensities(
                peptideIntensityListDuplicatePeptides, experimentToIdxMap, 0.1, 1
            ),
            [25.0, 10.0, 5.0],
        )

    def test_getLFQIntensities_multipleCharges(
        self, peptideIntensityListMultipleCharges, experimentToIdxMap
    ):
        """MaxLFQ considers different charge states of the same peptide as two different peptides for the sake of minPeptideRatiosLFQ"""
        np.testing.assert_almost_equal(
            _getLFQIntensities(
                peptideIntensityListMultipleCharges, experimentToIdxMap, 0.1, 2
            ),
            [40.0, 20.0, 10.0],
        )

    def test_getLFQIntensities_2fractions(
        self, peptideIntensityListMultipleFractions, experimentToIdxMap
    ):
        np.testing.assert_almost_equal(
            _getLFQIntensities(
                peptideIntensityListMultipleFractions, experimentToIdxMap, 0.1, 1
            ),
            [35.0, 10.0, 5.0],
        )

    def test_getLFQIntensities_silac(
        self, peptideIntensityListSilac, experimentToIdxMap
    ):
        np.testing.assert_almost_equal(
            _getLFQIntensities(
                peptideIntensityListSilac, experimentToIdxMap, 0.1, 1, False, 3
            ),
            [13.0, 17.0, 5.0, 3.0, 4.0, 3.0, 1.5, 2.5, 1.0],
        )


class TestStabilizeLargeRatios:
    def test_getMaxRatio(self):
        assert _getMaxRatio(5, 2) == 5 / 2

    def test_getMaxRatio_reverse(self):
        assert _getMaxRatio(2, 5) == 5 / 2

    def test_getMaxRatio_both_zero(self):
        with pytest.raises(AssertionError):
            _getMaxRatio(0, 0)

    def test_getLFQIntensities_3files(
        self, peptideIntensityListThreeFiles, experimentToIdxMap
    ):
        np.testing.assert_almost_equal(
            _getLFQIntensities(
                peptideIntensityListThreeFiles, experimentToIdxMap, 0.1, 1
            ),
            [25.0, 10.0, 5.0],
        )


class TestGetPeptideIntensities:
    def test_getPeptideIntensities(
        self, peptideIntensityListThreeFiles, experimentToIdxMap
    ):
        postErrProbCutoff = 0.1
        numSilacChannels = 0
        numExperiments = 3
        np.testing.assert_almost_equal(
            _getPeptideIntensities(
                peptideIntensityListThreeFiles,
                experimentToIdxMap,
                postErrProbCutoff,
                numSilacChannels,
                numExperiments,
            )[0][("APEPTIDE", 2)],
            [25.0, 10.0, 5.0],
        )


class TestCalculateRatios:
    def test_getLogMedianPeptideRatios(self, peptideIntensities):
        assert _getLogMedianPeptideRatios(peptideIntensities, 1) == {
            (0, 1): np.log(2.5),
            (0, 2): np.log(5.0),
            (1, 2): np.log(2.0),
        }


@pytest.fixture
def peptideIntensities():
    return {("APEPTIDE", 2): [25.0, 10.0, 5.0]}
