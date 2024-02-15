import pytest
import numpy as np

from picked_group_fdr.columns.sum_and_ibaq import _get_intensities
from picked_group_fdr.writers.base import _retain_only_identified_precursors


class TestGetIntensities:
    def test_getIntensities_3files(
        self, peptideIntensityListThreeFiles, experimentToIdxMap
    ):
        np.testing.assert_almost_equal(
            _get_intensities(
                peptideIntensityListThreeFiles,
                experimentToIdxMap,
                0.1,
                numSilacChannels=0,
            ),
            [25.0, 10.0, 5.0],
        )

    def test_getIntensities_5files(
        self, peptideIntensityListFiveFiles, experimentToIdxMapFiveFiles
    ):
        np.testing.assert_almost_equal(
            _get_intensities(
                peptideIntensityListFiveFiles,
                experimentToIdxMapFiveFiles,
                0.1,
                numSilacChannels=0,
            ),
            [25.0, 10.0, 5.0, 0.0, 7.0],
        )

    def test_getIntensities_missingValuesZeroes(
        self, peptideIntensityListMissingValuesZeroes, experimentToIdxMap
    ):
        np.testing.assert_almost_equal(
            _get_intensities(
                peptideIntensityListMissingValuesZeroes,
                experimentToIdxMap,
                0.1,
                numSilacChannels=0,
            ),
            [25.0, 0.0, 0.0],
        )

    def test_getIntensities_missingValuesNans(
        self, peptideIntensityListMissingValuesNans, experimentToIdxMap
    ):
        np.testing.assert_almost_equal(
            _get_intensities(
                peptideIntensityListMissingValuesNans,
                experimentToIdxMap,
                0.1,
                numSilacChannels=0,
            ),
            [25.0, 0.0, 0.0],
        )

    def test_getIntensities_2peptides(
        self, peptideIntensityListMultiplePeptides, experimentToIdxMap
    ):
        np.testing.assert_almost_equal(
            _get_intensities(
                peptideIntensityListMultiplePeptides,
                experimentToIdxMap,
                0.1,
                numSilacChannels=0,
            ),
            [40.0, 20.0, 10.0],
        )

    def test_getIntensities_duplicatePeptides(
        self, peptideIntensityListDuplicatePeptides, experimentToIdxMap
    ):
        np.testing.assert_almost_equal(
            _get_intensities(
                peptideIntensityListDuplicatePeptides,
                experimentToIdxMap,
                0.1,
                numSilacChannels=0,
            ),
            [35.0, 10.0, 5.0],
        )

    def test_getIntensities_2fractions(
        self, peptideIntensityListMultipleFractions, experimentToIdxMap
    ):
        np.testing.assert_almost_equal(
            _get_intensities(
                peptideIntensityListMultipleFractions,
                experimentToIdxMap,
                0.1,
                numSilacChannels=0,
            ),
            [35.0, 10.0, 5.0],
        )

    def test_getIntensities_silac(self, peptideIntensityListSilac, experimentToIdxMap):
        np.testing.assert_almost_equal(
            _get_intensities(
                peptideIntensityListSilac, experimentToIdxMap, 0.1, numSilacChannels=3
            ),
            [35.0, 13.0, 17.0, 5.0, 10.0, 3.0, 4.0, 3.0, 5.0, 1.5, 2.5, 1.0],
        )

    def test_getIntensities_matchBetweenRuns(
        self, peptideIntensityListMatchBetweenRuns, experimentToIdxMap
    ):
        postErrProbCutoff = 0.1
        peptideIntensityList = _retain_only_identified_precursors(
            peptideIntensityListMatchBetweenRuns, postErrProbCutoff
        )
        np.testing.assert_almost_equal(
            _get_intensities(
                peptideIntensityList,
                experimentToIdxMap,
                postErrProbCutoff,
                numSilacChannels=0,
            ),
            [25.0, 10.0, 5.0],
        )

    def test_getIntensities_matchBetweenRunsWithUnidentified(
        self, peptideIntensityListMatchBetweenRunsWithUnidentified, experimentToIdxMap
    ):
        postErrProbCutoff = 0.1
        peptideIntensityList = _retain_only_identified_precursors(
            peptideIntensityListMatchBetweenRunsWithUnidentified, postErrProbCutoff
        )
        np.testing.assert_almost_equal(
            _get_intensities(
                peptideIntensityList,
                experimentToIdxMap,
                postErrProbCutoff,
                numSilacChannels=0,
            ),
            [0.0, 0.0, 0.0],
        )
