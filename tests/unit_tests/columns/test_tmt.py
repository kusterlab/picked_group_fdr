import numpy as np

from picked_group_fdr.columns import tmt


def test_getTmtIntensities(peptideIntensityListTMT, experimentToIdxMap):
    num_tmt_channels = 1
    np.testing.assert_almost_equal(
        tmt._get_tmt_intensities(
            peptideIntensityListTMT, experimentToIdxMap, 0.1, num_tmt_channels
        ),
        [3.0, 7.0, 4.0, 3.0, 4.0, 5.0, 5.0, 6.0, 7.0],
    )
