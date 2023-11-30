import pytest
import numpy as np

from picked_group_fdr.quant.tmt import TMTIntensityColumns


def test_getTmtIntensities(peptideIntensityListTMT, experimentToIdxMap):
  pq = TMTIntensityColumns(num_tmt_channels=1)
  np.testing.assert_almost_equal(pq._getTmtIntensities(peptideIntensityListTMT, experimentToIdxMap, 0.1), [3.0, 7.0, 4.0, 3.0, 4.0, 5.0, 5.0, 6.0, 7.0])
