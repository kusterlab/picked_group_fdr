import pytest

import numpy as np

import picked_group_fdr.fdr as fdr


def test_count_below_threshold():
    assert (
        fdr.count_below_threshold([0.001, 0.002, 0.003, 0.004, 0.1, 0.2, 0.3], 0.01)
        == 4
    )


def test_count_below_threshold_with_decoy_labels():
    assert (
        fdr.count_below_threshold(
            [0.001, 0.002, 0.003, 0.004, 0.1, 0.2, 0.3],
            0.01,
            [False, False, True, False, False, False, True],
        )
        == 3
    )


def test_fdr_to_qvals():
    np.testing.assert_almost_equal(
        fdr.fdrs_to_qvals([0.1, 0.01, 0.001, 0.002]), [0.001, 0.001, 0.001, 0.002]
    )
