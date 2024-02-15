import pytest
import numpy as np

import picked_group_fdr.fdr as pq


def test_calcPostErrProbCutoff():
    postErrProbs = [0.12, 0.001, 0.02]
    psmQvalCutoff = 0.01
    assert pq.calc_post_err_prob_cutoff(postErrProbs, psmQvalCutoff) == 0.02


def test_calcPostErrProbCutoff_withNans():
    postErrProbs = [np.nan, 0.001, 0.02, np.inf]
    psmQvalCutoff = 0.01
    assert pq.calc_post_err_prob_cutoff(postErrProbs, psmQvalCutoff) == 0.02
