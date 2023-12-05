from __future__ import annotations

import logging

import numpy as np
import triqler.qvality

from .base import ProteinGroupColumns

# for type hints only
from .. import results

logger = logging.getLogger(__name__)


class ProteinProbabilityColumns(ProteinGroupColumns):
    def append_headers(
        self,
        protein_group_results: results.ProteinGroupResults,
    ) -> None:
        protein_group_results.append_header("Protein Probability")

    def append_columns(
        self,
        protein_group_results: results.ProteinGroupResults,
        post_err_prob_cutoff: float,
    ) -> None:
        """_summary_

        NOTE: this expects the protein group results to be sorted by descending score.
        """
        logger.info("Doing quantification: Computing protein probabilities")
        target_scores, decoy_scores = list(), list()
        for pgr in protein_group_results:
            if pgr.reverse == "+":
                decoy_scores.append(pgr.score)
            else:
                target_scores.append(pgr.score)

        target_scores = np.array(target_scores)
        decoy_scores = np.array(decoy_scores)

        target_scores, decoy_scores = transform_scores(target_scores, decoy_scores)

        _, protein_post_err_probs = triqler.qvality.getQvaluesFromScores(
            target_scores,
            decoy_scores,
            includePEPs=True,
            includeDecoys=True,
            tdcInput=True,
        )

        for pgr, protein_post_err_prob in zip(
            protein_group_results, protein_post_err_probs
        ):
            pgr.append("%.3f" % (1.0 - protein_post_err_prob,))


def transform_scores(target_scores, decoy_scores):
    """Apply logit or log transform to scores to reduce IRLS convergence problems.
    
    see BaseSpline::setData in
    # https://github.com/percolator/percolator/blob/master/src/BaseSpline.cpp
    """
    all_scores = np.concatenate((target_scores, decoy_scores))
    min_score = np.min(all_scores)
    max_score = np.max(all_scores)
    if np.min(all_scores) >= 0.0 and np.max(all_scores) <= 1.0:
        target_scores, decoy_scores = logit_transform(
            target_scores, decoy_scores, min_score, max_score
        )
    elif np.min(all_scores) >= 0.0:
        target_scores, decoy_scores = log_transform(
            target_scores, decoy_scores, min_score
        )
    return target_scores, decoy_scores


def logit_transform(target_scores, decoy_scores, min_score, max_score):
    delta_low = 0.0 if min_score > 0.0 else 1e-20
    delta_high = 0.0 if max_score < 1.0 else 1e-10

    target_scores = (target_scores + delta_low) * (1 - delta_low - delta_high)
    target_scores = np.log(target_scores / (1 - target_scores))

    decoy_scores = (decoy_scores + delta_low) * (1 - delta_low - delta_high)
    decoy_scores = np.log(decoy_scores / (1 - decoy_scores))
    return target_scores, decoy_scores


def log_transform(target_scores, decoy_scores, min_score):
    delta_low = 0.0 if min_score > 0.0 else 1e-20

    target_scores = target_scores + delta_low
    target_scores = np.log(target_scores)

    decoy_scores = decoy_scores + delta_low
    decoy_scores = np.log(decoy_scores)
    return target_scores, decoy_scores
