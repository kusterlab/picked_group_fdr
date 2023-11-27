from typing import List, Dict
import logging

import numpy as np
import triqler.qvality

from .. import helpers
from .base import ProteinGroupColumns

# for type hints only
from .precursor_quant import PrecursorQuant
from ..results import ProteinGroupResults

logger = logging.getLogger(__name__)


class ProteinProbabilityColumns(ProteinGroupColumns):
    def append_headers(
        self,
        protein_group_results: ProteinGroupResults,
        experiments: List[str],
    ) -> None:
        protein_group_results.append_header("Protein Probability")

    def append_columns(
        self,
        protein_group_results: ProteinGroupResults,
        experiment_to_idx_map: Dict[str, int],
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
        _, protein_post_err_probs = triqler.qvality.getQvaluesFromScores(
            target_scores,
            decoy_scores,
            includePEPs=True,
            includeDecoys=True,
            tdcInput=True,
        )

        for pgr, protein_post_err_prob in zip(protein_group_results, protein_post_err_probs):
            pgr.append("%.3f" % (1.0 - protein_post_err_prob, ))