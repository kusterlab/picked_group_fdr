from __future__ import annotations

from typing import List
import logging

from .. import helpers
from .base import ProteinGroupColumns

# for type hints only
from .. import precursor_quant
from .. import results

logger = logging.getLogger(__name__)


class TopPeptideProbabilityColumns(ProteinGroupColumns):
    def append_headers(
        self,
        protein_group_results: results.ProteinGroupResults,
    ) -> None:    
        protein_group_results.append_header("Top Peptide Probability")

    def append_columns(
        self,
        protein_group_results: results.ProteinGroupResults,
        post_err_prob_cutoff: float,
    ) -> None:
        logger.info("Doing quantification: Top Peptide Probability")
        for pgr in protein_group_results:
            top_peptide_probability = _top_peptide_probability(
                pgr.precursorQuants
            )
            pgr.append("%.3f" % (top_peptide_probability, ))


def _top_peptide_probability(
    precursor_list: List[precursor_quant.PrecursorQuant],
) -> float:
    top_peptide_probability = 0.0
    for precursor in precursor_list:
        if (
            not helpers.is_mbr(precursor.post_err_prob) and
            1.0 - precursor.post_err_prob > top_peptide_probability
        ):
            top_peptide_probability = 1.0 - precursor.post_err_prob
    return top_peptide_probability
