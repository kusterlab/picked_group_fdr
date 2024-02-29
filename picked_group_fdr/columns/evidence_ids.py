from __future__ import annotations

from typing import List
import logging

from .. import helpers
from .base import ProteinGroupColumns

# for type hints only
from ..precursor_quant import PrecursorQuant
from .. import results


logger = logging.getLogger(__name__)


class EvidenceIdsColumns(ProteinGroupColumns):
    def append_headers(
        self,
        proteinGroupResults: results.ProteinGroupResults,
    ) -> None:
        proteinGroupResults.append_header("Evidence IDs")

    def append_columns(
        self,
        proteinGroupResults: results.ProteinGroupResults,
        postErrProbCutoff: float,
    ) -> None:
        logger.info("Doing quantification: adding evidence ids")
        for pgr in proteinGroupResults:
            evidenceIds = _collect_evidence_ids(
                pgr.precursorQuants, postErrProbCutoff
            )
            pgr.append(";".join(map(str, evidenceIds)))


def _collect_evidence_ids(
    precursor_list: List[PrecursorQuant], post_err_prob_cutoff: float
):
    evidence_ids = list()
    for precursor in precursor_list:
        if (
            helpers.is_mbr(precursor.post_err_prob)
            or precursor.post_err_prob <= post_err_prob_cutoff
        ):
            evidence_ids.append(precursor.evidence_id)
    return sorted(evidence_ids)
