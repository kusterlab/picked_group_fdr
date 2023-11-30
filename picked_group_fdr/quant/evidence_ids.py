from typing import List, Dict
import logging

from .. import helpers
from .base import ProteinGroupColumns

# for type hints only
from .precursor_quant import PrecursorQuant
from ..results import ProteinGroupResults


logger = logging.getLogger(__name__)


class EvidenceIdsColumns(ProteinGroupColumns):
    def append_headers(
        self,
        proteinGroupResults: ProteinGroupResults,
        experiments: List[str],
    ) -> None:
        proteinGroupResults.append_header("Evidence IDs")

    def append_columns(
        self,
        proteinGroupResults: ProteinGroupResults,
        experimentToIdxMap: Dict[str, int],
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
