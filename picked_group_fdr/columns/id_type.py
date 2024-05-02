from __future__ import annotations

from typing import List, Dict
import logging

from .. import helpers
from .base import ProteinGroupColumns

# for type hints only
from .. import precursor_quant
from .. import results


logger = logging.getLogger(__name__)


class IdentificationTypeColumns(ProteinGroupColumns):
    def append_headers(
        self,
        protein_group_results: results.ProteinGroupResults,
    ) -> None:
        for experiment in protein_group_results.experiments:
            protein_group_results.append_header("Identification type " + experiment)

    def append_columns(
        self,
        protein_group_results: results.ProteinGroupResults,
        post_err_prob_cutoff: float,
    ) -> None:
        logger.info("Doing quantification: adding identification types")
        experiment_to_idx_map = protein_group_results.get_experiment_to_idx_map()
        for pgr in protein_group_results:
            id_type = _identification_type_per_experiment(
                pgr.precursorQuants, experiment_to_idx_map, post_err_prob_cutoff
            )
            pgr.extend(id_type)


def _identification_type_per_experiment(
    peptideIntensityList: List[precursor_quant.PrecursorQuant],
    experimentToIdxMap: Dict[str, int],
    postErrProbCutoff: float,
):
    idType = [""] * len(experimentToIdxMap)
    for precursor in peptideIntensityList:
        idx = experimentToIdxMap[precursor.experiment]
        if helpers.is_mbr(precursor.post_err_prob) and idType[idx] != "By MS/MS":
            idType[idx] = "By matching"
        elif precursor.post_err_prob <= postErrProbCutoff:
            idType[idx] = "By MS/MS"
    return idType
