from __future__ import annotations

from typing import List, Dict
import logging

from .. import helpers
from .base import ProteinGroupColumns

# for type hints only
from .precursor_quant import PrecursorQuant
from .. import results


logger = logging.getLogger(__name__)


class IdentificationTypeColumns(ProteinGroupColumns):
    def append_headers(
        self,
        proteinGroupResults: results.ProteinGroupResults,
        experiments: List[str],
    ) -> None:
        for experiment in experiments:
            proteinGroupResults.append_header("Identification type " + experiment)

    def append_columns(
        self,
        proteinGroupResults: results.ProteinGroupResults,
        experimentToIdxMap: Dict[str, int],
        postErrProbCutoff: float,
    ) -> None:
        logger.info("Doing quantification: adding identification types")
        for pgr in proteinGroupResults:
            idType = _identification_type_per_experiment(
                pgr.precursorQuants, experimentToIdxMap, postErrProbCutoff
            )
            pgr.extend(idType)


def _identification_type_per_experiment(
    peptideIntensityList: List[PrecursorQuant],
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
