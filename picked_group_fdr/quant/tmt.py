from __future__ import annotations

from typing import List
import logging

import numpy as np

from .. import helpers
from .base import ProteinGroupColumns

# for type hints only
from .precursor_quant import PrecursorQuant
from .. import results


logger = logging.getLogger(__name__)


class TMTIntensityColumns(ProteinGroupColumns):
    numTmtChannels: int

    def __init__(self, num_tmt_channels: int):
        self.numTmtChannels = num_tmt_channels

    def append_headers(
        self, proteinGroupResults: results.ProteinGroupResults, experiments: List[str]
    ):
        for experiment in experiments:
            for i in range(1, self.numTmtChannels + 1):
                proteinGroupResults.append_header(
                    "Reporter intensity corrected " + str(i) + " " + experiment
                )

            for i in range(1, self.numTmtChannels + 1):
                proteinGroupResults.append_header(
                    "Reporter intensity " + str(i) + " " + experiment
                )

            for i in range(1, self.numTmtChannels + 1):
                proteinGroupResults.append_header(
                    "Reporter intensity count " + str(i) + " " + experiment
                )

    def append_columns(
        self,
        proteinGroupResults: results.ProteinGroupResults,
        experimentToIdxMap,
        postErrProbCutoff,
    ):
        logger.info("Doing quantification: TMT intensity")
        for pgr in proteinGroupResults:
            intensities = self._getTmtIntensities(
                pgr.precursorQuants, experimentToIdxMap, postErrProbCutoff
            )
            pgr.extend(intensities)

    def _getTmtIntensities(
        self,
        peptideIntensityList: List[PrecursorQuant],
        experimentToIdxMap,
        postErrProbCutoff,
    ):
        intensities = [
            np.zeros(self.numTmtChannels * 3) for i in range(len(experimentToIdxMap))
        ]
        for precursor in peptideIntensityList:
            if (
                helpers.is_mbr(precursor.post_err_prob)
                or precursor.post_err_prob <= postErrProbCutoff
            ):
                intensities[
                    experimentToIdxMap[precursor.experiment]
                ] += precursor.tmt_intensities
        return np.concatenate(intensities)
