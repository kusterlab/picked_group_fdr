from __future__ import annotations

from typing import Dict, List
import logging

import numpy as np

from .. import helpers
from .base import ProteinGroupColumns

# for type hints only
from .. import precursor_quant
from .. import results


logger = logging.getLogger(__name__)


class TMTIntensityColumns(ProteinGroupColumns):
    def is_valid(self, protein_group_results: results.ProteinGroupResults):
        return protein_group_results.num_tmt_channels > 0
    
    def append_headers(self, protein_group_results: results.ProteinGroupResults):
        for experiment in protein_group_results.experiments:
            for i in range(1, protein_group_results.num_tmt_channels + 1):
                protein_group_results.append_header(
                    "Reporter intensity corrected " + str(i) + " " + experiment
                )

            for i in range(1, protein_group_results.num_tmt_channels + 1):
                protein_group_results.append_header(
                    "Reporter intensity " + str(i) + " " + experiment
                )

            for i in range(1, protein_group_results.num_tmt_channels + 1):
                protein_group_results.append_header(
                    "Reporter intensity count " + str(i) + " " + experiment
                )

    def append_columns(
        self,
        protein_group_results: results.ProteinGroupResults,
        post_err_prob_cutoff,
    ):
        logger.info("Doing quantification: TMT intensity")
        experiment_to_idx_map = protein_group_results.get_experiment_to_idx_map()
        for pgr in protein_group_results:
            intensities = _get_tmt_intensities(
                pgr.precursorQuants,
                experiment_to_idx_map,
                post_err_prob_cutoff,
                protein_group_results.num_tmt_channels,
            )
            pgr.extend(intensities)


def _get_tmt_intensities(
    precursor_list: List[precursor_quant.PrecursorQuant],
    experiment_to_idx_map: Dict[str, int],
    post_err_prob_cutoff: float,
    num_tmt_channels: int,
):
    intensities = [
        np.zeros(num_tmt_channels * 3) for i in range(len(experiment_to_idx_map))
    ]
    for precursor in precursor_list:
        if (
            helpers.is_mbr(precursor.post_err_prob)
            or precursor.post_err_prob <= post_err_prob_cutoff
        ):
            intensities[
                experiment_to_idx_map[precursor.experiment]
            ] += precursor.tmt_intensities
    return np.concatenate(intensities)
