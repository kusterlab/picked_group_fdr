from __future__ import annotations

from typing import List, Dict
import logging

from .. import helpers
from .base import ProteinGroupColumns

# for type hints only
from .. import precursor_quant
from .. import results

logger = logging.getLogger(__name__)


class SpectralCountColumns(ProteinGroupColumns):
    def append_headers(
        self,
        protein_group_results: results.ProteinGroupResults,
    ) -> None:
        protein_group_results.append_header("Combined Spectral Count")
        for experiment in protein_group_results.experiments:
            protein_group_results.append_header("Spectral count " + experiment)

    def append_columns(
        self,
        protein_group_results: results.ProteinGroupResults,
        post_err_prob_cutoff: float,
    ) -> None:
        experiment_to_idx_map = protein_group_results.get_experiment_to_idx_map()
        
        logger.info("Doing quantification: Spectral count")
        for pgr in protein_group_results:
            pepCounts = _spectral_counts_per_experiment(
                pgr.precursorQuants, experiment_to_idx_map, post_err_prob_cutoff
            )
            pgr.append(sum(pepCounts))
            pgr.extend(pepCounts)


def _spectral_counts_per_experiment(
    precursor_list: List[precursor_quant.PrecursorQuant],
    experiment_to_idx_map: Dict[str, int],
    post_err_prob_cutoff: float,
):
    spectral_counts = [0] * len(experiment_to_idx_map)
    for precursor in precursor_list:
        if (
            not helpers.is_mbr(precursor.post_err_prob)
            and precursor.post_err_prob <= post_err_prob_cutoff
        ):
            spectral_counts[experiment_to_idx_map[precursor.experiment]] += 1
    return spectral_counts