from typing import List, Dict
import logging

from .. import helpers
from .base import ProteinGroupColumns

# for type hints only
from .precursor_quant import PrecursorQuant
from ..results import ProteinGroupResults

logger = logging.getLogger(__name__)


class SpectralCountColumns(ProteinGroupColumns):
    def append_headers(
        self,
        protein_group_results: ProteinGroupResults,
        experiments: List[str],
    ) -> None:
        protein_group_results.append_header("Combined Spectral Count")
        for experiment in experiments:
            protein_group_results.append_header("Spectral count " + experiment)

    def append_columns(
        self,
        protein_group_results: ProteinGroupResults,
        experiment_to_idx_map: Dict[str, int],
        post_err_prob_cutoff: float,
    ) -> None:
        logger.info("Doing quantification: Spectral count")
        for pgr in protein_group_results:
            pepCounts = _spectral_counts_per_experiment(
                pgr.precursorQuants, experiment_to_idx_map, post_err_prob_cutoff
            )
            pgr.append(sum(pepCounts))
            pgr.extend(pepCounts)


def _spectral_counts_per_experiment(
    precursor_list: List[PrecursorQuant],
    experiment_to_idx_map: Dict[str, int],
    post_err_prob_cutoff: float,
):
    spectral_counts = [0] * len(experiment_to_idx_map)
    for precursor in precursor_list:
        if (
            not helpers.isMbr(precursor.post_err_prob)
            and precursor.post_err_prob <= post_err_prob_cutoff
        ):
            spectral_counts[experiment_to_idx_map[precursor.experiment]] += 1
    return spectral_counts