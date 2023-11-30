from typing import List, Dict
import logging

from .. import helpers
from .base import ProteinGroupColumns

# for type hints only
from .precursor_quant import PrecursorQuant
from ..results import ProteinGroupResults

logger = logging.getLogger(__name__)


class UniquePeptideCountColumns(ProteinGroupColumns):
    def append_headers(
        self,
        protein_group_results: ProteinGroupResults,
        experiments: List[str],
    ) -> None:
        protein_group_results.append_header("Combined Total Peptides")
        for experiment in experiments:
            protein_group_results.append_header("Unique peptides " + experiment)

    def append_columns(
        self,
        protein_group_results: ProteinGroupResults,
        experiment_to_idx_map: Dict[str, int],
        post_err_prob_cutoff: float,
    ) -> None:
        logger.info("Doing quantification: Count unique peptides")
        for pgr in protein_group_results:
            pepCount = _unique_peptide_counts_combined(
                pgr.precursorQuants, post_err_prob_cutoff
            )
            pgr.append(pepCount)

            pepCounts = _unique_peptide_counts_per_experiment(
                pgr.precursorQuants, experiment_to_idx_map, post_err_prob_cutoff
            )
            pgr.extend(pepCounts)


def _unique_peptide_counts_per_experiment(
    precursor_list: List[PrecursorQuant],
    experiment_to_idx_map: Dict[str, int],
    post_err_prob_cutoff: float,
):
    uniquePeptides = [set() for _ in range(len(experiment_to_idx_map))]
    for precursor in precursor_list:
        if (
            helpers.is_mbr(precursor.post_err_prob)
            or precursor.post_err_prob <= post_err_prob_cutoff
        ):
            uniquePeptides[experiment_to_idx_map[precursor.experiment]].add(
                precursor.peptide
            )
    return list(map(len, uniquePeptides))


def _unique_peptide_counts_combined(
    precursor_list: List[PrecursorQuant],
    post_err_prob_cutoff: float,
):
    uniquePeptides = set()
    for precursor in precursor_list:
        if (
            helpers.is_mbr(precursor.post_err_prob)
            or precursor.post_err_prob <= post_err_prob_cutoff
        ):
            uniquePeptides.add(precursor.peptide)
    return len(uniquePeptides)
