from __future__ import annotations

from typing import List, Dict
import logging

from .. import helpers
from .base import ProteinGroupColumns

# for type hints only
from .. import precursor_quant
from .. import results

logger = logging.getLogger(__name__)


class UniquePeptideCountColumns(ProteinGroupColumns):
    """Counts unique peptides per experiment and in total for each protein.

    Modified peptides are considered unique peptides themselves, e.g. 
    AP(ox)EPTIDE and APEPTIDE would count as 2 unique peptides.
    """    
    def append_headers(
        self,
        protein_group_results: results.ProteinGroupResults,
    ) -> None:
        protein_group_results.append_header("Combined Total Peptides")
        for experiment in protein_group_results.experiments:
            protein_group_results.append_header("Unique peptides " + experiment)

    def append_columns(
        self,
        protein_group_results: results.ProteinGroupResults,
        post_err_prob_cutoff: float,
    ) -> None:
        experiment_to_idx_map = protein_group_results.get_experiment_to_idx_map()
        
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
    precursor_list: List[precursor_quant.PrecursorQuant],
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
    precursor_list: List[precursor_quant.PrecursorQuant],
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
