from __future__ import annotations

from typing import List, Dict, Set
import logging

import numpy as np

from .. import helpers
from .base import ProteinGroupColumns

# for type hints only
from .. import precursor_quant
from .. import results

logger = logging.getLogger(__name__)


class SequenceCoverageColumns(ProteinGroupColumns):
    protein_sequences: Dict[str, str]

    def __init__(self, protein_sequences):
        self.protein_sequences = protein_sequences

    def append_headers(
        self,
        protein_group_results: results.ProteinGroupResults,
    ) -> None:
        # TODO: properly implement this, right now it is the same as Unique sequence coverage [%]
        protein_group_results.append_header("Sequence coverage [%]")
        # TODO: properly implement this, right now it is the same as Unique sequence coverage [%]
        protein_group_results.append_header("Unique + razor sequence coverage [%]")
        protein_group_results.append_header("Unique sequence coverage [%]")
        # TODO: add these columns
        # Mol. weight [kDa]
        # Sequence length
        # Sequence lengths

        for experiment in protein_group_results.experiments:
            protein_group_results.append_header("Sequence coverage [%] " + experiment)

    def append_columns(
        self,
        protein_group_results: results.ProteinGroupResults,
        post_err_prob_cutoff: float,
    ) -> None:
        logger.info("Doing quantification: sequence coverage")
        experiment_to_idx_map = protein_group_results.get_experiment_to_idx_map()
        for pgr in protein_group_results:
            sequence_coverages = self.get_sequence_coverages(
                pgr.precursorQuants,
                experiment_to_idx_map,
                post_err_prob_cutoff,
                pgr.proteinIds,
            )
            pgr.extend(sequence_coverages)

    def get_sequence_coverages(
        self,
        precursor_list: List[precursor_quant.PrecursorQuant],
        experiment_to_idx_map: Dict[str, int],
        post_err_prob_cutoff: float,
        protein_ids: str,
    ) -> List[str]:
        peptide_set_per_experiment = self.unique_peptides_per_experiment(
            precursor_list, experiment_to_idx_map, post_err_prob_cutoff
        )
        return self.calculate_sequence_coverages(
            peptide_set_per_experiment, protein_ids
        )

    def unique_peptides_per_experiment(
        self,
        precursor_list: List[precursor_quant.PrecursorQuant],
        experiment_to_idx_map: Dict[str, int],
        post_err_prob_cutoff: float,
    ) -> List[Set[str]]:
        uniquePeptides = [set() for _ in range(len(experiment_to_idx_map))]
        for precursor in precursor_list:
            if (
                helpers.is_mbr(precursor.post_err_prob)
                or precursor.post_err_prob <= post_err_prob_cutoff
            ):
                uniquePeptides[experiment_to_idx_map[precursor.experiment]].add(
                    helpers.remove_modifications(precursor.peptide)
                )
        return uniquePeptides

    def calculate_sequence_coverages(
        self, peptide_set_per_experiment: List[Set[str]], proteinIds: str
    ) -> List[float]:
        leading_protein = proteinIds.split(";")[0]
        protein_sequence = self.protein_sequences.get(leading_protein, "")
        coverage_total = np.zeros(len(protein_sequence))
        coverage_experiment_ratios = list()
        for peptides in peptide_set_per_experiment:
            if len(peptides) == 0:
                coverage_experiment_ratios.append(0.0)
                continue

            coverage_experiment = np.zeros(len(protein_sequence))
            for clean_peptide in peptides:
                pos = protein_sequence.find(clean_peptide)
                coverage_total[pos : pos + len(clean_peptide)] = 1
                coverage_experiment[pos : pos + len(clean_peptide)] = 1

            coverage_experiment_ratios.append(
                self.calculate_coverage_ratio(coverage_experiment)
            )

        coverageTotalRatio = self.calculate_coverage_ratio(coverage_total)
        allCoverageRatios = [
            coverageTotalRatio,
            coverageTotalRatio,
            coverageTotalRatio,
        ] + coverage_experiment_ratios
        return self.format_as_percentage(allCoverageRatios)

    @staticmethod
    def format_as_percentage(l: List[float]) -> List[str]:
        return ["%.1f" % (x * 100) for x in l]

    @staticmethod
    def calculate_coverage_ratio(coverage: List[int]) -> float:
        if len(coverage) > 0:
            return coverage.sum() / len(coverage)
        else:
            return 0.0
