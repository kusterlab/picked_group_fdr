from __future__ import annotations

import collections
from typing import Dict, List

from .base import ProteinGroupsWriter
from .. import columns
from ..protein_groups import ProteinGroups

from .fragpipe import fragpipe_format_extra_columns

# for type hints only
from .. import results
from .. import protein_annotation as pa


FRAGPIPE_COMBINED_PROTEIN_OUTPUT_DICT = {
    "Protein": "Protein",
    "Protein group": "Protein IDs",
    "Protein ID": "Protein ID",
    "Entry Name": "Entry Name",
    "Gene": "Gene",
    "Protein Length": "Length",
    "Organism": "Organism",
    "Protein Existence": "Protein Existence",
    "Description": "Protein Description",
    "Protein Probability": "Protein Probability",
    "Top Peptide Probability": "Top Peptide Probability",
    "Combined Total Peptides": "Combined Total Peptides",
    "Combined Spectral Count": "Combined Spectral Count",
    "Combined Unique Spectral Count": "Combined Spectral Count",
    "Combined Total Spectral Count": "Combined Spectral Count",
    # the following columns are not part of Fragpipe's format but were added upon popular demand
    "Majority protein IDs": "Majority protein IDs",
    "Q-value": "Q-value",
    "Score": "Score",
    "Reverse": "Reverse",
    "Potential contaminant": "Potential contaminant",
}


class FragPipeCombinedProteinWriter(ProteinGroupsWriter):
    def __init__(
        self,
        protein_groups: ProteinGroups,
        protein_annotations: Dict[str, pa.ProteinAnnotation],
        skip_lfq: bool = False,
        min_peptide_ratios_lfq: int = 1,
        stabilize_large_ratios_lfq: bool = True,
        num_threads: int = 1,
        protein_group_fdr_threshold: float = 0.01,
    ) -> None:
        self.protein_groups = protein_groups
        self.protein_annotations = protein_annotations
        self.skip_lfq = skip_lfq
        self.min_peptide_ratios_lfq = min_peptide_ratios_lfq
        self.stabilize_large_ratios_lfq = stabilize_large_ratios_lfq
        self.num_threads = num_threads

        # only used for reporting number of protein groups at the given threshold
        self.protein_group_fdr_threshold = protein_group_fdr_threshold

    def get_header_dict(
        self, protein_group_results: results.ProteinGroupResults
    ) -> Dict[str, str]:
        """Adds experiment specific headers.

        - <Experiment> Spectral Count
        - <Experiment> Unique Spectral Count
        - <Experiment> Total Spectral Count
        - <Experiment> Intensity
        - <Experiment> MaxLFQ Intensity

        - <Experiment> iBAQ Intensity
        - <Experiment> Unique Peptide Count
        """
        experiments = protein_group_results.experiments
        header_dict = FRAGPIPE_COMBINED_PROTEIN_OUTPUT_DICT.copy()
        for experiment in experiments:
            header_dict[f"{experiment} Spectra Count"] = f"Spectral count {experiment}"

        for experiment in experiments:
            header_dict[f"{experiment} Unique Spectra Count"] = (
                f"Spectral count {experiment}"
            )

        for experiment in experiments:
            header_dict[f"{experiment} Total Spectra Count"] = (
                f"Spectral count {experiment}"
            )

        for experiment in experiments:
            header_dict[f"{experiment} Intensity"] = f"Intensity {experiment}"

        if len(experiments) > 1 and not self.skip_lfq:
            for experiment in experiments:
                header_dict[f"{experiment} MaxLFQ Intensity"] = (
                    f"LFQ Intensity {experiment}"
                )

        # the iBAQ and peptide counts columns are not part of Fragpipe's format but were added upon popular demand
        for experiment in experiments:
            header_dict[f"{experiment} iBAQ Intensity"] = f"iBAQ {experiment}"

        for experiment in experiments:
            header_dict[f"{experiment} Unique Peptide Count"] = (
                f"Unique peptides {experiment}"
            )

        header_dict["Indistinguishable Proteins"] = "Indistinguishable Proteins"

        return header_dict

    def get_columns(self) -> List[columns.ProteinGroupColumns]:
        num_ibaq_peptides_per_protein = collections.defaultdict(lambda: 1)

        output_columns = [
            columns.FragpipeProteinAnnotationsColumns(
                self.protein_groups, self.protein_annotations
            ),
            columns.ProteinProbabilityColumns(),
            columns.TopPeptideProbabilityColumns(),
            columns.UniquePeptideCountColumns(),
            columns.SpectralCountColumns(),
            columns.SummedIntensityAndIbaqColumns(
                num_ibaq_peptides_per_protein, self.protein_group_fdr_threshold
            ),
            columns.IndistinguishableProteinsColumns(),
        ]
        if not self.skip_lfq:
            output_columns += [
                columns.LFQIntensityColumns(
                    self.min_peptide_ratios_lfq,
                    self.stabilize_large_ratios_lfq,
                    self.num_threads,
                    self.protein_group_fdr_threshold,
                )
            ]
        return output_columns

    def get_extra_columns_formatter(self):
        return fragpipe_format_extra_columns
