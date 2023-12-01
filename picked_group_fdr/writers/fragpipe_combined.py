from __future__ import annotations

import collections
from typing import Dict, List

from .. import columns
from ..protein_annotation import ProteinAnnotation
from ..protein_groups import ProteinGroups

from .fragpipe import fragpipe_format_extra_columns


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
    "Combined Total Peptides": "Combined Total Peptides",
    "Combined Spectral Count": "Combined Spectral Count",
    "Combined Unique Spectral Count": "Combined Spectral Count",
    "Combined Total Spectral Count": "Combined Spectral Count",
}


class FragPipeCombinedProteinWriter:
    def __init__(
        self,
        protein_groups: ProteinGroups,
        protein_annotations: Dict[str, ProteinAnnotation],
        min_peptide_ratios_lfq: int = 1,
        stabilize_large_ratios_lfq: bool = True,
        num_threads: int = 1,
    ) -> None:
        self.protein_groups = protein_groups
        self.protein_annotations = protein_annotations
        self.min_peptide_ratios_lfq = min_peptide_ratios_lfq
        self.stabilize_large_ratios_lfq = stabilize_large_ratios_lfq
        self.num_threads = num_threads

    def get_header_dict(self, experiments: List[str]) -> Dict[str, str]:
        """Adds experiment specific headers.

        - <Experiment> Spectral Count
        - <Experiment> Unique Spectral Count
        - <Experiment> Total Spectral Count
        - <Experiment> Intensity
        - <Experiment> MaxLFQ Intensity
        """
        header_dict = FRAGPIPE_COMBINED_PROTEIN_OUTPUT_DICT.copy()
        for experiment in experiments:
            header_dict[f"{experiment} Spectra Count"] = f"Spectral count {experiment}"

        for experiment in experiments:
            header_dict[
                f"{experiment} Unique Spectra Count"
            ] = f"Spectral count {experiment}"

        for experiment in experiments:
            header_dict[
                f"{experiment} Total Spectra Count"
            ] = f"Spectral count {experiment}"

        for experiment in experiments:
            header_dict[f"{experiment} Intensity"] = f"Intensity {experiment}"

        for experiment in experiments:
            header_dict[f"{experiment} MaxLFQ Intensity"] = f"LFQ Intensity {experiment}"

        header_dict["Indistinguishable Proteins"] = "Indistinguishable Proteins"

        return header_dict
        
    def get_columns(self) -> List[columns.ProteinGroupColumns]:
        num_ibaq_peptides_per_protein = collections.defaultdict(lambda: 1)

        return [
            columns.FragpipeProteinAnnotationsColumns(
                self.protein_groups, self.protein_annotations
            ),
            columns.ProteinProbabilityColumns(),
            columns.TopPeptideProbabilityColumns(),
            columns.UniquePeptideCountColumns(),
            columns.SpectralCountColumns(),
            columns.SummedIntensityAndIbaqColumns(num_ibaq_peptides_per_protein),
            columns.LFQIntensityColumns(
                self.min_peptide_ratios_lfq,
                self.stabilize_large_ratios_lfq,
                self.num_threads,
            ),
            columns.IndistinguishableProteinsColumns(),
        ]
    
    def get_extra_columns_formatter(self):
        return fragpipe_format_extra_columns