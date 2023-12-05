from __future__ import annotations

import collections
from typing import Dict, List

from .base import ProteinGroupsWriter
from .fragpipe import fragpipe_format_extra_columns

from .. import columns
from ..protein_groups import ProteinGroups

# for type hints only
from .. import results
from .. import protein_annotation as pa

FRAGPIPE_PROTEIN_OUTPUT_DICT = {
    "Protein": "Protein",
    "Protein ID": "Protein ID",
    "Entry Name": "Entry Name",
    "Gene": "Gene",
    "Length": "Length",
    "Organism": "Organism",
    "Protein Description": "Protein Description",
    "Protein Existence": "Protein Existence",
    "Coverage": "Sequence coverage [%]",
    "Protein Probability": "Protein Probability",
    "Top Peptide Probability": "Top Peptide Probability",
    "Total Peptides": "Unique peptides 1",
    "Unique Peptides": "Unique peptides 1",
    "Razor Peptides": "Unique peptides 1",
    "Total Spectral Count": "Spectral count 1",
    "Unique Spectral Count": "Spectral count 1",
    "Razor Spectral Count": "Spectral count 1",
    "Total Intensity": "Intensity 1",
    "Unique Intensity": "Intensity 1",
    "Razor Intensity": "Intensity 1",
    "Razor Assigned Modifications": "Razor Assigned Modifications",
    "Razor Observed Modifications": "Razor Observed Modifications",
    "Indistinguishable Proteins": "Indistinguishable Proteins",
}


class FragPipeSingleProteinWriter(ProteinGroupsWriter):
    def __init__(
        self,
        protein_groups: ProteinGroups,
        protein_annotations: Dict[str, pa.ProteinAnnotation],
        protein_sequences: Dict[str, str],
    ):
        self.protein_groups = protein_groups
        self.protein_annotations = protein_annotations
        self.protein_sequences = protein_sequences

    def get_header_dict(
        self, protein_group_results: results.ProteinGroupResults
    ) -> Dict[str, str]:
        return FRAGPIPE_PROTEIN_OUTPUT_DICT

    def get_columns(self) -> List[columns.ProteinGroupColumns]:
        num_ibaq_peptides_per_protein = collections.defaultdict(lambda: 1)
        return [
            columns.FragpipeProteinAnnotationsColumns(
                self.protein_groups, self.protein_annotations
            ),
            columns.SequenceCoverageColumns(self.protein_sequences),
            columns.ProteinProbabilityColumns(),
            columns.TopPeptideProbabilityColumns(),
            columns.UniquePeptideCountColumns(),
            columns.SpectralCountColumns(),
            columns.SummedIntensityAndIbaqColumns(num_ibaq_peptides_per_protein),
            columns.ModificationsColumns(),
            columns.IndistinguishableProteinsColumns(),
        ]

    def get_extra_columns_formatter(self):
        return fragpipe_format_extra_columns
