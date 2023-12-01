from __future__ import annotations

from typing import Any, Dict, List

# for type hints only
from .. import columns


class MaxQuantProteinGroupsWriter:
    def __init__(
        self,
        num_ibaq_peptides_per_protein: Dict[str, int],
        protein_sequences: Dict[str, str],
        min_peptide_ratios_lfq: int,
        stabilize_large_ratios_lfq: bool,
        num_threads: int,
        params: Dict[str, Any],
    ) -> None:
        self.num_ibaq_peptides_per_protein = num_ibaq_peptides_per_protein
        self.protein_sequences = protein_sequences
        self.min_peptide_ratios_lfq = min_peptide_ratios_lfq
        self.stabilize_large_ratios_lfq = stabilize_large_ratios_lfq
        self.num_threads = num_threads
        self.params = params

    def get_columns(self) -> List[columns.ProteinGroupColumns]:
        return [
            columns.UniquePeptideCountColumns(),
            columns.IdentificationTypeColumns(),
            columns.SummedIntensityAndIbaqColumns(self.num_ibaq_peptides_per_protein),
            columns.SequenceCoverageColumns(self.protein_sequences),
            columns.EvidenceIdsColumns(),
            columns.LFQIntensityColumns(
                self.min_peptide_ratios_lfq,
                self.stabilize_large_ratios_lfq,
                self.num_threads,
            ),
            columns.TMTIntensityColumns(),
            columns.TriqlerIntensityColumns(self.params),
        ]
