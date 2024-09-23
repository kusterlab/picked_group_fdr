from __future__ import annotations

from typing import Any, Dict, List

from .base import ProteinGroupsWriter

# for type hints only
from .. import protein_annotation as pa
from .. import columns


class MaxQuantProteinGroupsWriter(ProteinGroupsWriter):
    def __init__(
        self,
        num_ibaq_peptides_per_protein: Dict[str, int],
        protein_annotations: Dict[str, pa.ProteinAnnotation],
        protein_sequences: Dict[str, str],
        skip_lfq: bool,
        min_peptide_ratios_lfq: int,
        stabilize_large_ratios_lfq: bool,
        num_threads: int,
        params: Dict[str, Any],
        protein_group_fdr_threshold: float = 0.01,
    ) -> None:
        self.num_ibaq_peptides_per_protein = num_ibaq_peptides_per_protein
        self.protein_annotations = protein_annotations
        self.protein_sequences = protein_sequences
        self.skip_lfq = skip_lfq
        self.min_peptide_ratios_lfq = min_peptide_ratios_lfq
        self.stabilize_large_ratios_lfq = stabilize_large_ratios_lfq
        self.num_threads = num_threads
        self.params = params

        # only used for reporting number of protein groups at the given threshold
        self.protein_group_fdr_threshold = protein_group_fdr_threshold

    def get_columns(self) -> List[columns.ProteinGroupColumns]:
        output_columns = [
            columns.ProteinAnnotationsColumns(self.protein_annotations),
            columns.UniquePeptideCountColumns(),
            columns.IdentificationTypeColumns(),
            columns.SummedIntensityAndIbaqColumns(
                self.num_ibaq_peptides_per_protein, self.protein_group_fdr_threshold
            ),
            columns.SequenceCoverageColumns(self.protein_sequences),
            columns.EvidenceIdsColumns(),
            columns.TMTIntensityColumns(),
            columns.TriqlerIntensityColumns(self.params),
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
