from __future__ import annotations

from typing import Dict, List

from .base import ProteinGroupsWriter
from .fragpipe import fragpipe_format_extra_columns

from .. import columns

# for type hints only
from .. import results
from .. import protein_annotation as pa

DIANN_PROTEIN_OUTPUT_DICT = {
    "Protein.Group": "Protein.Group",
    "Protein.Names": "Protein.Names",
    "Genes": "Genes",
    "First.Protein.Description": "First.Protein.Description",
    "N.Sequences": "Combined Total Peptides",
    "N.Proteotypic.Sequences": "Combined Total Peptides",
}


class DiannProteinGroupsWriter(ProteinGroupsWriter):
    def __init__(
        self,
        protein_annotations: Dict[str, pa.ProteinAnnotation],
        min_peptide_ratios_lfq: int,
        stabilize_large_ratios_lfq: bool,
        fast_lfq: bool,
        num_threads: int,
        protein_group_fdr_threshold: float = 0.01,
    ) -> None:
        self.protein_annotations = protein_annotations
        self.min_peptide_ratios_lfq = min_peptide_ratios_lfq
        self.stabilize_large_ratios_lfq = stabilize_large_ratios_lfq
        self.fast_lfq = fast_lfq
        self.num_threads = num_threads

        # only used for reporting number of protein groups at the given threshold
        self.protein_group_fdr_threshold = protein_group_fdr_threshold

    def get_header_dict(
        self, protein_group_results: results.ProteinGroupResults
    ) -> Dict[str, str]:
        experiments = protein_group_results.experiments
        header_dict = DIANN_PROTEIN_OUTPUT_DICT.copy()
        return header_dict | {
            experiment: f"LFQ Intensity {experiment}" for experiment in experiments
        }

    def get_columns(self) -> List[columns.ProteinGroupColumns]:
        return [
            columns.DiannProteinAnnotationsColumns(self.protein_annotations),
            columns.UniquePeptideCountColumns(),
            columns.LFQIntensityColumns(
                min_peptide_ratios_lfq=self.min_peptide_ratios_lfq,
                stabilize_large_ratios_lfq=self.stabilize_large_ratios_lfq,
                fast_lfq=self.fast_lfq,
                num_threads=self.num_threads,
                protein_group_fdr_threshold=self.protein_group_fdr_threshold,
            ),
        ]

    def get_extra_columns_formatter(self):
        return fragpipe_format_extra_columns
