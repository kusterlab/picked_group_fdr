from __future__ import annotations

from typing import Dict, List

from .base import ProteinGroupsWriter
from .. import columns

# for type hints only
from .. import results
from .. import protein_annotation as pa


class MinimalProteinGroupsWriter(ProteinGroupsWriter):
    def __init__(self, protein_annotations: Dict[str, pa.ProteinAnnotation]) -> None:
        self.protein_annotations = protein_annotations

    def get_columns(self) -> List[columns.ProteinGroupColumns]:
        return [columns.ProteinAnnotationsColumns(self.protein_annotations)]

    def append_quant_columns(
        self,
        protein_group_results: results.ProteinGroupResults,
        post_err_probs: List,
        psm_fdr_cutoff: float,
    ):
        for c in self.get_columns():
            c.append(protein_group_results, None)
        return protein_group_results
