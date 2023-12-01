from __future__ import annotations

from typing import Dict, List

from ..protein_annotation import ProteinAnnotation
from .. import columns


class MinimalProteinGroupsWriter:
    def __init__(self, protein_annotations: Dict[str, ProteinAnnotation]) -> None:
        self.protein_annotations = protein_annotations

    def get_columns(self) -> List[columns.ProteinGroupColumns]:
        return [
            columns.ProteinAnnotationsColumns(self.protein_annotations)
        ]
