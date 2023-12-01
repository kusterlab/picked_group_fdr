from __future__ import annotations

from typing import Dict, List

from ..protein_annotation import ProteinAnnotation
from .. import columns


def get_minimal_protein_groups_columns(
    protein_annotations: Dict[str, ProteinAnnotation],
) -> List[columns.ProteinGroupColumns]:
    return [
        columns.ProteinAnnotationsColumns(protein_annotations)
    ]
