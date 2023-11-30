from __future__ import annotations

from typing import Dict, List

from ..protein_annotation import ProteinAnnotation
from .. import quant


def get_minimal_protein_groups_columns(
    protein_annotations: Dict[str, ProteinAnnotation],
) -> List[quant.ProteinGroupColumns]:
    columns: List[quant.ProteinGroupColumns] = [
        quant.ProteinAnnotationsColumns(protein_annotations)
    ]
    return columns
