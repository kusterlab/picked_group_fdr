from __future__ import annotations

from typing import Any, Dict, List

# for type hints only
from .. import columns


def get_mq_protein_groups_columns(
    num_ibaq_peptides_per_protein: Dict[str, int],
    protein_sequences: Dict[str, str],
    min_peptide_ratios_lfq: int,
    stabilize_large_ratios_lfq: bool,
    num_threads: int,
    params: Dict[str, Any],
) -> List[columns.ProteinGroupColumns]:
    return [
        columns.UniquePeptideCountColumns(),
        columns.IdentificationTypeColumns(),
        columns.SummedIntensityAndIbaqColumns(
            num_ibaq_peptides_per_protein
        ),
        columns.SequenceCoverageColumns(protein_sequences),
        columns.EvidenceIdsColumns(),
        columns.LFQIntensityColumns(
            min_peptide_ratios_lfq,
            stabilize_large_ratios_lfq,
            num_threads,
        ),
        columns.TMTIntensityColumns(),
        columns.TriqlerIntensityColumns(params)
    ]


