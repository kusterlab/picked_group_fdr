from __future__ import annotations

from typing import Any, Dict, List

# for type hints only
from .. import quant


def get_mq_protein_groups_columns(
    num_ibaq_peptides_per_protein: Dict[str, int],
    protein_sequences: Dict[str, str],
    min_peptide_ratios_lfq: int,
    stabilize_large_ratios_lfq: bool,
    num_threads: int,
    params: Dict[str, Any],
) -> List[quant.ProteinGroupColumns]:
    columns: List[quant.ProteinGroupColumns] = [
        quant.UniquePeptideCountColumns(),
        quant.IdentificationTypeColumns(),
        quant.SummedIntensityAndIbaqColumns(
            num_ibaq_peptides_per_protein
        ),
        quant.SequenceCoverageColumns(protein_sequences),
        quant.EvidenceIdsColumns(),
        quant.LFQIntensityColumns(
            min_peptide_ratios_lfq,
            stabilize_large_ratios_lfq,
            num_threads,
        ),
        quant.TMTIntensityColumns(),
        quant.TriqlerIntensityColumns(params)
    ]
    return columns


