from __future__ import annotations

import sys
from typing import Any, Dict, List

# for type hints only
from .. import quant

MQ_PROTEIN_ANNOTATION_HEADERS = ["Protein names", "Gene names", "Fasta headers"]


def get_mq_protein_groups_columns(
    num_ibaq_peptides_per_protein: Dict[str, int],
    protein_sequences: Dict[str, str],
    num_tmt_channels: int,
    num_silac_channels: int,
    min_peptide_ratios_lfq: int,
    stabilize_large_ratios_lfq: bool,
    num_threads: int,
    params: Dict[str, Any],
) -> List[quant.ProteinGroupColumns]:
    silac_channels = get_silac_channels(num_silac_channels)

    columns: List[quant.ProteinGroupColumns] = [
        quant.UniquePeptideCountColumns(),
        quant.IdentificationTypeColumns(),
        quant.SummedIntensityAndIbaqColumns(
            silac_channels, num_ibaq_peptides_per_protein
        ),
        quant.SequenceCoverageColumns(protein_sequences),
        quant.EvidenceIdsColumns(),
    ]

    if num_tmt_channels > 0:
        columns.append(quant.TMTIntensityColumns(num_tmt_channels))
    else:
        columns.append(
            quant.LFQIntensityColumns(
                silac_channels,
                min_peptide_ratios_lfq,
                stabilize_large_ratios_lfq,
                num_threads,
            )
        )
        # TODO: add SILAC functionality of Triqler
        if num_silac_channels == 0:
            columns.append(quant.TriqlerIntensityColumns(params))
    return columns


def get_silac_channels(num_silac_channels: int):
    silac_channels = list()
    if num_silac_channels == 3:
        silac_channels = ["L", "M", "H"]
    elif num_silac_channels == 2:
        silac_channels = ["L", "H"]
    elif num_silac_channels != 0:
        sys.exit("ERROR: Found a number of SILAC channels not equal to 2 or 3")
    return silac_channels
