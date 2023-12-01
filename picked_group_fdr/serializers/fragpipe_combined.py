from __future__ import annotations

import collections
from typing import Dict, List

from .. import columns
from ..protein_annotation import ProteinAnnotation
from ..protein_groups import ProteinGroups


FRAGPIPE_COMBINED_PROTEIN_OUTPUT_DICT = {
    "Protein": "Protein",
    "Protein group": "Protein IDs",
    "Protein ID": "Protein ID",
    "Entry Name": "Entry Name",
    "Gene": "Gene",
    "Protein Length": "Length",
    "Organism": "Organism",
    "Protein Existence": "Protein Existence",
    "Description": "Protein Description",
    "Combined Total Peptides": "Combined Total Peptides",
    "Combined Spectral Count": "Combined Spectral Count",
    "Combined Unique Spectral Count": "Combined Spectral Count",
    "Combined Total Spectral Count": "Combined Spectral Count",
}


def get_fragpipe_combined_protein_columns(
    protein_groups: ProteinGroups,
    protein_annotations: Dict[str, ProteinAnnotation],
    min_peptide_ratios_lfq: int = 1,
    stabilize_large_ratios_lfq: bool = True,
    num_threads: int = 1,
) -> List[columns.ProteinGroupColumns]:
    silac_channels = []
    num_ibaq_peptides_per_protein = collections.defaultdict(lambda: 1)

    return [
        columns.FragpipeProteinAnnotationsColumns(protein_groups, protein_annotations),
        columns.ProteinProbabilityColumns(),
        columns.TopPeptideProbabilityColumns(),
        columns.UniquePeptideCountColumns(),
        columns.SpectralCountColumns(),
        columns.SummedIntensityAndIbaqColumns(num_ibaq_peptides_per_protein),
        columns.LFQIntensityColumns(
            min_peptide_ratios_lfq,
            stabilize_large_ratios_lfq,
            num_threads,
        ),
        columns.IndistinguishableProteinsColumns(),
    ]


def get_fragpipe_combined_protein_headers(experiments: List[str]):
    """Adds experiment specific headers.

    - <Experiment> Spectral Count
    - <Experiment> Unique Spectral Count
    - <Experiment> Total Spectral Count
    - <Experiment> Intensity
    - <Experiment> MaxLFQ Intensity
    """
    header_dict = FRAGPIPE_COMBINED_PROTEIN_OUTPUT_DICT.copy()
    for experiment in experiments:
        header_dict[f"{experiment} Spectra Count"] = f"Spectral count {experiment}"

    for experiment in experiments:
        header_dict[
            f"{experiment} Unique Spectra Count"
        ] = f"Spectral count {experiment}"

    for experiment in experiments:
        header_dict[
            f"{experiment} Total Spectra Count"
        ] = f"Spectral count {experiment}"

    for experiment in experiments:
        header_dict[f"{experiment} Intensity"] = f"Intensity {experiment}"

    for experiment in experiments:
        header_dict[f"{experiment} MaxLFQ Intensity"] = f"LFQ Intensity {experiment}"

    header_dict["Indistinguishable Proteins"] = "Indistinguishable Proteins"

    return header_dict
