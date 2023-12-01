import collections
from typing import Dict, List

from .. import quant
from ..protein_annotation import ProteinAnnotation
from ..protein_groups import ProteinGroups


FRAGPIPE_PROTEIN_OUTPUT_DICT = {
    "Protein": "Protein",
    "Protein ID": "Protein ID",
    "Entry Name": "Entry Name",
    "Gene": "Gene",
    "Length": "Length",
    "Organism": "Organism",
    "Protein Description": "Protein Description",
    "Protein Existence": "Protein Existence",
    "Coverage": "Sequence coverage [%]",
    "Protein Probability": "Protein Probability",
    "Top Peptide Probability": "Top Peptide Probability",
    "Total Peptides": "Unique peptides 1",
    "Unique Peptides": "Unique peptides 1",
    "Razor Peptides": "Unique peptides 1",
    "Total Spectral Count": "Spectral count 1",
    "Unique Spectral Count": "Spectral count 1",
    "Razor Spectral Count": "Spectral count 1",
    "Total Intensity": "Intensity 1",
    "Unique Intensity": "Intensity 1",
    "Razor Intensity": "Intensity 1",
    "Razor Assigned Modifications": "Razor Assigned Modifications",
    "Razor Observed Modifications": "Razor Observed Modifications",
    "Indistinguishable Proteins": "Indistinguishable Proteins",
}


def get_fragpipe_protein_tsv_columns(
    protein_groups: ProteinGroups,
    protein_annotations: Dict[str, ProteinAnnotation],
    protein_sequences: Dict[str, str],
):
    silac_channels = []
    num_ibaq_peptides_per_protein = collections.defaultdict(lambda: 1)

    columns: List[quant.ProteinGroupColumns] = [
        quant.FragpipeProteinAnnotationsColumns(protein_groups, protein_annotations),
        quant.SequenceCoverageColumns(protein_sequences),
        quant.ProteinProbabilityColumns(),
        quant.TopPeptideProbabilityColumns(),
        quant.UniquePeptideCountColumns(),
        quant.SpectralCountColumns(),
        quant.SummedIntensityAndIbaqColumns(num_ibaq_peptides_per_protein),
        quant.ModificationsColumns(),
        quant.IndistinguishableProteinsColumns(),
    ]

    return columns
