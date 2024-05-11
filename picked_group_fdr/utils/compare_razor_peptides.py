from picked_group_fdr import digest, helpers
from picked_group_fdr.parsers import psm
from picked_group_fdr.scoring import logger
from picked_group_fdr.scoring_strategy import ProteinScoringStrategy


def compare_razor_peptides(
    mq_evidence_file, peptide_to_protein_map, protein_groups, score_type
):
    """Compares the chosen protein by MaxQuant according to the razor peptide rule with our implementation of the razor peptide rule"""
    score_type = ProteinScoringStrategy("multPEP razor")
    for peptide_row in psm.parse_evidence_file_single(
        mq_evidence_file, score_type=score_type
    ):
        modified_peptide, tmp_proteins, _, _ = peptide_row

        proteins = digest.get_proteins(
            peptide_to_protein_map, helpers.remove_modifications(modified_peptide)
        )

        leading_proteins = protein_groups.get_leading_proteins(proteins)
        if len(leading_proteins) > 1:
            predicted_razor = score_type.filter_proteins(
                proteins
            )  # filtering for razor peptide approach
            logger.debug(
                f"{tmp_proteins[0] == predicted_razor[0]} {tmp_proteins[0]} {predicted_razor[0]}"
            )
