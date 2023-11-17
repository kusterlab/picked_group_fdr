import logging

from . import tsv

logger = logging.getLogger(__name__)


def is_fragpipe_psm_file(headers):
    return "peptideprophet probability" in map(str.lower, headers)


def parse_fragpipe_psm_file(reader, headers, get_proteins, score_type="PEP"):
    pept_col = tsv.get_column_index(headers, "Peptide")
    score_col = tsv.get_column_index(headers, "SpectralSim") # could also use Hyperscore
    post_err_prob_col = tsv.get_column_index(headers, "PeptideProphet Probability")
    protein_col = tsv.get_column_index(headers, "Protein")
    other_proteins_col = tsv.get_column_index(headers, "Mapped Proteins")

    if score_type.get_score_column() == "pep":
        score_col = post_err_prob_col

    logger.info("Parsing FragPipe psm.tsv file")
    for line_idx, row in enumerate(reader):
        if line_idx % 500000 == 0:
            logger.info(f"    Reading line {line_idx}")

        peptide = row[pept_col]
        experiment = 1
        score = float(row[score_col])
        if score_type.get_score_column() == "pep":
            score = 1 - score + 1e-16

        proteins = [row[protein_col]]
        if len(row[other_proteins_col]) > 0:
            proteins += row[other_proteins_col].split(", ")

        proteins = get_proteins(peptide, proteins)
        if proteins:
            yield peptide, proteins, experiment, score


def parse_fragpipe_psm_file_for_peptide_remapping(reader, headers):
    protein_col = tsv.get_column_index(headers, "Protein")
    other_proteins_col = tsv.get_column_index(headers, "Mapped Proteins")

    logger.info("Parsing FragPipe psm.tsv file")
    for line_idx, row in enumerate(reader):
        if line_idx % 500000 == 0:
            logger.info(f"    Reading line {line_idx}")

        proteins = [row[protein_col]]
        if len(row[other_proteins_col]) > 0:
            proteins += row[other_proteins_col].split(", ")

        yield row, proteins

        


