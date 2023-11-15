import csv
from typing import Dict, List
import logging

import numpy as np

from .. import digest, helpers
from . import tsv
from . import maxquant
from . import percolator

# for type hints only
from ..scoring import ProteinScoringStrategy

logger = logging.getLogger(__name__)


# csv.field_size_limit(sys.maxsize)
csv.field_size_limit(2147483647)


def parse_protein_groups_file_multiple(
    protein_groups_files: List[str], are_decoy_file: List[bool], **kwargs
):
    for protein_groups_file, is_decoy_file in zip(protein_groups_files, are_decoy_file):
        yield from parse_protein_groups_file_single(
            protein_groups_file, is_decoy_file=is_decoy_file, **kwargs
        )


def parse_protein_groups_file_single(
    protein_groups_file: str,
    protein_column: str = "Protein IDs",
    score_column: str = "Score",
    is_decoy_file: bool = False,
):
    """Parse protein groups file from MaxQuant or ProteomeDiscoverer.

    For PD with Chimerys, the column names are:
    - protein_column="Accession"
    - score_column="Score CHIMERY CHIMERYS"

    Args:
        protein_groups_file (_type_): _description_
        protein_column (str, optional): _description_. Defaults to 'Protein IDs'.
        score_column (str, optional): _description_. Defaults to 'Score'.
        is_decoy_file (bool, optional): _description_. Defaults to False.

    Yields:
        _type_: _description_
    """
    delimiter = tsv.get_delimiter(protein_groups_file)

    reader = tsv.get_tsv_reader(protein_groups_file, delimiter)
    headers = next(reader)  # save the header

    score_col = tsv.get_column_index(headers, score_column)
    protein_col = tsv.get_column_index(headers, protein_column)

    logger.info(f"Parsing proteinGroups file: {protein_groups_file}")
    for row in reader:
        proteins = list(map(str.strip, row[protein_col].split(";")))
        if is_decoy_file:
            proteins = [f"REV__{p}" for p in proteins]
        score = -100.0
        if len(row[score_col]) > 0:
            score = float(row[score_col])
        yield proteins, score


def parse_peptides_files_multiple(
    peptides_files: List[str], are_decoy_file: List[bool], **kwargs
):
    for peptide_file, is_decoy_file in zip(peptides_files, are_decoy_file):
        yield from parse_peptides_file_single(
            peptide_file, is_decoy_file=is_decoy_file, **kwargs
        )


def parse_peptides_file_single(
    peptides_file: str,
    peptide_column: str = "Modified sequence",
    protein_column: str = "Protein IDs",
    score_column: str = "Score",
    is_decoy_file: bool = False,
):
    """Parse peptide-level file from MaxQuant (evidence/msms.txt) or ProteomeDiscoverer.

    For PD with Chimerys, the column names are:
    - protein_column="Accession"
    - score_column="Score CHIMERY CHIMERYS"

    Args:
        peptides_file (_type_): evidence.txt or msms.txt
        peptide_column (str, optional): column name for peptide sequence. Defaults to 'Protein IDs'.
        protein_column (str, optional): column name for protein identifiers. Defaults to 'Protein IDs'.
        score_column (str, optional): column name for score. Defaults to 'Score'.
        is_decoy_file (bool, optional): if this file only contains decoy peptides. Defaults to False.

    Yields:
        (str, str, str, float): Peptide, Proteins, Experiment, Score
    """
    delimiter = tsv.get_delimiter(peptides_file)

    reader = tsv.get_tsv_reader(peptides_file, delimiter)
    headers = next(reader)  # save the header

    peptide_col = tsv.get_column_index(headers, peptide_column)
    score_col = tsv.get_column_index(headers, score_column)
    if not is_decoy_file:
        protein_col = tsv.get_column_index(headers, protein_column)

    experiment = "Experiment1"

    logger.info(f"Parsing peptides file: {peptides_file}")
    for row in reader:
        proteins = []
        if is_decoy_file:
            proteins = ["REV__protein"]
        elif len(row[protein_col]) > 0:
            proteins = list(map(str.strip, row[protein_col].split(";")))

        score = -100.0
        if len(row[score_col]) > 0:
            score = float(row[score_col])
            if np.isnan(score):
                continue
        yield row[peptide_col], proteins, experiment, score


def parse_evidence_file_multiple(
    evidence_files: List[str],
    peptide_to_protein_maps: List[Dict],
    score_type: ProteinScoringStrategy,
    for_quantification: bool = False,
    suppress_missing_peptide_warning: bool = False,
):
    for evidence_file, peptide_to_protein_map in zip(
        evidence_files, peptide_to_protein_maps
    ):
        yield from parse_evidence_file_single(
            evidence_file,
            peptide_to_protein_map,
            score_type,
            for_quantification,
            suppress_missing_peptide_warning,
        )


def parse_evidence_file_single(
    evidence_file: str,
    peptide_to_protein_map: Dict,
    score_type: ProteinScoringStrategy,
    for_quantification: bool = False,
    suppress_missing_peptide_warning: bool = False,
):
    delimiter = tsv.get_delimiter(evidence_file)
    reader = tsv.get_tsv_reader(evidence_file, delimiter)
    headers = next(reader)

    get_proteins = get_peptide_to_protein_mapper(
        peptide_to_protein_map, score_type, suppress_missing_peptide_warning
    )

    if percolator.is_percolator_file(headers):
        yield from percolator.parse_percolator_out_file(
            reader, headers, get_proteins, score_type
        )
    else:
        # convert headers to lowercase since MQ changes the capitalization frequently
        headers = list(map(str.lower, headers))
        if evidence_file.endswith(".csv"):
            headers = [x.replace(".", " ") for x in headers]
        yield from maxquant.parse_mq_evidence_file(
            reader, headers, get_proteins, score_type, for_quantification
        )


def get_peptide_to_protein_mapper(
    peptide_to_protein_map: Dict,
    score_type: ProteinScoringStrategy,
    suppress_missing_peptide_warning: bool,
):
    def get_proteins(peptide, tmp_proteins):
        if score_type.remaps_peptides_to_proteins():
            proteins = digest.get_proteins(peptide_to_protein_map, peptide)
            if len(proteins) == 0:
                if (
                    not helpers.is_contaminant(tmp_proteins)
                    and not suppress_missing_peptide_warning
                ):
                    logger.warning(f"Missing peptide: {peptide} {tmp_proteins}")
                return None
        else:
            proteins = tmp_proteins

        # filtering for razor peptide approach
        proteins = score_type.filter_proteins(proteins)

        return helpers.remove_decoy_proteins_from_target_peptides(proteins)

    return get_proteins
