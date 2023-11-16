from typing import List
import logging

import numpy as np

from . import tsv


logger = logging.getLogger(__name__)


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
