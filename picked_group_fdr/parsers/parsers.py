from typing import Dict, List
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from . import tsv

logger = logging.getLogger(__name__)


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


def add_triqler_group_params(experimental_design: pd.DataFrame, params: Dict):
    # Create formatted group names, e.g. 1:group1
    conditions = experimental_design["Condition"].unique()
    groupNames = [
        f"{groupIdx + 1}:{groupName}" for groupIdx, groupName in enumerate(conditions)
    ]
    experimental_design["Condition"] = experimental_design["Condition"].map(
        lambda x: groupNames[conditions.tolist().index(x)]
    )

    experiments = experimental_design["Experiment"].unique().tolist()

    for _, condition, experiment, _ in experimental_design.itertuples(
        index=False, name=None
    ):
        if condition not in params["groupLabels"]:
            params["groupLabels"].append(condition)
            params["groups"].append([])

        params["groups"][params["groupLabels"].index(condition)].append(
            experiments.index(experiment)
        )
    return params


def get_file_mapping(file_info_list: pd.DataFrame):
    dict_with_keys = file_info_list.set_index("Name")[
        ["Experiment", "Fraction"]
    ].to_dict(orient="index")
    return {key: tuple(value.values()) for key, value in dict_with_keys.items()}


def parse_triqler_file_list(input_file: str):
    column_names = ["Name", "Condition", "Experiment", "Fraction"]
    df = pd.read_csv(input_file, sep="\t", header=None, names=column_names)
    return normalize_experimental_design(df)


def parse_mq_experimental_design(input_file: str):
    df = pd.read_csv(input_file, sep="\t")
    if "Condition" not in df.columns:
        df["Condition"] = df["Experiment"]
    df = df[["Name", "Condition", "Experiment", "Fraction"]]
    return normalize_experimental_design(df)


def normalize_experimental_design(df: pd.DataFrame) -> pd.DataFrame:
    # Extract file names without extensions
    df["Name"] = df["Name"].apply(lambda x: Path(x).stem)
    df["Experiment"] = df["Experiment"].fillna(df["Name"])
    df["Fraction"] = df["Fraction"].fillna(-1)
    df["Condition"] = df["Condition"].fillna(df["Experiment"])

    return df


