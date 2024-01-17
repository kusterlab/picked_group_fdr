import re
from typing import List

import numpy as np


def _all_contain(protein_group: List[str], prefix) -> bool:
    return len([1 for x in protein_group if prefix in x]) == len(protein_group)


def is_contaminant(protein_group: List[str]) -> bool:
    return _all_contain(protein_group, "CON__")


def is_decoy(protein_group: List[str]) -> bool:
    return _all_contain(protein_group, "REV__") or _all_contain(protein_group, "rev_")


def is_obsolete(protein_group: List[str]) -> bool:
    return _all_contain(protein_group, "OBSOLETE__")


def is_mbr(post_err_prob: float) -> bool:
    return np.isnan(post_err_prob)  # match-between-runs have NaNs in the PEP column


def is_shared_peptide(protein_group_idxs: List[int]) -> bool:
    return len(protein_group_idxs) > 1 or protein_group_idxs == [-1]


def remove_decoy_proteins_from_target_peptides(proteins: List[str]) -> List[str]:
    """
    shared peptides with at least one target protein are only associated with
    their target proteins. This prevents mixed target-decoy protein groups and
    undesired situations where a decoy protein in a target protein group
    outcompetes its target counterpart
    """
    is_target_peptide = not is_decoy(proteins)
    new_proteins = list()
    for protein in proteins:
        if is_target_peptide and (
            protein.startswith("REV__") or protein.startswith("rev_")
        ):
            continue
        else:
            new_proteins.append(protein)
    return new_proteins


def clean_peptide(peptide: str, remove_flanks: bool = True):
    """
    Removes modifications in [] or () as well as the _ that MaxQuant puts in front and after the sequence

    Warning! Dirty hack: last removal of ) is needed because some MQ versions have nested parentheses for modifications
    """
    if remove_flanks:
        peptide = peptide[1:-1]
    return re.sub(r"\[[^]]*\]", "", re.sub(r"\([^)]*\)", "", peptide)).replace(")", "")


def chunks(lst: List, n: int):
    """
    Splits a list into sub lists of length n by taking n consecutive elements at a time
    """
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def string_to_bool(v: str):
    return v.lower() in ("yes", "true", "t", "1")
