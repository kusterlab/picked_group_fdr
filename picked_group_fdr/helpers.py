import re
from typing import List, Set

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


def is_shared_peptide(protein_group_idxs: Set[int]) -> bool:
    """Checks if a peptide is assigned to multiple protein groups.

    Args:
        protein_group_idxs (Set[int]): Set of protein group indices, i.e. the index of
            the protein group in ProteinGroups object.

    Returns:
        bool: True if peptide is assigned to multiple protein groups.
    """
    return len(protein_group_idxs) > 1


def is_missing_in_protein_groups(protein_group_idxs: Set[int]) -> bool:
    """Checks if a peptide could not be assigned to any protein group.

    If a protein is not assigned to any protein group, a protein_group_idx of -1 is
    assigned (see ProteinGroups.get_protein_group_idxs). This for example happens for
    proteins that got outcompeted in the picked protein competition.

    Args:
        protein_group_idxs (Set[int]): Set of protein group indices, i.e. the index of
            the protein group in ProteinGroups object.

    Returns:
        bool: True if peptide is assigned to no protein groups
    """
    return len(protein_group_idxs) == 0 or protein_group_idxs == {-1}

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


def remove_modifications(modified_peptide: str):
    """
    Removes modifications in [] or ()

    Warning! Dirty hack: last removal of ) is needed because some MQ versions have nested parentheses for modifications
    """
    return re.sub(r"\[[^]]*\]", "", re.sub(r"\([^)]*\)", "", modified_peptide)).replace(")", "")


def chunks(lst: List, n: int):
    """
    Splits a list into sub lists of length n by taking n consecutive elements at a time
    """
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def string_to_bool(v: str):
    return v.lower() in ("yes", "true", "t", "1")
