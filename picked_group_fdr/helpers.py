import re
from typing import List, Dict, Optional

import numpy as np


def allContain(proteinGroup: List[str], prefix) -> bool:
    return len([1 for x in proteinGroup if prefix in x]) == len(proteinGroup)


def is_contaminant(proteinGroup: List[str]) -> bool:
    return allContain(proteinGroup, "CON__")

    
def isDecoy(proteinGroup: List[str]) -> bool:
    return allContain(proteinGroup, "REV__") or allContain(proteinGroup, "rev_")


def isObsolete(proteinGroup: List[str]) -> bool:
    return allContain(proteinGroup, "OBSOLETE__")


def isMbr(postErrProb: float) -> bool:
    return np.isnan(postErrProb) # match-between-runs have NaNs in the PEP column


def is_shared_peptide(proteinGroupIdxs: List[int]) -> bool:
    return len(proteinGroupIdxs) > 1 or proteinGroupIdxs == [-1]


def remove_decoy_proteins_from_target_peptides(proteins: List[str]) -> List[str]:
    """
    shared peptides with at least one target protein are only associated with 
    their target proteins. This prevents mixed target-decoy protein groups and 
    undesired situations where a decoy protein in a target protein group 
    outcompetes its target counterpart
    """
    isTargetPeptide = not isDecoy(proteins)
    new_proteins = list()
    for protein in proteins:
        if isTargetPeptide and (protein.startswith("REV__") or protein.startswith("rev_")):
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
    return re.sub(r'\[[^]]*\]', '', re.sub(r'\([^)]*\)', '', peptide)).replace(")", "") 


def chunks(lst: List, n: int):
    """
    Splits a list into sub lists of length n by taking n consecutive elements at a time
    """
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def string_to_bool(v: str):
    return v.lower() in ("yes", "true", "t", "1")