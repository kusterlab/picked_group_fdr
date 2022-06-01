import re
from typing import List, Dict, Optional

import numpy as np


def isContaminant(proteinGroup: List[str]) -> bool:
  return len([1 for x in proteinGroup if "CON__" in x]) == len(proteinGroup)

  
def isDecoy(proteinGroup: List[str]) -> bool:
  return len([1 for x in proteinGroup if "REV__" in x]) == len(proteinGroup)


def isMbr(postErrProb: float) -> bool:
  return np.isnan(postErrProb) # match-between-runs have NaNs in the PEP column


def isSharedPeptide(proteinGroupIdxs: List[int]) -> bool:
  return len(proteinGroupIdxs) > 1


def removeDecoyProteinsFromTargetPeptides(proteins: List[str]) -> List[str]:
  """
  shared peptides with at least one target protein are only associated with 
  their target proteins. This prevents mixed target-decoy protein groups and 
  undesired situations where a decoy protein in a target protein group 
  outcompetes its target counterpart
  """
  isTargetPeptide = not isDecoy(proteins)
  new_proteins = list()
  for protein in proteins:
    if isTargetPeptide and protein.startswith("REV__"): 
      continue
    else:
      new_proteins.append(protein)
  return new_proteins

  
def cleanPeptide(peptide: str, removeFlanks: bool = True):
  """
  Removes modifications in [] or () as well as the _ that MaxQuant puts in front and after the sequence
  
  Warning! Dirty hack: last removal of ) is needed because some MQ versions have nested parentheses for modifications
  """
  if removeFlanks:
    peptide = peptide[1:-1]
  return re.sub(r'\[[^]]*\]', '', re.sub(r'\([^)]*\)', '', peptide)).replace(")", "") 


def chunks(lst: List, n: int):
  """
  Splits a list into sub lists of length n by taking n consecutive elements at a time
  """
  for i in range(0, len(lst), n):
    yield lst[i:i + n]

