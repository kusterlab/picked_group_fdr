from typing import List, Dict, Set
import logging

from . import parsers


logger = logging.getLogger(__name__)


def markEntrapmentProteins(peptideToProteinMap: Dict, mqProteinGroupsFile: str) -> None:
  """
  necessary for simulated datasets with entrapments 
  since the peptide-protein map is made without knowledge of which proteins 
  are entrapment proteins we update the peptideToProteinMap by adding a
  "_entrapment" suffix to the simulated entrapment proteins
  """
  if not mqProteinGroupsFile:
    return
  
  entrapmentProteins = getEntrapmentProteins(mqProteinGroupsFile)
  if len(entrapmentProteins) > 0:
    logger.info("Updating peptide to protein map with entrapment protein identifiers")
    for peptide, proteins in peptideToProteinMap.items():
      peptideToProteinMap[peptide] = markEntrapment(proteins, entrapmentProteins)


def getEntrapmentProteins(mqProteinGroupsFile: str) -> Set[str]:
  entrapmentProteins = set()
  for proteinGroup, proteinScore in parsers.parseMqProteinGroupsFile(mqProteinGroupsFile):
    for protein in proteinGroup:
      if "_entrapment" in protein:
        entrapmentProteins.add(protein)
  return entrapmentProteins


def markEntrapment(proteins: List[str], entrapmentProteins: Set[str]) -> List[str]:
  return list(map(lambda protein : protein + "_entrapment" if protein + "_entrapment" in entrapmentProteins else protein, proteins))


def isEntrapment(proteinGroup: List[str]) -> bool:
  return len([1 for x in proteinGroup if "_entrapment" in x or "Random_" in x or "REV__" in x]) == len(proteinGroup)
  #return len([1 for x in proteinGroup if "Random_" in x or "REV__" in x]) == len(proteinGroup)

