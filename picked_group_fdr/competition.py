from abc import ABC, abstractmethod
from typing import List, Tuple

import numpy as np

from .protein_groups import ProteinGroups
from . import helpers

class ProteinCompetitionStrategy(ABC):
  @abstractmethod
  def _add_seen_proteins(self, proteins: List[str]) -> None:
    pass
  
  @abstractmethod
  def _is_protein_seen(self, proteins: List[str]) -> bool:
    pass

  @abstractmethod
  def short_description(self) -> str:
    pass
  
  @abstractmethod
  def long_description(self) -> str:
    pass
    
  @abstractmethod
  def reset(self) -> None:
    pass
    
  def do_competition(self, proteinGroups, proteinGroupScores, scoreType) -> Tuple[ProteinGroups, List]:
    scoreGroupPairs = list(zip(proteinGroupScores, proteinGroups))
    
    # shuffle the list of proteinGroups before sorting by protein score. This avoids 
    # biases between groups with equal score, e.g. when all target protein groups 
    # are listed above decoy protein groups.
    np.random.shuffle(scoreGroupPairs)
    
    scoreGroupPairs = sorted(scoreGroupPairs, key = lambda x : scoreType.calculate_score(x[0]), reverse = True)
    scoreGroupPairs = filter(lambda x : len(x[0]) > 0, scoreGroupPairs)
    proteinGroupScores, proteinGroups = zip(*scoreGroupPairs)
    
    filteredProteinGroups = ProteinGroups()
    filteredProteinGroupScores = []
    for proteinGroup, proteinGroupScoreList in zip(proteinGroups, proteinGroupScores):
      if self._is_protein_seen(proteinGroup) or helpers.isContaminant(proteinGroup):
        continue
      
      self._add_seen_proteins(proteinGroup)
      filteredProteinGroups.append(proteinGroup)
      filteredProteinGroupScores.append(proteinGroupScoreList)
    
    self.reset() # clears list of seen proteins
    
    return filteredProteinGroups, filteredProteinGroupScores


class PickedStrategy(ProteinCompetitionStrategy):
  seenProteins: set
  
  def __init__(self):
    self.seenProteins = set()
  
  def _add_seen_proteins(self, proteins: List[str]) -> None:
    self.seenProteins.add(self._get_protein_group_string(proteins))

  def _is_protein_seen(self, proteins: List[str]) -> bool:
    return self._get_protein_group_string(proteins) in self.seenProteins
  
  def _get_protein_group_string(self, proteins: List[str]) -> str:
    return ";".join(map(lambda x : x.replace("REV__",""), proteins))

  def short_description(self) -> str:
    return "pT"
  
  def long_description(self) -> str:
    return "picked target-decoy strategy"
  
  def reset(self) -> None:
    self.seenProteins = set()


class PickedGroupStrategy(ProteinCompetitionStrategy):
  seenProteins: set
  
  def __init__(self):
    self.seenProteins = set()
  
  def _add_seen_proteins(self, proteins: List[str]) -> None:
    for protein in proteins:
      self.seenProteins.add(protein.replace("REV__",""))

  def _is_protein_seen(self, proteins: List[str]) -> bool:
    return sum([1 for x in proteins if x.replace("REV__","") in self.seenProteins]) > 0
  
  def short_description(self) -> str:
    return "pgT"
  
  def long_description(self) -> str:
    return "picked group target-decoy strategy"
  
  def reset(self) -> None:
    self.seenProteins = set()


class ClassicStrategy(ProteinCompetitionStrategy):
  def _add_seen_proteins(self, proteins: List[str]) -> None:
    pass

  def _is_protein_seen(self, proteins: List[str]) -> bool:
    return False
  
  def short_description(self) -> str:
    return "cT"
  
  def long_description(self) -> str:
    return "classic target-decoy strategy"
  
  def reset(self) -> None:
    pass
  
