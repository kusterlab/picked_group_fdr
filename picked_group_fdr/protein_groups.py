from typing import List, Set, Dict
import logging

from . import parsers
from . import helpers
from . import entrapment


logger = logging.getLogger(__name__)


class ProteinGroups:
  protein_groups: List[List[str]]
  protein_to_group_idx_map: Dict[str, int]
  valid_idx: bool
  
  def __init__(self, protein_groups: List[List[str]] = None):
    if protein_groups is None: # https://docs.python-guide.org/writing/gotchas/
      self.protein_groups = []
    else:
      self.protein_groups = protein_groups
    self.protein_to_group_idx_map = dict()
    self.valid_idx = False
  
  @classmethod
  def from_mq_protein_groups_file(cls, mqProteinGroupsFile: str):
    protein_groups = [proteinGroup for proteinGroup, _ in parsers.parseMqProteinGroupsFile(mqProteinGroupsFile)]
    return cls(protein_groups)
  
  @classmethod
  def from_observed_peptide_map(cls, protein_to_peptides_dict):
    protein_groups = [[protein] for protein in protein_to_peptides_dict.keys()]
    return cls(protein_groups)
    
  @classmethod
  def from_protein_group_results(cls, proteinGroupsResults):
    protein_groups = [pgr.proteinIds.split(";") for pgr in proteinGroupsResults]
    return cls(protein_groups)
  
  def __iter__(self):
    return iter(self.protein_groups)
  
  def __next__(self):
    return next(self.protein_groups)
  
  def __len__(self):
    return len(self.protein_groups)
  
  def __eq__(self, other):
    return sorted(self.protein_groups) == sorted(other.protein_groups)
  
  def append(self, protein_group: List[str]):
    self.protein_groups.append(protein_group)
    self.valid_idx = False
  
  def size(self) -> int:
    return len(self.protein_groups)
  
  def get_all_proteins(self):
    return set([protein for protein_group in self.protein_groups for protein in protein_group])
  
  def get_leading_proteins(self, proteins: List[str]) -> Set[str]:
    """Returns the leading protein for each of the protein groups that the proteins belong to"""
    getLeadingProtein = lambda protein : self.get_protein_group(protein)[0]
    return set(map(getLeadingProtein, proteins))
  
  def create_index(self):
    """ Creates a map of proteins to their index in the proteinGroups list

    :param proteinGroups: list of protein groups, where each entry is itself a 
      list of proteins in that protein group
    :returns: dictionary of protein -> index in proteinGroups list
    """
    self.protein_to_group_idx_map = dict()
    for proteinGroupIdx, proteinGroup in enumerate(self.protein_groups):
      for protein in proteinGroup:
        self.protein_to_group_idx_map[protein] = proteinGroupIdx
    self.valid_idx = True
  
  def get_protein_group(self, protein: str, check_idx_valid: bool = True):
    return self.protein_groups[self._get_protein_group_idx(protein, check_idx_valid)]
  
  def _get_protein_group_idx(self, protein: str, check_idx_valid: bool = True):
    if not self.valid_idx and check_idx_valid:
      raise Exception("Trying to get group index while index is invalid")
    return self.protein_to_group_idx_map[protein]
  
  def get_protein_group_idxs(self, proteins: List[str], check_idx_valid: bool = True) -> Set[int]:
    if not self.valid_idx and check_idx_valid:
      raise Exception("Trying to get group index while index is invalid")
    
    proteinGroupIdxs = set()
    for protein in proteins:  
      if protein in self.protein_to_group_idx_map:
        proteinGroupIdxs.add(self._get_protein_group_idx(protein, check_idx_valid))
    return proteinGroupIdxs
  
  def merge_groups(self, superSetProtein: str, protein: str):
    """Merge two protein groups"""
    superSetGroupIdx = self._get_protein_group_idx(superSetProtein, check_idx_valid = False)
    self.protein_groups[superSetGroupIdx].extend(self.get_protein_group(protein, check_idx_valid = False))
    self.protein_groups[self._get_protein_group_idx(protein, check_idx_valid = False)] = []
    self.valid_idx = False
  
  def add_unseen_protein_groups(self, proteinGroups):
    """Add protein groups from another ProteinGroups object that do not have
    proteins in common with a protein group in the current ProteinGroups object"""
    identifiedNewProteins = self.get_all_proteins()
    for proteinGroup in proteinGroups:
      newProteinGroup = list(filter(lambda protein : protein not in identifiedNewProteins, proteinGroup))
      if len(newProteinGroup) > 0:
        self.append(newProteinGroup)
    self.create_index()
    
  def remove_empty_groups(self):
    self.protein_groups = list(filter(lambda x : len(x) > 0, self.protein_groups))
    self.create_index()
  
  def print(self, minSize: int = 0):  
    for pg in self.protein_groups:
      if len(pg) >= minSize:
        logger.info(f'{";".join(pg)} {"entrapmentGroup" if entrapment.isEntrapment(pg) else ""} {"decoyGroup" if helpers.isDecoy(pg) else ""}')

