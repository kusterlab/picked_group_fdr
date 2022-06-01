from abc import ABC, abstractmethod
from typing import List, Set, Dict, Tuple
import logging

import numpy as np

from . import parsers
from . import helpers
from .protein_groups import ProteinGroups
from .observed_peptides import ObservedPeptides
from .results import ProteinGroupResult


logger = logging.getLogger(__name__)


class ProteinGroupingStrategy(ABC):
  @abstractmethod
  def needs_peptide_to_protein_map(self):
    pass
  
  @abstractmethod
  def group_proteins(self, peptideInfoList: Dict[float, List[str]], mqProteinGroupsFile: str) -> ProteinGroups:
    pass
  
  def get_rescue_steps(self):
    return [False]
    
  def rescue_protein_groups(self, 
                            peptideInfoList: Dict[float, List[str]], 
                            proteinFdrResults: List[ProteinGroupResult], 
                            oldProteinGroups: ProteinGroups) -> ProteinGroups:
    pass

  @abstractmethod
  def short_description(self, rescue_step):
    pass
  
  @abstractmethod
  def long_description(self, rescue_step):
    pass


class NoGrouping(ProteinGroupingStrategy):
  def needs_peptide_to_protein_map(self):
    return True
  
  def group_proteins(self, peptideInfoList: Dict[float, List[str]], mqProteinGroupsFile: str) -> ProteinGroups:
    """
    Creates a single protein group for each protein in peptideInfoList
    
    :param peptideInfoList: dictionary of peptide -> (score, proteins)
    :param mqProteinGroupsFile: not used
    :returns: protein groups object
    """
    proteinGroups = ProteinGroups()
    seenProteins = set()
    for _, (_, proteins) in peptideInfoList.items():
      for protein in proteins:
        if protein not in seenProteins:
          seenProteins.add(protein)
          proteinGroups.append([protein])
    proteinGroups.create_index()
    return proteinGroups
  
  def short_description(self, rescue_step):
    return 'nG'

  def long_description(self, rescue_step):
    return 'no protein grouping'


class SubsetGrouping(ProteinGroupingStrategy):
  def needs_peptide_to_protein_map(self):
    return True
  
  def group_proteins(self, peptideInfoList: Dict[float, List[str]], mqProteinGroupsFile: str) -> ProteinGroups:
    logger.info("Grouping proteins (subset strategy):")
    
    observedPeptides = ObservedPeptides()
    observedPeptides.create(peptideInfoList)
    
    proteinGroups = observedPeptides.generate_protein_groups()
    return proteinGroups
    
  def short_description(self, rescue_step):
    return 'sG'

  def long_description(self, rescue_step):
    return 'subset protein grouping'

class MQNativeGrouping(ProteinGroupingStrategy):
  """Use grouping provided by a MaxQuant proteinGroups.txt file"""
  def needs_peptide_to_protein_map(self):
    return False
  
  def group_proteins(self, peptideInfoList: Dict[float, List[str]], mqProteinGroupsFile: str) -> ProteinGroups:
    if not mqProteinGroupsFile:
      raise ValueError("Missing MQ protein groups file input --mq_protein_groups")
    
    proteinGroups = ProteinGroups.from_mq_protein_groups_file(mqProteinGroupsFile)
    proteinGroups.create_index()
    return proteinGroups
  
  def short_description(self, rescue_step):
    return 'sG'

  def long_description(self, rescue_step):
    return 'subset protein grouping'

class RescuedGrouping:
  score_cutoff: float
  
  def rescue_protein_groups(self, 
                            peptideInfoList: Dict[float, List[str]], 
                            proteinFdrResults: List[ProteinGroupResult], 
                            oldProteinGroups: ProteinGroups) -> ProteinGroups:
    self._calculate_rescue_score_cutoff(proteinFdrResults)
    peptideInfoListFiltered = self._filter_peptide_list_by_score_cutoff(peptideInfoList)
    return self.merge_with_rescued_protein_groups(peptideInfoListFiltered, oldProteinGroups)
  
  def get_rescue_steps(self):
    return [False, True]
  
  def short_description(self, rescue_step):
    if rescue_step:
      return 'rsG'
    else:
      return 'sG'
  
  def long_description(self, rescue_step):
    if rescue_step:
      return 'rescued subset protein grouping'
    else:
      return 'subset protein grouping'
  
  def get_rescued_protein_groups(self, peptideInfoList: Dict[float, List[str]]):
    """Generates new protein grouping by "rescuing" protein groups that were split
    due to low-scoring PSMs uniquely mapping to particular isoforms
    
    :param peptideInfoList: dictionary of peptide -> (score, proteins)
    :returns: updated list of protein groups
    """
    
    logger.info("Redoing protein grouping using peptides below equivalent 1% protein FDR threshold")
    
    observedPeptides = ObservedPeptides()
    observedPeptides.create(peptideInfoList)
    
    newProteinGroups = observedPeptides.generate_protein_groups()
    
    connectedProteinGraphs = observedPeptides.get_connected_proteins(newProteinGroups)
    newProteinGroups = connectedProteinGraphs.decouple_connected_proteins(newProteinGroups)
    
    return newProteinGroups
  
  def merge_with_rescued_protein_groups(self, peptideInfoList: Dict[float, List[str]], proteinGroups: ProteinGroups):  
    newProteinGroups = self.get_rescued_protein_groups(peptideInfoList)
    newProteinGroups.add_unseen_protein_groups(proteinGroups) 
       
    logger.info(f"#protein groups before rescue: {len(proteinGroups)}")
    logger.info(f"#protein groups after rescue: {len(newProteinGroups)}")
    
    return newProteinGroups
    
  def _filter_peptide_list_by_score_cutoff(self, peptideInfoList: Dict[float, List[str]]):
    return { peptide : (score, proteins) for peptide, (score, proteins) in peptideInfoList.items() if score < self.score_cutoff }
  
  def _calculate_rescue_score_cutoff(self, proteinFdrResults: List[ProteinGroupResult]):
    """Calculate PEP corresponding to 1% protein-level FDR. N.B. this only works if the protein score is bestPEP!"""
    proteinScores, reportedQvals = zip(*[(x.score, x.qValue) for x in proteinFdrResults])
    identifiedProteinScores = [x for x, y in zip(proteinScores, reportedQvals) if y < 0.01]
    
    if len(identifiedProteinScores) == 0:
      raise ValueError("Could not calculate rescuing threshold as no proteins were found below 1% FDR")

    self.score_cutoff = np.power(10, min(identifiedProteinScores)*-1)
    logger.info(f"Rescuing threshold: protein score = {'{0:.3g}'.format(min(identifiedProteinScores))}, peptide PEP = {'{0:.3g}'.format(self.score_cutoff)}")


class RescuedSubsetGrouping(RescuedGrouping, SubsetGrouping):
  def long_description(self, rescue_step):
    if rescue_step:
      return 'rescued subset protein grouping'
    else:
      return 'subset protein grouping'


class RescuedMQNativeGrouping(RescuedGrouping, MQNativeGrouping):
  def long_description(self, rescue_step):
    if rescue_step:
      return 'rescued subset protein grouping'
    else:
      return 'subset protein grouping'


class PseudoGeneGrouping(ProteinGroupingStrategy):
  """
  Groups proteins that have at least one peptide in common, forming pseudo-genes.
  This is necessary for consistent quantification across runs in the presence of 
  isoform specific peptides only detected in a subset of the runs.
  
  For example:
  
  Run 1: 
  peptide1 -> (isoformA1, isoformA2)
  peptide2 -> (isoformA1)
  peptide3 -> (isoformA2)
  
  Run 2:
  peptide1 -> (isoformA1, isoformA2)
  
  The regular protein grouping would create 2 protein groups, one with isoformA1
  and one with isoformA2. However, since Run 2 only has shared peptides for both
  protein groups, it would report neither as detected/quantified even though
  evidence exists for the corresponding gene. By grouping both isoforms in one
  group, it can be detected/quantified in both runs.
  """
  
  def needs_peptide_to_protein_map(self):
    return True
  
  def group_proteins(self, peptideInfoList: Dict[float, List[str]], mqProteinGroupsFile: str) -> ProteinGroups:
    observedPeptides = ObservedPeptides()
    observedPeptides.create(peptideInfoList)   
    
    newProteinGroups = observedPeptides.generate_protein_groups()
    
    connectedProteinGraphs = observedPeptides.get_connected_proteins(newProteinGroups, exclude_identified_proteins = False) 
    newProteinGroups = connectedProteinGraphs.get_connected_proteins(newProteinGroups)
    
    return newProteinGroups
  
  def short_description(self, rescue_step):
    return 'gG'
  
  def long_description(self, rescue_step):
    return 'pseudo-gene grouping'
