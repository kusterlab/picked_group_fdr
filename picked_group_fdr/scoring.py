from abc import ABC, abstractmethod
from typing import Dict
import logging

import numpy as np

from . import parsers
from . import helpers
from . import fdr
from .observed_peptides import ObservedPeptides


logger = logging.getLogger(__name__)


"""
scoreType:  MQ_protein = MaxQuant protein score from proteinGroups.txt (see below)
            multPEP = emulate MaxQuant protein score with multPEP with dividing of PEPs by fixed value
            bestPEP = best scoring peptide based on PEP
            Andromeda = raw Andromeda score from evidence.txt

razor: True = use Occam's razor to decide which protein to assign to a shared peptide

useSharedPeptides:  True = use peptides shared by multiple proteins / protein groups
                    False = discard peptides shared by multiple proteins / protein groups
"""


class ProteinScore(ABC):
  @abstractmethod
  def calculate_score(self, scorePeptidePairs):
    pass
  
  @abstractmethod
  def can_do_protein_group_rescue(self):
    pass
  
  @abstractmethod
  def get_score_column(self, percolator_input):
    pass
    
  @abstractmethod
  def short_description(self):
    pass
  
  def optimize_hyperparameters(self, proteinGroups, proteinGroupScores):
    pass


class MQProteinScore(ProteinScore):
  mq_protein_groups_file: str
  
  def __init__(self, mq_protein_groups_file):
    self.mq_protein_groups_file = mq_protein_groups_file
  
  def calculate_score(self, scorePeptidePairs):
    raise NotImplementedError
  
  def can_do_protein_group_rescue(self):
    return False
  
  def get_score_column(self, percolator_input):
    return None
  
  def short_description(self):
    return 'm'
  
  def long_description(self):
    return 'multiplication of'
    
  def get_protein_scores_from_file(self):
    proteinGroupScores = list()
    for proteinGroup, proteinScore in parsers.parseMqProteinGroupsFile(self.mq_protein_groups_file):
      proteinGroupScores.append([(proteinScore, "NA", proteinGroup)])
    return proteinGroupScores


class BestAndromedaScore(ProteinScore):
  def calculate_score(self, scorePeptidePairs):
    return max([y[0] for y in scorePeptidePairs]) if len(scorePeptidePairs) > 0 else -100.0
  
  def can_do_protein_group_rescue(self):
    return False
  
  def get_score_column(self, percolator_input):
    return 'score'
  
  def short_description(self):
    return 'b'

  def long_description(self):
    return 'best'

class BestPEPScore(ProteinScore):
  def calculate_score(self, scorePeptidePairs):
    #return max([-1*y[0] for y in scorePeptidePairs]) if len(scorePeptidePairs) > 0 else -100.0
    return max([-1*np.log10(y[0] + np.nextafter(0,1)) for y in scorePeptidePairs]) if len(scorePeptidePairs) > 0 else -100.0
  
  def can_do_protein_group_rescue(self):
    return True
  
  def get_score_column(self, percolator_input):
    if percolator_input:
      return 'posterior_error_prob'
    else:
      return 'pep'
  
  def short_description(self):
    return 'b'

  def long_description(self):
    return 'best'

class MultPEPScore(ProteinScore):
  div: float
  
  def __init__(self):
    self.div = 1.0
  
  def optimize_hyperparameters(self, proteinGroups, proteinGroupScores):
    proteinScoreTuples = list()
    for proteinGroup, proteinGroupScoreList in zip(proteinGroups, proteinGroupScores):
      multPEP, numPeptides = self._get_score_and_num_peptides(proteinGroupScoreList)
      if numPeptides > 0 and not np.isnan(multPEP):
        proteinScoreTuples.append([multPEP, numPeptides, helpers.isDecoy(proteinGroup)])
    
    proteinScoreTuples = np.array(proteinScoreTuples)
    min_range, max_range = 0.1, 1.0
    for step_size in [1e-1, 1e-2, 1e-3, 1e-4]:
      numIdentifiedTargets, div = self._get_optimal_div(proteinScoreTuples, np.arange(min_range, max_range, step_size))
      min_range, max_range = div - step_size*0.9, div + step_size*0.9
    
    self.div = div
    logger.info(f"Optimal division factor for multPEP score: {self.div:.4f} (#targets = {numIdentifiedTargets})")
  
  def _get_optimal_div(self, proteinScoreTuples, div_range: np.array):  
    mostIdentifiedTargets = (-np.inf, 1.0)
    for div in div_range:
      proteinScores = proteinScoreTuples[:,0] + proteinScoreTuples[:,1] * np.log10(div)
      
      sortIdxs = np.argsort(proteinScores)[::-1]
      numDecoys = proteinScoreTuples[:,2][sortIdxs].cumsum()
      numTargets = (-proteinScoreTuples[:,2] + 1)[sortIdxs].cumsum()
      
      fdrs = np.divide(numDecoys+1, numTargets+1)
      qvals = fdr.fdrsToQvals(fdrs)
      numIdentifiedTargets = fdr.countBelowThreshold(qvals, 0.01, proteinScoreTuples[:,2][sortIdxs])
      if numIdentifiedTargets > mostIdentifiedTargets[0]:
        mostIdentifiedTargets = (numIdentifiedTargets, div)
    return mostIdentifiedTargets

  def calculate_score(self, scorePeptidePairs):
    multPEP, numPeptides = self._get_score_and_num_peptides(scorePeptidePairs)
    multPEP += np.log10(self.div) * numPeptides
    if numPeptides == 0 or np.isnan(multPEP):
      multPEP = -100.0
    return multPEP
  
  def _get_score_and_num_peptides(self, scorePeptidePairs):
    multPEP = 0.0
    seenPeptides = set()
    for PEP, peptide, _ in sorted(scorePeptidePairs):
      if peptide not in seenPeptides:
        seenPeptides.add(peptide)
        multPEP -= np.log10(PEP + np.nextafter(0,1))
    return multPEP, len(seenPeptides)
    
  def can_do_protein_group_rescue(self):
    return True
  
  def get_score_column(self, percolator_input):
    if percolator_input:
      return 'posterior_error_prob'
    else:
      return 'pep'
  
  def short_description(self):
    return 'm'
  
  def long_description(self):
    return 'multiplication of'


class ScoreOrigin(ABC):
  @abstractmethod
  def get_evidence_file(self, args):
    pass
  
  @abstractmethod
  def remaps_peptides_to_proteins(self):
    pass
  
  @abstractmethod
  def can_do_quantification(self):
    pass


class PercolatorInput(ScoreOrigin):
  def get_evidence_file(self, args):
    return args.perc_evidence
  
  def remaps_peptides_to_proteins(self):
    return False
  
  def can_do_quantification(self):
    return False
  
  def short_description(self):
    return 'p'
  
  def long_description(self):
    return 'Percolator'


class PercolatorInputRemapped(PercolatorInput):  
  def remaps_peptides_to_proteins(self):
    return True
    

class MaxQuantInput(ScoreOrigin):
  def get_evidence_file(self, args):
    return args.mq_evidence
  
  def remaps_peptides_to_proteins(self):
    return True
  
  def can_do_quantification(self):
    return True
  
  def short_description(self):
    return 'm'
  
  def long_description(self):
    return 'MaxQuant'


class ProteinScoringStrategy:
  use_proteotypicity: bool
  use_razor: bool
  use_shared_peptides: bool
  protein_score: ProteinScore
  score_origin: ScoreOrigin
  peptide_counts_per_protein: Dict[str, int]
  
  def __init__(self, score_description, mq_protein_groups_file = ""):
    if "multPEP" in score_description:
      self.protein_score = MultPEPScore()
    elif "bestPEP" in score_description:
      self.protein_score = BestPEPScore()
    elif "Andromeda" in score_description:
      self.protein_score = BestAndromedaScore()
    elif "MQ_protein" in score_description:
      self.protein_score = MQProteinScore(mq_protein_groups_file)
    else:
      raise NotImplementedError
    
    self.use_razor = "razor" in score_description
    self.use_proteotypicity = "proteotypicity" in score_description
    self.use_shared_peptides = "with_shared" in score_description
    
    if "Perc" in score_description:
      # Add "remap" to the scoreType if the fasta database used for protein grouping is different from the one used by Percolator
      if "remap" in score_description:
        self.score_origin = PercolatorInputRemapped()
      else:
        self.score_origin = PercolatorInput()
    else:
      self.score_origin = MaxQuantInput()
    
  def get_evidence_file(self, args):
    return self.score_origin.get_evidence_file(args)
  
  def remaps_peptides_to_proteins(self):
    return self.score_origin.remaps_peptides_to_proteins()
  
  def can_do_quantification(self):
    return self.score_origin.can_do_quantification()
  
  def optimize_hyperparameters(self, proteinGroups, proteinGroupScores):
    return self.protein_score.optimize_hyperparameters(proteinGroups, proteinGroupScores)
  
  def calculate_score(self, scorePeptidePairs):
    return self.protein_score.calculate_score(scorePeptidePairs)
  
  def get_score_column(self):
    return self.protein_score.get_score_column(self.score_origin.short_description() == 'p')
  
  def can_do_protein_group_rescue(self):
    return self.protein_score.can_do_protein_group_rescue()
  
  def short_description(self):
    return self.protein_score.short_description() + self.score_origin.short_description() + 'P'
  
  def short_description_razor(self):
    return "rS" if self.use_razor else "dS"
  
  def long_description(self):
    return f"{self.protein_score.long_description()} {self.score_origin.long_description()} PEP"
  
  def long_description_razor(self):
    return "razor peptides" if self.use_razor else "discard shared peptides"
    
  def filter_proteins(self, proteins):
    if self.use_razor:
      return self._retain_protein_with_most_observed_peptides(proteins)
    else:
      return proteins
  
  def set_peptide_counts_per_protein(self, peptideInfoList):
    if self.use_razor:
      observedPeptides = ObservedPeptides()
      observedPeptides.create(peptideInfoList)
      self.peptide_counts_per_protein = observedPeptides.get_peptide_counts_per_protein()
  
  def _retain_protein_with_most_observed_peptides(self, proteins):
    """Retains only the protein with the most observed peptides, ties are broken based on alphabetical order"""
    numPeptidesPerProteinPairs = [(self.peptide_counts_per_protein.get(protein, 0), protein) for protein in proteins]
    pairWithMostObservedPeptides = sorted(numPeptidesPerProteinPairs, reverse = True)[0]
    proteinWithMostObservedPeptides = pairWithMostObservedPeptides[1]
    return [proteinWithMostObservedPeptides]
  
  def collect_peptide_scores_per_protein(self,
      proteinGroups, peptideInfoList, suppressMissingProteinWarning = False):
    """Groups peptides with associated scores by protein
    
    :param proteinGroups: ProteinGroups object
    :param peptideInfoList: Dict of peptide -> (score, proteins)
    :param suppressMissingProteinWarning: suppresses the warning for missing proteins
      in the proteinGroups object. This is set during the rescuing grouping procedure
      since some protein groups will have been filtered out in the rescuing step.
    :returns: lists of (score, peptide, proteins) tuples per protein group
    """
    if not self.get_score_column():
      return self.protein_score.get_protein_scores_from_file()
    
    logger.info("Assigning peptides to protein groups")
    sharedPeptides, uniquePeptides = 0, 0
    proteinGroupScores = [list() for _ in range(len(proteinGroups))]
    
    for peptide, (score, proteins) in peptideInfoList.items():
      proteins = self.filter_proteins(proteins) # filtering for razor peptide approach
      
      proteinGroupIdxs = proteinGroups.get_protein_group_idxs(proteins)
      if len(proteinGroupIdxs) == 0 and not suppressMissingProteinWarning:
        raise Exception(f"Could not find any of the proteins {proteins} in the ProteinGroups object, check if the identifier format is the same. 1st protein group in ProteinGroups object: {proteinGroups.protein_groups[0]}")
      
      if not self.use_shared_peptides and len(proteinGroupIdxs) > 1: # ignore shared peptides
        sharedPeptides += 1
      else:
        uniquePeptides += 1
        for proteinGroupIdx in proteinGroupIdxs:
          proteinGroupScores[proteinGroupIdx].append((score, peptide, proteins))
    
    logger.info(f"Shared peptides: {sharedPeptides}; Unique peptides: {uniquePeptides}")
    return proteinGroupScores
  

def compareRazorPeptides(mqEvidenceFile, peptideToProteinMap, proteinGroups, scoreType):
  """Compares the chosen protein by MaxQuant according to the razor peptide rule with our implementation of the razor peptide rule"""
  scoreType = ProteinScoringStrategy("multPEP razor")
  for peptideRow in parsers.parseMqEvidenceFile(mqEvidenceFile, scoreType = scoreType):
    peptide, tmp_proteins, _, _ = peptideRow
    
    proteins = digest.getProteins(peptideToProteinMap, helpers.cleanPeptide(peptide))
    
    leadingProteins = proteinGroups.get_leading_proteins(proteins)
    if len(leadingProteins) > 1:
      predictedRazor = scoreType.filter_proteins(proteins) # filtering for razor peptide approach
      logger.debug(f"{tmp_proteins[0] == predictedRazor[0]} {tmp_proteins[0]} {predictedRazor[0]}")
