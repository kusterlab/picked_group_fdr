import collections
import copy
import logging
from dataclasses import dataclass, field
from typing import List

import numpy as np

from . import parsers
from . import helpers


logger = logging.getLogger(__name__)


@dataclass
class ProteinGroupResult:
  proteinIds: str = ""
  majorityProteinIds: str = ""
  peptideCountsUnique: int = 0
  proteinNames: str = ""
  geneNames: str = ""
  fastaHeaders: str = ""
  bestPeptide: str = ""
  numberOfProteins: int = 0
  qValue: float = np.nan
  score: float = np.nan
  reverse: str = ""
  potentialContaminant: str = ""
  precursorQuants: List[float] = field(default_factory=list)
  extraColumns: List[float] = field(default_factory=list)
  
  def extend(self, e):
    self.extraColumns.extend(e)
  
  def append(self, a):
    self.extraColumns.append(a)
  
  @classmethod
  def from_mq_protein_groups(cls, row, cols):
    _get_field = lambda x : row[cols[x]] if x in cols else ""
    
    proteinIds = _get_field('Protein IDs')
    majorityProteinIds = _get_field('Majority protein IDs')
    peptideCountsUnique = _get_field("Peptide counts (unique)")
    proteinNames = _get_field('Protein names')
    geneNames = _get_field('Gene names')
    fastaHeaders = _get_field('Fasta headers')
    bestPeptide = ""
    numberOfProteins = int(_get_field("Number of proteins"))
    qValue = float(_get_field('Q-value'))
    score = float(_get_field('Score'))
    reverse = _get_field('Reverse')
    potentialContaminant = _get_field('Potential contaminant')
    return cls(proteinIds, majorityProteinIds, peptideCountsUnique, 
               proteinNames, geneNames, fastaHeaders, bestPeptide, 
               numberOfProteins, qValue, score, reverse, potentialContaminant, [], [])
    
  @classmethod
  def from_protein_group(cls, proteinGroup, peptideScores, reportedFdr, proteinScore, scoreCutoff, proteinAnnotations):
    proteinIds = ";".join(proteinGroup)
    bestPeptide = sorted([(p[0], p[1]) for p in peptideScores])[0][1]
    numUniquePeptidesPerProtein = cls._get_peptide_counts(peptideScores, scoreCutoff)
    peptideCountsUnique = ";".join([str(numUniquePeptidesPerProtein[p]) for p in proteinGroup])
    majorityProteinIds = ";".join([p for p in proteinGroup if numUniquePeptidesPerProtein[p] >= max(numUniquePeptidesPerProtein.values()) / 2])
    numberOfProteins = len(proteinGroup)
    
    proteinNames, geneNames, fastaHeaders = list(), list(), list()
    for p in proteinGroup:
      if p in proteinAnnotations:
        proteinName, geneName, fastaHeader = proteinAnnotations[p]
        if proteinName not in proteinNames:
          proteinNames.append(proteinName)
        if geneName not in geneNames:
          geneNames.append(geneName)
        if fastaHeader not in fastaHeaders:
          fastaHeaders.append(fastaHeader)
    
    proteinNames = ";".join(proteinNames)
    geneNames = ";".join(geneNames)
    fastaHeaders = ";".join(fastaHeaders)
    
    qValue = reportedFdr
    score = proteinScore
    reverse = '+' if helpers.isDecoy(proteinGroup) else ''
    potentialContaminant = '+' if helpers.isContaminant(proteinGroup) else ''
    return cls(proteinIds, majorityProteinIds, peptideCountsUnique, 
               proteinNames, geneNames, fastaHeaders, bestPeptide, 
               numberOfProteins, qValue, score, reverse, potentialContaminant, [], [])
  
  @staticmethod
  def _get_peptide_counts(scorePeptidePairs, scoreCutoff):
    proteinPeptideCount = collections.defaultdict(int)
    seenPeptides = set()
    for PEP, peptide, proteins in sorted(scorePeptidePairs):
      if PEP > scoreCutoff:
        break
      if peptide not in seenPeptides:
        seenPeptides.add(peptide)
        for protein in proteins:
          proteinPeptideCount[protein] += 1
    return proteinPeptideCount
    
  def to_list(self):
    return [self.proteinIds, 
            self.majorityProteinIds, 
            self.peptideCountsUnique,
            self.proteinNames,
            self.geneNames,
            self.fastaHeaders,
            self.bestPeptide,
            self.numberOfProteins,
            self.qValue,
            self.score,
            self.reverse,
            self.potentialContaminant] + ["%.0f" % (x) if not type(x) == str else x for x in self.extraColumns]


ProteinGroupHeaders = [
      "Protein IDs", 
      "Majority protein IDs", 
      "Peptide counts (unique)", 
      "Protein names", 
      "Gene names", 
      "Fasta headers", 
      "Best peptide", 
      "Number of proteins", 
      "Q-value", 
      "Score", 
      "Reverse", 
      "Potential contaminant"]


class ProteinGroupResults:
  headers = List[str]
  protein_group_results: List[ProteinGroupResult] = field(default_factory=list)
  
  def __init__(self, protein_group_results=None):
    if protein_group_results is None: # https://docs.python-guide.org/writing/gotchas/
      self.protein_group_results = []
    else:  
      self.protein_group_results = protein_group_results
    self.headers = copy.deepcopy(ProteinGroupHeaders)
    
  def __iter__(self):
    return iter(self.protein_group_results)
  
  def __next__(self):
    return next(self.protein_group_results)
  
  def __len__(self):
    return len(self.protein_group_results)
  
  def __getitem__(self, indices):
    return self.protein_group_results[indices]
  
  def append_header(self, header):
    self.headers.append(header)
  
  def write(self, output_file):
    writer = parsers.getTsvWriter(output_file)
    writer.writerow(self.headers)
    for proteinRow in self.protein_group_results:
      writer.writerow(proteinRow.to_list())

  @classmethod
  def from_mq_protein_proteins_file(cls, mqProteinGroupsFile):
    if mqProteinGroupsFile.endswith('.csv'):
      delimiter = ','
    else:
      delimiter = '\t'
      
    reader = parsers.getTsvReader(mqProteinGroupsFile, delimiter)
    headers = next(reader)
    
    cols = { x : headers.index(x) for x in ProteinGroupHeaders if x in headers }
    
    logger.info("Parsing MaxQuant proteinGroups.txt file")
    protein_group_results = list()
    for row in reader:
      protein_group_results.append(ProteinGroupResult.from_mq_protein_groups(row, cols))
    return cls(protein_group_results)
  
  @classmethod
  def from_protein_groups(cls, proteinGroups, proteinGroupScores, proteinScores, reportedQvals, scoreCutoff, proteinAnnotations):
    protein_group_results = list()
    for proteinGroup, peptideScores, proteinScore, reportedFdr in zip(proteinGroups, proteinGroupScores, proteinScores, reportedQvals):
      pgr = ProteinGroupResult.from_protein_group(proteinGroup, peptideScores, reportedFdr, proteinScore, scoreCutoff, proteinAnnotations)
      protein_group_results.append(pgr)
    
    return cls(protein_group_results)

