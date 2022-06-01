import collections
import logging

from . import helpers
from . import entrapment
from . import parsers
from .results import ProteinGroupResult


logger = logging.getLogger(__name__)


def getPeptideToProteotypicityFromFile(peptideToProteotypicityFile):
  logger.info("Reading peptide to proteotypicity map")
  peptideToProteotypicityMap = dict()
  proteinToMaxProteotypicityMap = collections.defaultdict(lambda : -1000)
  proteinToNumInSilicoPeptidesMap = collections.defaultdict(int)
  reader = parsers.getTsvReader(peptideToProteotypicityFile)
  next(reader)
  for i, row in enumerate(reader):
    if (i+1) % 1000000 == 0:
      logger.info(f"Processing peptide {i+1}")
    
    #peptide, proteins, proteotypicity = row[0], row[1].split(";"), float(row[2]) + 10.0
    peptide, proteins, proteotypicity = row[0], row[4].split(";"), float(row[3]) + 1.0
    peptideToProteotypicityMap[peptide] = proteotypicity
    for protein in proteins:
      proteinToMaxProteotypicityMap[protein] = max([proteotypicity, proteinToMaxProteotypicityMap[protein]])
      proteinToNumInSilicoPeptidesMap[protein] += 1
  
  return (peptideToProteotypicityMap, proteinToMaxProteotypicityMap, proteinToNumInSilicoPeptidesMap)


def getProteotypicityScore(peptideToProteotypicityMap, peptideScorePairs, proteins, scoreCutoff):
  peptideToProteotypicityMap, proteinToMaxProteotypicityMap, proteinToNumInSilicoPeptidesMap = peptideToProteotypicityMap
  
  maxProteinProteotypicity = max(proteinToMaxProteotypicityMap[p] for p in proteins)
  numInsilicoPeptides = max(proteinToNumInSilicoPeptidesMap[p] for p in proteins)
  observedProteotypicities = [peptideToProteotypicityMap.get(p, 0.0) for s, p, _ in peptideScorePairs if s < scoreCutoff]
  if len(observedProteotypicities) > 0:
    maxObservedProteotypicity = max(observedProteotypicities)
  else:
    maxObservedProteotypicity = 0.0
  
  #if maxObservedProteotypicity == 0.0:
  #  logger.info(peptideScorePairs)
  return maxObservedProteotypicity, maxProteinProteotypicity, numInsilicoPeptides


def getProteotypicities(peptideToProteotypicityMap, peptideScorePairs, scoreCutoff):
  peptideToProteotypicityMap, proteinToMaxProteotypicityMap, proteinToNumInSilicoPeptidesMap = peptideToProteotypicityMap
  
  uniquePeptides = set([p for s, p, _ in peptideScorePairs if s < scoreCutoff])
  
  return sorted([(peptideToProteotypicityMap.get(peptide, 0.0), peptide) for peptide in uniquePeptides], reverse = True)


def calculateProteotypicityScores(proteinGroups, proteinGroupScores, peptideToProteotypicityMap, scoreType, scoreCutoff):
  logger.info("Computing proteotypicity scores")
  
  ptWriter = parsers.getTsvWriter("protein_percolator_input.tab")
  ptWriter.writerow(["PSMid", "label", "bestPEP", "proteotypicityRatio", "numObservedPeptides", "numInSilicoPeptides", "maxObservedProteotypicity", "maxProteinProteotypicity", "proteotypicityDiff",  "peptide", "proteins"])
  for proteinGroup, proteinGroupScoreList in zip(proteinGroups, proteinGroupScores):
    proteinScore = scoreType.calculate_score(proteinGroupScoreList)
    if proteinScore == -100.0:
      break
    
    maxObservedProteotypicity, maxProteinProteotypicity, numInsilicoPeptides = getProteotypicityScore(peptideToProteotypicityMap, proteinGroupScoreList, proteinGroup, scoreCutoff)
    proteotypicities = getProteotypicities(peptideToProteotypicityMap, proteinGroupScoreList, scoreCutoff)
    
    proteotypicityRatio = maxObservedProteotypicity / maxProteinProteotypicity
    proteotypicityDiff = maxProteinProteotypicity - maxObservedProteotypicity
    
    numUniquePeptides = ProteinGroupResult._get_peptide_counts(proteinGroupScoreList, scoreCutoff)
    if len(numUniquePeptides) > 0:
      numUniquePeptides = max(numUniquePeptides.values())
    else:
      numUniquePeptides = 0
    
    targetDecoyLabel = -1 if helpers.isDecoy(proteinGroup) else 1
    ptWriter.writerow([";".join(proteinGroup), targetDecoyLabel, # metadata
                       proteinScore, proteotypicityRatio, numUniquePeptides, numInsilicoPeptides, maxObservedProteotypicity, maxProteinProteotypicity, proteotypicityDiff, # features
                       "-.APEPTIDE.-", entrapment.isEntrapment(proteinGroup)] + proteotypicities) # fill peptide and proteins columns with useful metadata

