from typing import List, Dict
import logging

import numpy as np

from .. import helpers
from .base import ProteinGroupColumns


logger = logging.getLogger(__name__)


class SequenceCoverageColumns(ProteinGroupColumns):
  proteinSequences: Dict[str, str]
  
  def __init__(self, proteinSequences):
    self.proteinSequences = proteinSequences
  
  def append_headers(self, proteinGroupResults, experiments):
    proteinGroupResults.append_header('Sequence coverage [%]')
    proteinGroupResults.append_header('Unique + razor sequence coverage [%]')
    proteinGroupResults.append_header('Unique sequence coverage [%]')
    # TODO: add these columns
    # Mol. weight [kDa]
    # Sequence length
    # Sequence lengths
    
    for experiment in experiments:
      proteinGroupResults.append_header('Sequence coverage [%] ' + experiment)
  
  def append_columns(self, proteinGroupResults, experimentToIdxMap, postErrProbCutoff):
    logger.info("Doing quantification: sequence coverage")
    for pgr in proteinGroupResults:
      sequenceCoverages = self.get_sequence_coverages(pgr.precursorQuants, experimentToIdxMap, postErrProbCutoff, pgr.proteinIds)
      pgr.extend(sequenceCoverages)
  
  def get_sequence_coverages(self, peptideIntensityList, experimentToIdxMap, postErrProbCutoff, proteinIds):
    peptidesPerExperiment = self.unique_peptides_per_experiment(peptideIntensityList, experimentToIdxMap, postErrProbCutoff)
    return self.calculate_sequence_coverages(peptidesPerExperiment, proteinIds)
      
  def unique_peptides_per_experiment(self, peptideIntensityList, experimentToIdxMap, postErrProbCutoff):
    uniquePeptides = [set() for _ in range(len(experimentToIdxMap))]
    for precursor in peptideIntensityList:
      if helpers.isMbr(precursor.postErrProb) or precursor.postErrProb <= postErrProbCutoff:
        uniquePeptides[experimentToIdxMap[precursor.experiment]].add(helpers.cleanPeptide(precursor.peptide))
    return uniquePeptides

  def calculate_sequence_coverages(self, peptidesPerExperiment, proteinIds):
    proteinId = proteinIds.split(";")[0]
    proteinSequence = self.proteinSequences.get(proteinId, "")
    coverageTotal = np.zeros(len(proteinSequence))
    coverageExperimentRatios = list()
    for peptides in peptidesPerExperiment:
      if len(peptides) == 0:
        coverageExperimentRatios.append(0.0)
        continue
      
      coverageExperiment = np.zeros(len(proteinSequence))
      for cleanPeptide in peptides:
        pos = proteinSequence.find(cleanPeptide)
        coverageTotal[pos:pos+len(cleanPeptide)] = 1
        coverageExperiment[pos:pos+len(cleanPeptide)] = 1
      
      coverageExperimentRatios.append(self.calculate_coverage_ratio(coverageExperiment))
      
    coverageTotalRatio = self.calculate_coverage_ratio(coverageTotal)
    allCoverageRatios = [coverageTotalRatio, 
                         coverageTotalRatio, 
                         coverageTotalRatio] + coverageExperimentRatios
    return self.format_as_percentage(allCoverageRatios)
  
  @staticmethod
  def format_as_percentage(l):
    return ['%.1f' % (x*100) for x in l]
  
  @staticmethod
  def calculate_coverage_ratio(coverage):
    if len(coverage) > 0:
      return coverage.sum() / len(coverage)
    else:
      return 0.0
