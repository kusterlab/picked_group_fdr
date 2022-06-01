from typing import List, Dict
import logging

from .. import helpers
from .base import ProteinGroupColumns


logger = logging.getLogger(__name__)


class UniquePeptideCountColumns(ProteinGroupColumns):
  def append_headers(self, proteinGroupResults, experiments):
    for experiment in experiments:
      proteinGroupResults.append_header('Unique peptides ' + experiment)
  
  def append_columns(self, proteinGroupResults, experimentToIdxMap, postErrProbCutoff):
    logger.info("Doing quantification: Count unique peptides")
    for pgr in proteinGroupResults:
      pepCounts = self.uniquePeptideCountsPerExperiment(pgr.precursorQuants, experimentToIdxMap, postErrProbCutoff)
      pgr.extend(pepCounts)
  
  @staticmethod
  def uniquePeptideCountsPerExperiment(peptideIntensityList, experimentToIdxMap, postErrProbCutoff):
    uniquePeptides = [set() for _ in range(len(experimentToIdxMap))]
    for precursor in peptideIntensityList:
      if helpers.isMbr(precursor.postErrProb) or precursor.postErrProb <= postErrProbCutoff:
        uniquePeptides[experimentToIdxMap[precursor.experiment]].add(precursor.peptide)
    return list(map(len, uniquePeptides))
