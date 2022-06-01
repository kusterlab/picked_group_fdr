from typing import List, Dict
import logging

from .. import helpers
from .base import ProteinGroupColumns


logger = logging.getLogger(__name__)


class IdentificationTypeColumns(ProteinGroupColumns):
  def append_headers(self, proteinGroupResults, experiments):
    for experiment in experiments:
      proteinGroupResults.append_header('Identification type ' + experiment)
  
  def append_columns(self, proteinGroupResults, experimentToIdxMap, postErrProbCutoff):
    logger.info("Doing quantification: adding identification types")
    for pgr in proteinGroupResults:
      idType = self._identificationTypePerExperiment(pgr.precursorQuants, experimentToIdxMap, postErrProbCutoff)
      pgr.extend(idType)
  
  def _identificationTypePerExperiment(self, peptideIntensityList, experimentToIdxMap, postErrProbCutoff):
    idType = [""]*len(experimentToIdxMap)
    for precursor in peptideIntensityList:
      idx = experimentToIdxMap[precursor.experiment]
      if helpers.isMbr(precursor.postErrProb) and idType[idx] != "By MS/MS":
        idType[idx] = "By matching"
      elif precursor.postErrProb <= postErrProbCutoff:
        idType[idx] = "By MS/MS"
    return idType
