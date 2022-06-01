from typing import List, Dict
import logging

import numpy as np

from .. import helpers
from .base import ProteinGroupColumns


logger = logging.getLogger(__name__)


class TMTIntensityColumns(ProteinGroupColumns):
  numTmtChannels: int
  
  def __init__(self, numTmtChannels):
    self.numTmtChannels = numTmtChannels
  
  def append_headers(self, proteinGroupResults, experiments):
    for experiment in experiments:
      for i in range(1, self.numTmtChannels+1):
        proteinGroupResults.append_header('Reporter intensity corrected ' + str(i) + " " + experiment)
      
      for i in range(1, self.numTmtChannels+1):
        proteinGroupResults.append_header('Reporter intensity ' + str(i) + " " + experiment)
      
      for i in range(1, self.numTmtChannels+1):
        proteinGroupResults.append_header('Reporter intensity count ' + str(i) + " " + experiment)
  
  def append_columns(self, proteinGroupResults, experimentToIdxMap, postErrProbCutoff):
    logger.info("Doing quantification: TMT intensity")
    for pgr in proteinGroupResults:
      intensities = self._getTmtIntensities(pgr.precursorQuants, experimentToIdxMap, postErrProbCutoff)
      pgr.extend(intensities)
  
  def _getTmtIntensities(self, peptideIntensityList, experimentToIdxMap, postErrProbCutoff):
    intensities = [np.zeros(self.numTmtChannels*3) for i in range(len(experimentToIdxMap))]
    for precursor in peptideIntensityList:
      if helpers.isMbr(precursor.postErrProb) or precursor.postErrProb <= postErrProbCutoff:
        intensities[experimentToIdxMap[precursor.experiment]] += precursor.tmtIntensities
    return np.concatenate(intensities)

