from typing import List, Dict
import logging

from .. import helpers
from .base import ProteinGroupColumns


logger = logging.getLogger(__name__)


class EvidenceIdsColumns(ProteinGroupColumns):
  def append_headers(self, proteinGroupResults, experiments):
    proteinGroupResults.append_header('Evidence IDs')
  
  def append_columns(self, proteinGroupResults, experimentToIdxMap, postErrProbCutoff):
    logger.info("Doing quantification: adding evidence ids")
    for pgr in proteinGroupResults:
      evidenceIds = self.collect_evidence_ids(pgr.precursorQuants, postErrProbCutoff)
      pgr.append(";".join(map(str, evidenceIds)))
  
  @staticmethod
  def collect_evidence_ids(peptideIntensityList, postErrProbCutoff):
    evidenceIds = list()
    for precursor in peptideIntensityList:
      if helpers.isMbr(precursor.postErrProb) or precursor.postErrProb <= postErrProbCutoff:
        evidenceIds.append(precursor.evidenceId)
    return sorted(evidenceIds)


