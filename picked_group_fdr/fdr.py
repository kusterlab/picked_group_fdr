import logging
from typing import List, Dict, Optional

import numpy as np

from . import helpers
from . import entrapment


logger = logging.getLogger(__name__)


def calculateProteinFDRs(proteinGroups, proteinGroupScores, scoreType):
  logger.info("Calculating protein group-level FDRs")
  numDecoys, numEntrapments, numTargets = 0, 0, 0
  proteinGroupInfoList = list()
  for proteinGroup, proteinGroupScoreList in zip(proteinGroups, proteinGroupScores):
    proteinScore = scoreType.calculate_score(proteinGroupScoreList)
    if proteinScore == -100.0:
      break
    
    if helpers.isDecoy(proteinGroup):
      numDecoys += 1
    else:
      numTargets += 1
      if entrapment.isEntrapment(proteinGroup):
        numEntrapments += 1
    reportedFdr = (numDecoys + 1) / (numTargets + 1)
    observedFdr = (numEntrapments + 1) / (numTargets + 1)
    proteinGroupInfoList.append((reportedFdr, observedFdr, helpers.isDecoy(proteinGroup), proteinScore))
    
  logger.info(f"Decoys: {numDecoys}, Entrapments: {numEntrapments}, Pool: {numTargets - numEntrapments}")
  
  if len(proteinGroupInfoList) == 0:
    raise Exception("No proteins with scores found, make sure that protein identifiers are consistent in the evidence and fasta files")
  
  reportedFdrs, observedFdrs, decoyLabels, proteinScores = zip(*proteinGroupInfoList)
  reportedQvals, observedQvals = fdrsToQvals(reportedFdrs), fdrsToQvals(observedFdrs)
  logger.info(f"#Targets at 1% decoy FDR: {countBelowThreshold(reportedQvals, 0.01, decoyLabels)}")
  if numEntrapments > 1:
    logger.info(f"#Targets at 1% entrapment FDR: {countBelowThreshold(observedFdrs, 0.01, decoyLabels)}")
    logger.info(f"Decoy FDR at 1% entrapment FDR: {'%.2g' % (reportedQvals[countBelowThreshold(observedFdrs, 0.01)])}")
    logger.info(f"Entrapment FDR at 1% decoy FDR: {'%.2g' % (observedFdrs[countBelowThreshold(reportedQvals, 0.01)])}")
    
    #printReportedAndEntrapmentFDRs(reportedQvals, observedQvals)
  
  return reportedQvals, observedQvals, proteinScores


def printReportedAndEntrapmentFDRs(reportedQvals, observedQvals):
  import csv
  writer = csv.writer(open('protein_fdr_calibration.txt', 'w'), delimiter = '\t')
  for reportedQval, observedQval in zip(reportedQvals, observedQvals):
    writer.writerow([reportedQval, observedQval])


def fdrsToQvals(fdrs: List[float]):
  """
  Makes a list of FDRs monotonically increasing (sometimes referred to as q-values after monotonization)
  """
  qvals = [0] * len(fdrs)
  if len(fdrs) > 0:
    qvals[len(fdrs)-1] = fdrs[-1]
    for i in range(len(fdrs)-2, -1, -1):
      qvals[i] = min(qvals[i+1], fdrs[i])
  return qvals


def countBelowThreshold(qvals: List[float], qvalThreshold: float, decoyLabels: Optional[List[bool]] = None):
  """
  Counts number of q-values below a threshold, if decoyLabels are provided, only the targets are counted
  """
  if decoyLabels is None:
    return len([1 for x in qvals if x < qvalThreshold])
  else:
    return len([1 for x, isDecoy in zip(qvals, decoyLabels) if x < qvalThreshold and not isDecoy])
 
