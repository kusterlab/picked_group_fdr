from typing import List, Dict
import itertools
import collections
import warnings
import logging

import numpy as np
import triqler
import triqler.parsers
import triqler.hyperparameters
import triqler.qvality
import triqler.triqler

from .. import helpers
from .base import ProteinGroupColumns


logger = logging.getLogger(__name__)


class TriqlerIntensityColumns(ProteinGroupColumns):  
  def __init__(self, params):
    self.params = params
  
  def append_headers(self, proteinGroupResults, experiments):
    numGroups = len(self.params["groups"])
    if numGroups <= 1:
      logger.info("Skipping quantification with Triqler, less than 2 conditions found in --file_list_file input")
      return
    
    for experiment in experiments:
      proteinGroupResults.append_header('Triqler Relative Quant ' + experiment)
    
    for groupId1, groupId2 in itertools.combinations(self.params["groupLabels"], 2):
      proteinGroupResults.append_header('Triqler log2_fc ' + str(groupId1) + ' vs ' + str(groupId2))
      proteinGroupResults.append_header('Triqler PEP_diff_exp_' + str(self.params['foldChangeEval']) + ' ' + str(groupId1) + ' vs ' + str(groupId2))
  
  def append_columns(self, proteinGroupResults, experimentToIdxMap, postErrProbCutoff):
    numGroups = len(self.params["groups"])
    if numGroups <= 1:
      return
    
    logger.info("Doing quantification: Triqler relative intensity")
    
    pickedProteinOutputRows = list()
    proteinGroupResultsTriqler = list()
    decoyScores, targetScores = list(), list()
    for pgr in proteinGroupResults:
      leadingProtein = pgr.proteinIds.split(";")[0]
      quantRows = self._addPeptideQuantRows(pgr.precursorQuants, experimentToIdxMap, 
                                            leadingProtein, self.params["minSamples"])
      if len(quantRows) > 0:
        if pgr.reverse == '+':
          decoyScores.append(pgr.score)
        else:
          targetScores.append(pgr.score)
        
        # pickedProteinOutputRow = (combinedPEP, protein, quantRows, numPeptides)
        # only quantRows is actually processed by Triqler
        # TODO: clean up Triqler so we only need to input a list of quantRows
        pickedProteinOutputRows.append([0.0, leadingProtein, quantRows, len(quantRows)])
        proteinGroupResultsTriqler.append(pgr)
      else:
        numComparisons = int(numGroups*(numGroups-1)/2)
        pgr.extend([0.0]*len(experimentToIdxMap) + [0.0, 1.0]*numComparisons)
    
    proteinIdPEPs = self._calculateProteinIdPEPs(targetScores, decoyScores)
    
    posteriors = self._getTriqlerPosteriors(pickedProteinOutputRows, proteinIdPEPs)
    
    for pgr, triqlerResults, proteinIdPEP in zip(proteinGroupResultsTriqler, posteriors, proteinIdPEPs):
      bayesQuantRow, muGroupDiffs, probsBelowFoldChange, posteriorDists = triqlerResults
      
      # format as string, as the default formatting is as an integer
      pgr.extend(["%.2f" % x for x in bayesQuantRow])
      for groupId1, groupId2 in itertools.combinations(range(numGroups), 2):
        pgr.append("%.2f" % muGroupDiffs[(groupId1, groupId2)])
        
        # combine the error probabilities of the protein is either not 
        # differentially expressed or not correctly identified
        foldChangePEP = probsBelowFoldChange[(groupId1, groupId2)]
        differentialExpressionPEP = 1.0 - (1.0 - foldChangePEP) * (1.0 - proteinIdPEP)
        pgr.append("%.2f" % differentialExpressionPEP)
    # TODO: print Triqler summary results to the command prompt
  
  def _getTriqlerPosteriors(self, pickedProteinOutputRows, proteinIdPEPs):
    peptideQuantRows = list()
    for row in pickedProteinOutputRows:
      _, _, quantRows, _ = row
      peptideQuantRows.extend(quantRows)
    
    with warnings.catch_warnings():
      warnings.simplefilter(self.params['warningFilter'])
      triqler.hyperparameters.fitPriors(peptideQuantRows, self.params)
      posteriors = triqler.triqler.getPosteriors(pickedProteinOutputRows, proteinIdPEPs, self.params)
    return posteriors

  def _calculateProteinIdPEPs(self, targetScores, decoyScores):  
    targetScores = np.array(targetScores)
    decoyScores = np.array(decoyScores)
    _, proteinIdPEPs = triqler.qvality.getQvaluesFromScores(
        targetScores, decoyScores, 
        includePEPs = True, includeDecoys = True, tdcInput = True)  
    return proteinIdPEPs
  
  def _addPeptideQuantRows(self, peptideIntensityList, experimentToIdxMap, 
                           leadingProtein, minSamples):
    peptideIntensities = collections.defaultdict(lambda : [[0.0, 1.01] for _ in range(len(experimentToIdxMap))])
    prevExpFrac = (None, None)
    prevPrecursor = (None, None)
    # for each (peptide, charge, experiment, fraction) tuple, 
    # sort the lowest PEP on top, break ties by higher intensity
    orderByPEPAndIntensity = lambda x : (x.peptide, x.charge, x.experiment, x.fraction, x.postErrProb if not np.isnan(x.postErrProb) else 1.01, -1*x.intensity)
    for precursor in sorted(peptideIntensityList, key = orderByPEPAndIntensity):
      if precursor.intensity > 0.0:
        expIdx = experimentToIdxMap[precursor.experiment]
        if prevExpFrac != (precursor.experiment, precursor.fraction) or prevPrecursor != (precursor.peptide, precursor.charge):
          peptideIntensities[(precursor.peptide, precursor.charge)][expIdx][0] += precursor.intensity
          if not helpers.isMbr(precursor.postErrProb):
            currentPostErrProb = peptideIntensities[(precursor.peptide, precursor.charge)][expIdx][1]
            peptideIntensities[(precursor.peptide, precursor.charge)][expIdx][1] = min([currentPostErrProb, precursor.postErrProb])
        prevExpFrac = (precursor.experiment, precursor.fraction)
        prevPrecursor = (precursor.peptide, precursor.charge)
    
    quantRows = list()
    for (peptide, charge), quantPEPpairs in peptideIntensities.items():
      quants, identificationPEPs = zip(*quantPEPpairs)
      maxPEP = max([x for x in identificationPEPs if x <= 1.0])
      identificationPEPs = [x if x <= 1.0 else maxPEP for x in identificationPEPs]
      if np.count_nonzero(quants) >= minSamples:
        quantRows.append(triqler.parsers.PeptideQuantRow(
            min(identificationPEPs), 
            charge, 
            0, # feature group idx
            0, # spectrum idx
            [0.0]*len(experimentToIdxMap), # match-between-runs PEPs
            quants, 
            identificationPEPs, 
            peptide, 
            [leadingProtein]))
    
    return quantRows
