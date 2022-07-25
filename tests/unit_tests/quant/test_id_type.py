import pytest
import numpy as np

from picked_group_fdr.quant.id_type import IdentificationTypeColumns
from picked_group_fdr.quantification import retainOnlyIdentifiedPrecursors


def test_getIdentificationTypePerExperiment_matchBetweenRuns(peptideIntensityListMatchBetweenRuns, experimentToIdxMap):
  pq = IdentificationTypeColumns()
  postErrProbCutoff = 0.1
  peptideIntensityList = retainOnlyIdentifiedPrecursors(peptideIntensityListMatchBetweenRuns, postErrProbCutoff)
  assert pq._identificationTypePerExperiment(peptideIntensityList, experimentToIdxMap, postErrProbCutoff) == ["By MS/MS", "By MS/MS", "By matching"]


def test_getIdentificationTypePerExperiment_matchBetweenRunsWithUnidentified(peptideIntensityListMatchBetweenRunsWithUnidentified, experimentToIdxMap):
  pq = IdentificationTypeColumns()
  postErrProbCutoff = 0.1
  peptideIntensityList = retainOnlyIdentifiedPrecursors(peptideIntensityListMatchBetweenRunsWithUnidentified, postErrProbCutoff)
  assert pq._identificationTypePerExperiment(peptideIntensityList, experimentToIdxMap, postErrProbCutoff) == ["", "", ""]


