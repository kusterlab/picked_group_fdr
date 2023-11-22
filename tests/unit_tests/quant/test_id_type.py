from picked_group_fdr.quant.id_type import _identification_type_per_experiment
from picked_group_fdr.quantification import retainOnlyIdentifiedPrecursors


def test_getIdentificationTypePerExperiment_matchBetweenRuns(peptideIntensityListMatchBetweenRuns, experimentToIdxMap):
  postErrProbCutoff = 0.1
  peptideIntensityList = retainOnlyIdentifiedPrecursors(peptideIntensityListMatchBetweenRuns, postErrProbCutoff)
  assert _identification_type_per_experiment(peptideIntensityList, experimentToIdxMap, postErrProbCutoff) == ["By MS/MS", "By MS/MS", "By matching"]


def test_getIdentificationTypePerExperiment_matchBetweenRunsWithUnidentified(peptideIntensityListMatchBetweenRunsWithUnidentified, experimentToIdxMap):
  postErrProbCutoff = 0.1
  peptideIntensityList = retainOnlyIdentifiedPrecursors(peptideIntensityListMatchBetweenRunsWithUnidentified, postErrProbCutoff)
  assert _identification_type_per_experiment(peptideIntensityList, experimentToIdxMap, postErrProbCutoff) == ["", "", ""]


