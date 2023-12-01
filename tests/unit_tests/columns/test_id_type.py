from picked_group_fdr.columns.id_type import _identification_type_per_experiment
from picked_group_fdr.serializers.base import retain_only_identified_precursors


def test_getIdentificationTypePerExperiment_matchBetweenRuns(peptideIntensityListMatchBetweenRuns, experimentToIdxMap):
  postErrProbCutoff = 0.1
  peptideIntensityList = retain_only_identified_precursors(peptideIntensityListMatchBetweenRuns, postErrProbCutoff)
  assert _identification_type_per_experiment(peptideIntensityList, experimentToIdxMap, postErrProbCutoff) == ["By MS/MS", "By MS/MS", "By matching"]


def test_getIdentificationTypePerExperiment_matchBetweenRunsWithUnidentified(peptideIntensityListMatchBetweenRunsWithUnidentified, experimentToIdxMap):
  postErrProbCutoff = 0.1
  peptideIntensityList = retain_only_identified_precursors(peptideIntensityListMatchBetweenRunsWithUnidentified, postErrProbCutoff)
  assert _identification_type_per_experiment(peptideIntensityList, experimentToIdxMap, postErrProbCutoff) == ["", "", ""]


