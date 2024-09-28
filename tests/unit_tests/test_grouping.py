import pytest
import numpy as np

import picked_group_fdr.grouping as grouping
from picked_group_fdr.protein_groups import ProteinGroups
from picked_group_fdr.results import ProteinGroupResult
from picked_group_fdr.scoring_strategy import ProteinScoringStrategy


class TestNoGrouping:
    def test_group_proteins(self, peptideInfoList):
        no_grouping = grouping.NoGrouping()
        assert no_grouping.group_proteins(peptideInfoList, "") == ProteinGroups(
            [["proteinA"], ["proteinB"], ["proteinC"]]
        )


class TestSubsetGrouping:
    def test_group_proteins(self, peptideInfoList):
        subset_grouping = grouping.SubsetGrouping()
        assert subset_grouping.group_proteins(peptideInfoList, "") == ProteinGroups(
            [["proteinA", "proteinB"], ["proteinC"]]
        )

    def test_group_proteins_each_unique(self, peptideInfoListRescue):
        subset_grouping = grouping.SubsetGrouping()
        assert subset_grouping.group_proteins(
            peptideInfoListRescue, ""
        ) == ProteinGroups([["proteinA"], ["proteinB"], ["proteinC"]])


class TestRescuedSubsetGrouping:
    def test_group_proteins(self, peptideInfoListRescue, proteinFdrResults):
        rescued_subset_grouping = grouping.RescuedSubsetGrouping()
        protein_groups_before_rescue = rescued_subset_grouping.group_proteins(
            peptideInfoListRescue, ""
        )
        score_type = ProteinScoringStrategy("bestPEP")
        protein_group_peptide_infos = score_type.collect_peptide_scores_per_protein(
            protein_groups_before_rescue,
            peptideInfoListRescue,
            peptide_qval_cutoff=0.01,
        )

        assert protein_groups_before_rescue == ProteinGroups(
            [["proteinA"], ["proteinB"], ["proteinC"]]
        )

        protein_group_fdr = 0.01
        protein_groups_after_rescue = rescued_subset_grouping.rescue_protein_groups(
            peptideInfoListRescue,
            proteinFdrResults,
            protein_group_fdr,
            protein_groups_before_rescue,
            protein_group_peptide_infos,
        )
        assert protein_groups_after_rescue == ProteinGroups(
            [["proteinA", "proteinB"], ["proteinC"]]
        )

    def test_scoring_threshold(self, proteinFdrResults):
        rescued_subset_grouping = grouping.RescuedSubsetGrouping()
        protein_group_fdr = 0.01
        rescued_subset_grouping._calculate_rescue_score_cutoff(
            proteinFdrResults, protein_group_fdr
        )
        assert rescued_subset_grouping.score_cutoff == 0.15


@pytest.fixture
def peptideInfoList():
    """
    dictionary of peptide -> (score, proteins)
    """
    return {
        "PEPTIDEA": (0.12, ["proteinA", "proteinB"]),
        "PEPTIDEB": (0.12, ["proteinA", "proteinC"]),
        "PEPTIDEC": (0.12, ["proteinC"]),
    }


@pytest.fixture
def peptideInfoListRescue():
    """
    dictionary of peptide -> (score, proteins)
    """
    return {
        "PEPTIDEA": (0.12, ["proteinA", "proteinB"]),
        "PEPTIDEB": (0.12, ["proteinA", "proteinC"]),
        "PEPTIDEC": (0.12, ["proteinC"]),
        "PEPTIDED": (0.35, ["proteinA"]),
        "PEPTIDED": (0.35, ["proteinB"]),
    }  # the last two peptides should be removed by the rescuing procedure


@pytest.fixture
def proteinFdrResults():
    """
    list of ProteinGroupResult
    """
    pgr = ProteinGroupResult()
    pgr.score = -1 * np.log10(0.15)
    pgr.qValue = 0.001

    pgr2 = ProteinGroupResult()
    pgr2.score = -1 * np.log10(0.16)
    pgr2.qValue = 0.02
    return [pgr, pgr2]
