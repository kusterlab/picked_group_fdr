import pytest

import picked_group_fdr.competition as competition
import picked_group_fdr.protein_groups as protein_groups
from picked_group_fdr.scoring import BestPEPScore


class TestPickedStrategy:
    def test_do_competition(self, proteinGroupsSimple, proteinGroupPeptideInfosSimple):
        score_type = BestPEPScore()
        picked_strategy = competition.PickedStrategy()
        filtered_protein_groups, filtered_protein_group_scores, _ = (
            picked_strategy.do_competition(
                proteinGroupsSimple, proteinGroupPeptideInfosSimple, score_type
            )
        )

        assert filtered_protein_groups.get_all_proteins() == set(
            ["proteinA", "proteinE", "REV__proteinA", "REV__proteinC"]
        )
        assert filtered_protein_group_scores == proteinGroupPeptideInfosSimple

    def test_do_competition_no_grouping(
        self, proteinGroupsNoGrouping, proteinGroupPeptideInfosNoGrouping
    ):
        score_type = BestPEPScore()
        picked_strategy = competition.PickedStrategy()
        filtered_protein_groups, filtered_protein_group_scores, _ = (
            picked_strategy.do_competition(
                proteinGroupsNoGrouping, proteinGroupPeptideInfosNoGrouping, score_type
            )
        )

        assert filtered_protein_groups.get_all_proteins() == set(
            ["proteinA", "proteinE", "REV__proteinC"]
        )
        assert filtered_protein_group_scores == [
            [(0.001, "", [])],
            [(0.001, "", [])],
            [(0.01, "", [])],
        ]

    def test_do_competition_unordered(
        self, proteinGroupsUnordered, proteinGroupPeptideInfosUnordered
    ):
        score_type = BestPEPScore()
        picked_strategy = competition.PickedStrategy()
        filtered_protein_groups, filtered_protein_group_scores, _ = (
            picked_strategy.do_competition(
                proteinGroupsUnordered, proteinGroupPeptideInfosUnordered, score_type
            )
        )

        assert filtered_protein_groups.get_all_proteins() == set(
            ["proteinA", "proteinE", "REV__proteinA", "REV__proteinC"]
        )
        assert filtered_protein_group_scores == [[(0.001, "", [])], [(0.01, "", [])]]


class TestPickedGroupStrategy:
    def test_do_competition(self, proteinGroupsSimple, proteinGroupPeptideInfosSimple):
        score_type = BestPEPScore()
        picked_strategy = competition.PickedGroupStrategy()
        filtered_protein_groups, filtered_protein_group_scores, _ = (
            picked_strategy.do_competition(
                proteinGroupsSimple, proteinGroupPeptideInfosSimple, score_type
            )
        )

        assert filtered_protein_groups.get_all_proteins() == set(
            ["proteinA", "proteinE"]
        )
        assert filtered_protein_group_scores == [
            [(0.001, "pepA", ["proteinA", "proteinE"])]
        ]

    def test_do_competition_no_grouping(
        self, proteinGroupsNoGrouping, proteinGroupPeptideInfosNoGrouping
    ):
        score_type = BestPEPScore()
        picked_strategy = competition.PickedGroupStrategy()
        filtered_protein_groups, filtered_protein_group_scores, _ = (
            picked_strategy.do_competition(
                proteinGroupsNoGrouping, proteinGroupPeptideInfosNoGrouping, score_type
            )
        )

        assert filtered_protein_groups.get_all_proteins() == set(
            ["proteinA", "proteinE", "REV__proteinC"]
        )
        assert filtered_protein_group_scores == [
            [(0.001, "", [])],
            [(0.001, "", [])],
            [(0.01, "", [])],
        ]

    def test_do_competition_unordered(
        self, proteinGroupsUnordered, proteinGroupPeptideInfosUnordered
    ):
        score_type = BestPEPScore()
        picked_strategy = competition.PickedGroupStrategy()
        filtered_protein_groups, filtered_protein_group_scores, _ = (
            picked_strategy.do_competition(
                proteinGroupsUnordered, proteinGroupPeptideInfosUnordered, score_type
            )
        )

        assert filtered_protein_groups.get_all_proteins() == set(
            ["REV__proteinA", "REV__proteinC"]
        )
        assert filtered_protein_group_scores == [[(0.001, "", [])]]


class TestClassicStrategy:
    def test_do_competition(self, proteinGroupsSimple, proteinGroupPeptideInfosSimple):
        score_type = BestPEPScore()
        classic_strategy = competition.ClassicStrategy()
        filtered_protein_groups, filtered_protein_group_scores, _ = (
            classic_strategy.do_competition(
                proteinGroupsSimple, proteinGroupPeptideInfosSimple, score_type
            )
        )

        assert filtered_protein_groups.get_all_proteins() == set(
            ["proteinA", "proteinE", "REV__proteinA", "REV__proteinC"]
        )
        assert filtered_protein_group_scores == proteinGroupPeptideInfosSimple

    def test_do_competition_no_grouping(
        self, proteinGroupsNoGrouping, proteinGroupPeptideInfosNoGrouping
    ):
        score_type = BestPEPScore()
        picked_strategy = competition.ClassicStrategy()
        filtered_protein_groups, filtered_protein_group_scores, _ = (
            picked_strategy.do_competition(
                proteinGroupsNoGrouping, proteinGroupPeptideInfosNoGrouping, score_type
            )
        )

        assert filtered_protein_groups.get_all_proteins() == set(
            ["proteinA", "proteinE", "REV__proteinA", "REV__proteinC"]
        )
        assert filtered_protein_group_scores == [
            [(0.001, "", [])],
            [(0.001, "", [])],
            [(0.01, "", [])],
            [(0.01, "", [])],
        ]

    def test_do_competition_unordered(
        self, proteinGroupsUnordered, proteinGroupPeptideInfosUnordered
    ):
        score_type = BestPEPScore()
        classic_strategy = competition.ClassicStrategy()
        filtered_protein_groups, filtered_protein_group_scores, _ = (
            classic_strategy.do_competition(
                proteinGroupsUnordered, proteinGroupPeptideInfosUnordered, score_type
            )
        )

        assert filtered_protein_groups.get_all_proteins() == set(
            ["proteinA", "proteinE", "REV__proteinA", "REV__proteinC"]
        )
        assert filtered_protein_group_scores == [[(0.001, "", [])], [(0.01, "", [])]]


@pytest.fixture
def proteinGroupsSimple():
    proteinGroups = protein_groups.ProteinGroups(
        [["proteinA", "proteinE"], ["REV__proteinA", "REV__proteinC"]]
    )
    proteinGroups.create_index()
    return proteinGroups


@pytest.fixture
def proteinGroupPeptideInfosSimple():
    # tuples of (score, peptide, proteins)
    return [
        [(0.001, "pepA", ["proteinA", "proteinE"])],
        [(0.01, "pepB", ["REV__proteinA", "REV__proteinC"])],
    ]


@pytest.fixture
def proteinGroupsNoGrouping():
    proteinGroups = protein_groups.ProteinGroups(
        [["proteinA"], ["proteinE"], ["REV__proteinA"], ["REV__proteinC"]]
    )
    proteinGroups.create_index()
    return proteinGroups


@pytest.fixture
def proteinGroupPeptideInfosNoGrouping():
    # tuples of (score, peptide, proteins)
    return [[(0.001, "", [])], [(0.001, "", [])], [(0.01, "", [])], [(0.01, "", [])]]


@pytest.fixture
def proteinGroupsUnordered():
    proteinGroups = protein_groups.ProteinGroups(
        [["proteinA", "proteinE"], ["REV__proteinA", "REV__proteinC"]]
    )
    proteinGroups.create_index()
    return proteinGroups


@pytest.fixture
def proteinGroupPeptideInfosUnordered():
    # tuples of (score, peptide, proteins)
    return [[(0.01, "", [])], [(0.001, "", [])]]
