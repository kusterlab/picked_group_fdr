import pytest
import numpy as np

import picked_group_fdr.protein_groups as protein_groups
from picked_group_fdr.scoring import MultPEPScore


class TestMultPEPScore:
    def test_optimize_hyperparameters(
        self, proteinGroupsSimple, proteinGroupPeptideInfosSimple
    ):
        score_type = MultPEPScore()
        protein_group_fdr_threshold = 0.01
        score_type.optimize_hyperparameters(
            proteinGroupsSimple,
            proteinGroupPeptideInfosSimple,
            protein_group_fdr_threshold,
        )

        np.testing.assert_almost_equal(score_type.div, 1e-4)


@pytest.fixture
def proteinGroupsSimple():
    proteinGroups = protein_groups.ProteinGroups(
        [
            ["proteinA", "proteinE"],
            ["REV__proteinA", "REV__proteinC"],
            ["REV__proteinB", "REV__proteinD"],
        ]
    )
    proteinGroups.create_index()
    return proteinGroups


@pytest.fixture
def proteinGroupPeptideInfosSimple():
    # tuples of (score, peptide, proteins)
    return [
        [(0.001, "A", []), (0.9, "B", []), (0.1, "A", [])],
        [(0.01, "", [])],
        [(0.009, "", [])],
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
