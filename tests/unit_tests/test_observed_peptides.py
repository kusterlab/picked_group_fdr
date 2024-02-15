import pytest
import numpy as np

import picked_group_fdr.observed_peptides as op
from picked_group_fdr.protein_groups import ProteinGroups
from picked_group_fdr.results import ProteinGroupResult


class TestCreate:
    def test_create(self, peptideInfoList):
        observed_peptides = op.ObservedPeptides()
        observed_peptides.create(peptideInfoList)

        assert observed_peptides.get_peptides("proteinA") == ["PEPTIDEA", "PEPTIDEB"]


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
