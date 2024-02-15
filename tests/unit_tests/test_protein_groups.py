import pytest

import picked_group_fdr.protein_groups as protein_groups


class TestProteinGrouping:
    def test_proteinToGroupMap(self, proteinGroups):
        assert proteinGroups.protein_to_group_idx_map == {
            "proteinA": 0,
            "proteinE": 0,
            "proteinB": 1,
            "proteinC": 2,
            "proteinD": 3,
            "proteinF": 3,
        }

    def test_getLeadingProteins(self, proteinGroups):
        proteins = ["proteinA", "proteinE", "proteinD"]
        assert proteinGroups.get_leading_proteins(proteins) == {"proteinA", "proteinD"}


@pytest.fixture
def proteinGroups():
    proteinGroups = protein_groups.ProteinGroups(
        [["proteinA", "proteinE"], ["proteinB"], ["proteinC"], ["proteinD", "proteinF"]]
    )
    proteinGroups.create_index()
    return proteinGroups
