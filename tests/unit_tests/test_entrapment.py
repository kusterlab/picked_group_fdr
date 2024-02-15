import pytest

import picked_group_fdr.entrapment as entrapment


def test_mark_entrapment(proteins, entrapment_proteins):
    assert entrapment.mark_entrapment(proteins, entrapment_proteins) == [
        "proteinA_entrapment",
        "proteinB",
        "proteinC_entrapment",
    ]


def test_is_entrapment_false():
    assert (
        entrapment.is_entrapment(
            ["proteinA_entrapment", "proteinB", "proteinC_entrapment"]
        )
        == False
    )


def test_is_entrapment_true():
    assert (
        entrapment.is_entrapment(
            ["proteinA_entrapment", "REV__proteinB", "proteinC_entrapment"]
        )
        == True
    )


@pytest.fixture
def proteins():
    return ["proteinA", "proteinB", "proteinC"]


@pytest.fixture
def entrapment_proteins():
    return set(["proteinA_entrapment", "proteinC_entrapment"])
