import pandas as pd
import pytest
from io import StringIO
from unittest.mock import patch

import picked_group_fdr.parsers.modifications as mod
import picked_group_fdr.parsers.parsers as parsers


class TestMaxQuantToProForma:
    def test_maxquant_to_internal_carbamidomethylation(self):
        assert mod.maxquant_mod_to_proforma("_ABCDEFGH_") == "ABC[UNIMOD:4]DEFGH"

    def test_maxquant_to_internal_tmt(self):
        fixed_mods = {"C": "C[UNIMOD:4]", "^_": "_[UNIMOD:737]-", "K": "K[UNIMOD:737]"}
        assert (
            mod.maxquant_mod_to_proforma("_ABCDEFGHK_", fixed_mods)
            == "[UNIMOD:737]-ABC[UNIMOD:4]DEFGHK[UNIMOD:737]"
        )

    def test_maxquant_to_internal_silac(self):
        fixed_mods = {"C": "C[UNIMOD:4]", "K": "K[UNIMOD:259]", "R": "R[UNIMOD:267]"}
        assert (
            mod.maxquant_mod_to_proforma("_ABCDEFGHRK_", fixed_mods)
            == "ABC[UNIMOD:4]DEFGHR[UNIMOD:267]K[UNIMOD:259]"
        )


class TestPrositToProForma:
    def test_prosit_to_internal_carbamidomethylation(self):
        assert mod.prosit_mod_to_proforma("ABC[UNIMOD:4]DEFGH") == "ABC[UNIMOD:4]DEFGH"

    def test_prosit_to_internal_tmt(self):
        assert (
            mod.prosit_mod_to_proforma("[UNIMOD:737]ABC[UNIMOD:4]DEFGHK[UNIMOD:737]")
            == "[UNIMOD:737]-ABC[UNIMOD:4]DEFGHK[UNIMOD:737]"
        )

    def test_prosit_c_terminal_mod(self):
        assert (
            mod.prosit_mod_to_proforma(
                "[UNIMOD:737]ABC[UNIMOD:4]DEFGHK[UNIMOD:737]",
                variable_mods={"[UNIMOD:737]$": "-[UNIMOD:737]"},
            )
            == "[UNIMOD:737]-ABC[UNIMOD:4]DEFGHK-[UNIMOD:737]"
        )


@pytest.fixture
def triqler_file_list():
    csv_string = """file1.raw\tgroup1
file2.RAW\tgroup2
file3.mzML\tgroup1
"""
    return pd.read_csv(
        StringIO(csv_string),
        sep="\t",
        header=None,
        names=["Name", "Condition", "Experiment", "Fraction"],
    )


def test_parse_triqler_file_list(triqler_file_list):
    expected_output = pd.DataFrame(
        {
            "Name": ["file1", "file2", "file3"],
            "Condition": ["group1", "group2", "group1"],
            "Experiment": ["file1", "file2", "file3"],
            "Fraction": [-1, -1, -1],
        }
    )

    with patch("pandas.read_csv", return_value=triqler_file_list):
        result = parsers.parse_triqler_file_list("fake_input_file.tsv")

    pd.testing.assert_frame_equal(result, expected_output, check_dtype=False)


@pytest.fixture
def mq_experimental_design():
    csv_string = """Name\tExperiment\tFraction
file1.raw\texperiment1\t1
file2.RAW\texperiment2\t1
file3.mzML\texperiment1\t2
"""
    return pd.read_csv(StringIO(csv_string), sep="\t")


def test_parse_mq_experimental_design(mq_experimental_design):
    expected_output = pd.DataFrame(
        {
            "Name": ["file1", "file2", "file3"],
            "Condition": ["experiment1", "experiment2", "experiment1"],
            "Experiment": ["experiment1", "experiment2", "experiment1"],
            "Fraction": [1, 1, 2],
        }
    )

    with patch("pandas.read_csv", return_value=mq_experimental_design):
        result = parsers.parse_mq_experimental_design("fake_input_file.tsv")

    pd.testing.assert_frame_equal(result, expected_output, check_dtype=False)


def test_add_triqler_group_params():
    params = {"groupLabels": [], "groups": []}

    expected_params = {
        "groupLabels": ["1:condition1", "2:condition2"],
        "groups": [[0, 2, 0], [1]],
    }

    columns = ["Name", "Condition", "Experiment", "Fraction"]
    data = [
        ("file1", "condition1", "experiment1", 1),
        ("file2", "condition2", "experiment2", 2),
        ("file3", "condition1", "experiment3", 3),
        ("file4", "condition1", "experiment1", 2),
    ]
    df = pd.DataFrame(data, columns=columns)

    result = parsers.add_triqler_group_params(df, params)

    assert result == expected_params


@pytest.fixture
def sample_file_info_dataframe():
    return pd.DataFrame(
        {
            "Name": ["file1", "file2", "file3"],
            "Experiment": ["experiment1", "experiment2", "experiment3"],
            "Fraction": [1, 2, 3],
        }
    )


# Test function using the sample_file_info_dataframe fixture
def test_get_file_mapping(sample_file_info_dataframe):
    expected_file_mapping = {
        "file1": ("experiment1", 1),
        "file2": ("experiment2", 2),
        "file3": ("experiment3", 3),
    }

    result = parsers.get_file_mapping(sample_file_info_dataframe)

    assert result == expected_file_mapping
