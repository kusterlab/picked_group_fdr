# Generation assisted by CodiumAI

import collections

import pytest

from picked_group_fdr.parsers.percolator import (
    parse_percolator_out_file_to_dict,
)


class TestParsePercolatorOutFileToDict:

    def test_correctly_parses_valid_percolator_output_file(self, mocker):
        mocker.patch("os.path.isfile", return_value=True)

        # Mocking the get_delimiter and get_tsv_reader functions
        mocker.patch("picked_group_fdr.parsers.tsv.get_delimiter", return_value="\t")
        mocker.patch(
            "picked_group_fdr.parsers.tsv.get_tsv_reader",
            return_value=iter(
                [
                    [
                        "PSMId",
                        "peptide",
                        "score",
                        "q-value",
                        "posterior_error_prob",
                        "proteinIds",
                    ],
                    ["file1_1_1_1", "_.PEPTIDE._", "0.9", "0.01", "0.001", "protein1"],
                ]
            ),
        )

        results_dict = collections.defaultdict(dict)
        fixed_mods_dict, results = parse_percolator_out_file_to_dict(
            "dummy_file.txt", results_dict
        )

        assert results == {"file1": {(1, "PEPTIDE"): (0.9, 0.001)}}
        assert fixed_mods_dict == {"C": "C[UNIMOD:4]"}

    def test_correctly_parses_valid_percolator_output_file_with_input_type_prosit(
        self, mocker
    ):
        mocker.patch("os.path.isfile", return_value=True)

        # Mocking the get_delimiter and get_tsv_reader functions
        mocker.patch("picked_group_fdr.parsers.tsv.get_delimiter", return_value="\t")
        mocker.patch(
            "picked_group_fdr.parsers.tsv.get_tsv_reader",
            return_value=iter(
                [
                    [
                        "PSMId",
                        "peptide",
                        "score",
                        "q-value",
                        "posterior_error_prob",
                        "proteinIds",
                    ],
                    ["file1-123-PEPTIDE-1-1", "_.PEPTIDE._", "0.9", "0.01", "0.001", "protein1"],
                ]
            ),
        )

        results_dict = collections.defaultdict(dict)
        fixed_mods_dict, results = parse_percolator_out_file_to_dict(
            "dummy_file.txt", results_dict, input_type="prosit"
        )

        assert results == {"file1": {(123, "PEPTIDE"): (0.9, 0.001)}}
        assert fixed_mods_dict == {"C": "C[UNIMOD:4]"}

    def test_correctly_identifies_and_applies_fixed_modifications_for_prosit_input_type(
        self, mocker
    ):
        mocker.patch("os.path.isfile", return_value=True)

        # Mocking the get_delimiter and get_tsv_reader functions
        mocker.patch("picked_group_fdr.parsers.tsv.get_delimiter", return_value="\t")
        mocker.patch(
            "picked_group_fdr.parsers.tsv.get_tsv_reader",
            return_value=iter(
                [
                    [
                        "PSMId",
                        "peptide",
                        "score",
                        "q-value",
                        "posterior_error_prob",
                        "proteinIds",
                    ],
                    ["file1-123-[UNIMOD:730]PEPTIDE-1-1", "_.[UNIMOD:730]PEPTIDE._", "0.9", "0.01", "0.001", "protein1"],
                ]
            ),
        )

        results_dict = collections.defaultdict(dict)
        fixed_mods_dict, results = parse_percolator_out_file_to_dict(
            "dummy_file.txt", results_dict, input_type="prosit"
        )

        assert results == {"file1": {(123, "[UNIMOD:730]PEPTIDE"): (0.9, 0.001)}}
        assert fixed_mods_dict == {'C': 'C[UNIMOD:4]', 'K': 'K[UNIMOD:730]', '^_': '_[UNIMOD:730]'}

    def test_correctly_identifies_and_applies_fixed_modifications_in_proforma_notation_for_prosit_input_type(
        self, mocker
    ):
        mocker.patch("os.path.isfile", return_value=True)

        # Mocking the get_delimiter and get_tsv_reader functions
        mocker.patch("picked_group_fdr.parsers.tsv.get_delimiter", return_value="\t")
        mocker.patch(
            "picked_group_fdr.parsers.tsv.get_tsv_reader",
            return_value=iter(
                [
                    [
                        "PSMId",
                        "peptide",
                        "score",
                        "q-value",
                        "posterior_error_prob",
                        "proteinIds",
                    ],
                    ["file1-123-[UNIMOD:730]-PEPTIDE-1-1", "_.[UNIMOD:730]-PEPTIDE._", "0.9", "0.01", "0.001", "protein1"],
                ]
            ),
        )

        results_dict = collections.defaultdict(dict)
        fixed_mods_dict, results = parse_percolator_out_file_to_dict(
            "dummy_file.txt", results_dict, input_type="prosit"
        )

        assert results == {"file1": {(123, "[UNIMOD:730]-PEPTIDE"): (0.9, 0.001)}}
        assert fixed_mods_dict == {'C': 'C[UNIMOD:4]', 'K': 'K[UNIMOD:730]', '^_': '_[UNIMOD:730]'}

    def test_raises_filenotfounderror_when_file_does_not_exist(self, mocker):
        mocker.patch("os.path.isfile", return_value=False)

        results_dict = collections.defaultdict(dict)

        with pytest.raises(FileNotFoundError):
            parse_percolator_out_file_to_dict("non_existent_file.txt", results_dict)
