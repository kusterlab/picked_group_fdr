# Generation assisted by CodiumAI

import collections

import pytest
from contextlib import nullcontext

from picked_group_fdr.parsers.percolator import (
    parse_percolator_out_file_to_dict,
    parse_prosit_psmid_and_peptide,
    parse_andromeda_psmid_and_peptide,
)
import picked_group_fdr.parsers.modifications as modifications


class TestParsePercolatorOutFileToDict:

    def test_correctly_parses_valid_percolator_output_file(self, mocker):
        mocker.patch("os.path.isfile", return_value=True)

        # Mocking the get_delimiter and get_tsv_reader functions
        mocker.patch("picked_group_fdr.parsers.tsv.get_delimiter", return_value="\t")
        mocker.patch(
            "picked_group_fdr.parsers.tsv.get_tsv_reader",
            return_value=nullcontext(
                iter(
                    [
                        [
                            "PSMId",
                            "peptide",
                            "score",
                            "q-value",
                            "posterior_error_prob",
                            "proteinIds",
                        ],
                        [
                            "file1_1_1_1",
                            "_.PEPTIDE._",
                            "0.9",
                            "0.01",
                            "0.001",
                            "protein1",
                        ],
                    ]
                )
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
            return_value=nullcontext(
                iter(
                    [
                        [
                            "PSMId",
                            "peptide",
                            "score",
                            "q-value",
                            "posterior_error_prob",
                            "proteinIds",
                        ],
                        [
                            "file1-123-PEPTIDE-1-1",
                            "_.PEPTIDE._",
                            "0.9",
                            "0.01",
                            "0.001",
                            "protein1",
                        ],
                    ]
                )
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
            return_value=nullcontext(
                iter(
                    [
                        [
                            "PSMId",
                            "peptide",
                            "score",
                            "q-value",
                            "posterior_error_prob",
                            "proteinIds",
                        ],
                        [
                            "file1-123-[UNIMOD:730]PEPTIDE-1-1",
                            "_.[UNIMOD:730]PEPTIDE._",
                            "0.9",
                            "0.01",
                            "0.001",
                            "protein1",
                        ],
                    ]
                )
            ),
        )

        results_dict = collections.defaultdict(dict)
        fixed_mods_dict, results = parse_percolator_out_file_to_dict(
            "dummy_file.txt", results_dict, input_type="prosit"
        )

        assert results == {"file1": {(123, "[UNIMOD:730]-PEPTIDE"): (0.9, 0.001)}}
        assert fixed_mods_dict == {
            "C": "C[UNIMOD:4]",
            "K": "K[UNIMOD:730]",
            "^_": "_[UNIMOD:730]-",
        }

    def test_correctly_identifies_and_applies_fixed_modifications_in_proforma_notation_for_prosit_input_type(
        self, mocker
    ):
        mocker.patch("os.path.isfile", return_value=True)

        # Mocking the get_delimiter and get_tsv_reader functions
        mocker.patch("picked_group_fdr.parsers.tsv.get_delimiter", return_value="\t")
        mocker.patch(
            "picked_group_fdr.parsers.tsv.get_tsv_reader",
            return_value=nullcontext(
                iter(
                    [
                        [
                            "PSMId",
                            "peptide",
                            "score",
                            "q-value",
                            "posterior_error_prob",
                            "proteinIds",
                        ],
                        [
                            "file1-123-[UNIMOD:730]-PEPTIDE-1-1",
                            "_.[UNIMOD:730]-PEPTIDE._",
                            "0.9",
                            "0.01",
                            "0.001",
                            "protein1",
                        ],
                    ]
                )
            ),
        )

        results_dict = collections.defaultdict(dict)
        fixed_mods_dict, results = parse_percolator_out_file_to_dict(
            "dummy_file.txt", results_dict, input_type="prosit"
        )

        assert results == {"file1": {(123, "[UNIMOD:730]-PEPTIDE"): (0.9, 0.001)}}
        assert fixed_mods_dict == {
            "C": "C[UNIMOD:4]",
            "K": "K[UNIMOD:730]",
            "^_": "_[UNIMOD:730]-",
        }

    def test_raises_filenotfounderror_when_file_does_not_exist(self, mocker):
        mocker.patch("os.path.isfile", return_value=False)

        results_dict = collections.defaultdict(dict)

        with pytest.raises(FileNotFoundError):
            parse_percolator_out_file_to_dict("non_existent_file.txt", results_dict)


class TestParsePrositPsmidAndPeptide:
    def test_with_scan_event(self):
        psm_id = "file1_with-multiple-dashes-123-[UNIMOD:730]-PEPTIDE-1-1"
        peptide = "[UNIMOD:730]-PEPTIDE"
        filename = "file1_with-multiple-dashes"
        convert_to_proforma = modifications.prosit_mod_to_proforma()
        raw_file, scan_number, converted_peptide = parse_prosit_psmid_and_peptide(
            psm_id, peptide, filename, convert_to_proforma
        )
        assert raw_file == "file1_with-multiple-dashes"
        assert scan_number == 123
        assert converted_peptide == "[UNIMOD:730]-PEPTIDE"

    def test_with_scan_event_without_filename(self):
        psm_id = "file1_with-multiple-dashes-123-[UNIMOD:730]-PEPTIDE-1-1"
        peptide = "[UNIMOD:730]-PEPTIDE"
        filename = ""
        convert_to_proforma = modifications.prosit_mod_to_proforma()
        raw_file, scan_number, converted_peptide = parse_prosit_psmid_and_peptide(
            psm_id, peptide, filename, convert_to_proforma
        )
        assert raw_file == "file1_with-multiple-dashes"
        assert scan_number == 123
        assert converted_peptide == "[UNIMOD:730]-PEPTIDE"

    def test_without_scan_event(self):
        psm_id = "file1_with-multiple-dashes-123-[UNIMOD:730]-PEPTIDE-1"
        peptide = "[UNIMOD:730]-PEPTIDE"
        filename = "file1_with-multiple-dashes"
        convert_to_proforma = modifications.prosit_mod_to_proforma()
        raw_file, scan_number, converted_peptide = parse_prosit_psmid_and_peptide(
            psm_id, peptide, filename, convert_to_proforma
        )
        assert raw_file == "file1_with-multiple-dashes"
        assert scan_number == 123
        assert converted_peptide == "[UNIMOD:730]-PEPTIDE"


class TestParseAndromedaPsmidAndPeptide:
    def test_valid_psm_id_and_modified_sequence(self):
        psm_id = "file1_1234_5678_9012"
        modified_sequence = "ACDEFGHIK[42]LM[16]NOP"
        expected_output = ("file1", 1234, "ACDEFGHIK(ac)LM(ox)NOP")
        result = parse_andromeda_psmid_and_peptide(psm_id, modified_sequence)
        assert result == expected_output
