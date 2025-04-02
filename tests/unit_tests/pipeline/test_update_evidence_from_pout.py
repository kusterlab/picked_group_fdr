# Generation assisted by CodiumAI

import pytest
from contextlib import nullcontext

import picked_group_fdr.parsers.maxquant as maxquant
from picked_group_fdr.pipeline.update_evidence_from_pout import update_evidence_single


class TestUpdateEvidenceSingle:

    def test_correctly_processes_evidence_file_and_writes_updated_rows(self, mocker):
        # Mock the necessary components
        mock_writer = mocker.patch("csv.writer")
        mocker.patch(
            "picked_group_fdr.parsers.tsv.get_tsv_reader",
            return_value=nullcontext(
                iter(
                    [
                        ["Header1", "Header2", "Header3"],
                        ["row1col1", "row1col2", "row1col3"],
                    ]
                )
            ),
        )
        # the first index is for the score column, the second for PEP
        mocker.patch(
            "picked_group_fdr.parsers.tsv.get_column_index", side_effect=[0, 1]
        )
        mocker.patch(
            "picked_group_fdr.parsers.maxquant.parse_evidence_file_for_percolator_matching",
            return_value=iter(
                [
                    (
                        ["row1col1", "row1col2", "row1col3"],
                        maxquant.EvidenceRow(
                            raw_file="file1.raw",
                            scannr=123,
                            peptide="_PEPTIDESEQ_",
                            score=95.5,
                            post_err_prob=0.05,
                            is_decoy=False,
                            is_contaminant=False,
                            id_type="MSMS",
                            labeling_state="NOT_SILAC",
                        ),
                    )
                ]
            ),
        )

        # Initialize the necessary variables
        evidence_file = "path/to/evidence_file.txt"
        first_headers = ["Header1", "Header2", "Header3"]
        fixed_mods = {"C": "C[UNIMOD:4]"}
        results_dict = {
            "file1.raw": {
                (123, "PEPTIDESEQ"): (0.1, 0.01),
                (2, "peptide2"): (0.2, 0.02),
            }
        }
        pout_input_type = "andromeda"
        suppress_missing_peptide_warning = False

        # Call the function
        updated_headers = update_evidence_single(
            evidence_file,
            mock_writer,
            first_headers,
            fixed_mods,
            results_dict,
            maxquant.parse_evidence_file_for_percolator_matching,
            pout_input_type,
            suppress_missing_peptide_warning,
        )

        # Assertions
        mock_writer.writerow.assert_any_call([0.1, 0.01, "row1col3"])
        assert updated_headers == ["header1", "header2", "header3"]

    def test_correctly_processes_evidence_file_and_writes_updated_rows_prosit_input(
        self, mocker
    ):
        # Mock the necessary components
        mock_writer = mocker.patch("csv.writer")
        mocker.patch(
            "picked_group_fdr.parsers.tsv.get_tsv_reader",
            return_value=nullcontext(
                iter(
                    [
                        ["Header1", "Header2", "Header3"],
                        ["row1col1", "row1col2", "row1col3"],
                    ]
                )
            ),
        )
        # the first index is for the score column, the second for PEP
        mocker.patch(
            "picked_group_fdr.parsers.tsv.get_column_index", side_effect=[0, 1]
        )
        mocker.patch(
            "picked_group_fdr.parsers.maxquant.parse_evidence_file_for_percolator_matching",
            return_value=iter(
                [
                    (
                        ["row1col1", "row1col2", "row1col3"],
                        maxquant.EvidenceRow(
                            raw_file="file1.raw",
                            scannr=123,
                            peptide="_PEPTIDESEQ_",
                            score=95.5,
                            post_err_prob=0.05,
                            is_decoy=False,
                            is_contaminant=False,
                            id_type="MSMS",
                            labeling_state="NOT_SILAC",
                        ),
                    )
                ]
            ),
        )

        # Initialize the necessary variables
        evidence_file = "path/to/evidence_file.txt"
        first_headers = ["Header1", "Header2", "Header3"]
        fixed_mods = {"C": "C[UNIMOD:4]", "K": "K[UNIMOD:730]", "^_": "_[UNIMOD:730]-"}
        results_dict = {
            "file1.raw": {
                (123, "[UNIMOD:730]-PEPTIDESEQ"): (0.1, 0.01),
                (2, "peptide2"): (0.2, 0.02),
            }
        }
        pout_input_type = "prosit"
        suppress_missing_peptide_warning = False

        # Call the function
        updated_headers = update_evidence_single(
            evidence_file,
            mock_writer,
            first_headers,
            fixed_mods,
            results_dict,
            maxquant.parse_evidence_file_for_percolator_matching,
            pout_input_type,
            suppress_missing_peptide_warning,
        )

        # Assertions
        mock_writer.writerow.assert_any_call([0.1, 0.01, "row1col3"])
        assert updated_headers == ["header1", "header2", "header3"]

    # handles missing or empty evidence file gracefully
    def test_handles_missing_or_empty_evidence_file_gracefully(self, mocker):
        # Mock the necessary components
        mocker.patch("builtins.open", side_effect=FileNotFoundError)
        mock_writer = mocker.patch("csv.writer")

        # Initialize the necessary variables
        evidence_file = "path/to/non_existent_evidence_file.txt"
        first_headers = []
        fixed_mods = {"mod1": "fixed_mod1", "mod2": "fixed_mod2"}
        results_dict = {
            "result1": {(1, "peptide1"): (0.1, 0.01), (2, "peptide2"): (0.2, 0.02)}
        }
        pout_input_type = "prosit"
        suppress_missing_peptide_warning = False

        # Call the function and assert it handles the missing file gracefully
        with pytest.raises(FileNotFoundError):
            update_evidence_single(
                evidence_file,
                mock_writer,
                first_headers,
                fixed_mods,
                results_dict,
                maxquant.parse_evidence_file_for_percolator_matching,
                pout_input_type,
                suppress_missing_peptide_warning,
            )

        # Ensure no writer operations were attempted
        mock_writer.writerow.assert_not_called()
