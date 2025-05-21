# Generation assisted by CodiumAI

from picked_group_fdr.parsers.fragpipe import parse_fragpipe_psm_file
from picked_group_fdr.scoring_strategy import ProteinScoringStrategy

import pytest
from contextlib import nullcontext


class TestParseFragpipePsmFile:

    # correctly parses a well-formed psm.tsv file with all required columns
    def test_correctly_parses_well_formed_psm_tsv(self, mocker):
        # Mock the csv reader
        mocker.patch(
            "picked_group_fdr.parsers.tsv.get_tsv_reader",
            return_value=nullcontext(
                iter(
                    [
                        [
                            "Spectrum",
                            "Peptide",
                            "Modified Peptide",
                            "SpectralSim",
                            "PeptideProphet Probability",
                            "Protein",
                            "Mapped Proteins",
                        ],
                        [
                            "Spectrum1",
                            "Peptide1",
                            "ModifiedPeptide1",
                            "0.95",
                            "0.99",
                            "Protein1",
                            "Protein2, Protein3",
                        ],
                        ["Spectrum2", "Peptide2", "", "0.85", "0.95", "Protein4", ""],
                    ]
                )
            ),
        )

        def get_proteins(peptide, proteins):
            return proteins

        evidence_file = "psm.tsv"
        score_column = "pep"

        results = list(parse_fragpipe_psm_file(evidence_file, get_proteins, score_column))

        assert results == [
            (
                "ModifiedPeptide1",
                ["Protein1", "Protein2", "Protein3"],
                1,
                pytest.approx(0.01),
            ),
            ("Peptide2", ["Protein4"], 1, pytest.approx(0.05)),
        ]

    # handles missing "Mapped Proteins" column gracefully
    def test_handles_missing_mapped_proteins_column(self, mocker):
        # Mock the csv reader
        mocker.patch(
            "picked_group_fdr.parsers.tsv.get_tsv_reader",
            return_value=nullcontext(
                iter(
                    [
                        [
                            "Spectrum",
                            "Peptide",
                            "Modified Peptide",
                            "SpectralSim",
                            "PeptideProphet Probability",
                            "Protein",
                        ],
                        [
                            "Spectrum1",
                            "Peptide1",
                            "ModifiedPeptide1",
                            "0.95",
                            "0.99",
                            "Protein1",
                        ],
                        ["Spectrum2", "Peptide2", "", "0.85", "0.95", "Protein4"],
                    ]
                )
            ),
        )

        def get_proteins(peptide, proteins):
            return proteins

        evidence_file = "psm.tsv"
        score_type = ProteinScoringStrategy("bestPEP")

        with pytest.raises(ValueError):
            list(parse_fragpipe_psm_file(evidence_file, get_proteins, score_type))
