import pytest
import numpy as np

from picked_group_fdr.columns.sequence_coverage import SequenceCoverageColumns
from picked_group_fdr.precursor_quant import PrecursorQuant


class TestGetSequenceCoverages:
    def test_get_sequence_coverages_3files(
        self, pq, peptideIntensityListThreeFiles, experimentToIdxMap
    ):
        np.testing.assert_equal(
            pq.get_sequence_coverages(
                peptideIntensityListThreeFiles, experimentToIdxMap, 0.1, "ProteinId1"
            ),
            SequenceCoverageColumns.format_as_percentage(
                [8.0 / 20.0, 8.0 / 20.0, 8.0 / 20.0, 8.0 / 20.0, 8.0 / 20.0, 8.0 / 20.0]
            ),
        )

    def test_get_sequence_coverages_2peptides(
        self, pq, peptideIntensityListMultiplePeptides, experimentToIdxMap
    ):
        np.testing.assert_equal(
            pq.get_sequence_coverages(
                peptideIntensityListMultiplePeptides,
                experimentToIdxMap,
                0.1,
                "ProteinId1",
            ),
            SequenceCoverageColumns.format_as_percentage(
                [
                    16.0 / 20.0,
                    16.0 / 20.0,
                    16.0 / 20.0,
                    16.0 / 20.0,
                    16.0 / 20.0,
                    16.0 / 20.0,
                ]
            ),
        )

    def test_get_sequence_coverages_4peptides(
        self,
        pq,
        peptideIntensityListMultiplePeptidesSequenceCoverage,
        experimentToIdxMap,
    ):
        np.testing.assert_equal(
            pq.get_sequence_coverages(
                peptideIntensityListMultiplePeptidesSequenceCoverage,
                experimentToIdxMap,
                0.1,
                "ProteinId1",
            ),
            SequenceCoverageColumns.format_as_percentage(
                [
                    20.0 / 20.0,
                    20.0 / 20.0,
                    20.0 / 20.0,
                    17.0 / 20.0,
                    18.0 / 20.0,
                    12.0 / 20.0,
                ]
            ),
        )


@pytest.fixture
def pq(proteinSequences):
    return SequenceCoverageColumns(proteinSequences)


@pytest.fixture
def proteinSequences():
    return {"ProteinId1": "APEPTIDEAAAABPEPTIDE"}


@pytest.fixture
def peptideIntensityListMultiplePeptidesSequenceCoverage():
    peptideIntensityList = list()
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE(ox)", 2, "file1", 1, 25.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file2", 1, 10.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file3", 1, 5.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("ABPEPTIDE", 2, "file1", 1, 15.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("A(ox)ABPEPTIDE", 2, "file2", 1, 0.0, 0.001, [], [], 1)
    )  # only identified by MS/MS; no MS1 feature available
    peptideIntensityList.append(
        PrecursorQuant("DEAAAA", 2, "file3", 1, 5.0, 0.001, [], [], 1)
    )  # tests overlap

    return peptideIntensityList
