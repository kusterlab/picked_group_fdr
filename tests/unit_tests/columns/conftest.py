import pytest

import numpy as np

from picked_group_fdr.precursor_quant import PrecursorQuant


@pytest.fixture
def experimentToIdxMap():
    experiments = ["file1", "file2", "file3"]
    experimentToIdxMap = dict([(v, k) for k, v in enumerate(experiments)])
    return experimentToIdxMap


@pytest.fixture
def experimentToIdxMapFiveFiles():
    experiments = ["file1", "file2", "file3", "file4", "file5"]
    experimentToIdxMap = dict([(v, k) for k, v in enumerate(experiments)])
    return experimentToIdxMap


@pytest.fixture
def peptideIntensityListThreeFiles():
    peptideIntensityList = list()
    # peptide, charge, experiment, fraction, intensity, postErrProb, tmtIntensities, silacIntensities
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file1", 1, 25.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file2", 1, 10.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file3", 1, 5.0, 0.001, [], [], 1)
    )

    return peptideIntensityList


@pytest.fixture
def peptideIntensityListFiveFiles():
    peptideIntensityList = list()
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file1", 1, 25.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file2", 1, 10.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file3", 1, 5.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file4", 1, 0.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file5", 1, 7.0, 0.001, [], [], 1)
    )

    return peptideIntensityList


@pytest.fixture
def peptideIntensityListMissingValuesZeroes():
    peptideIntensityList = list()
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file1", 1, 25.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file2", 1, 0.0, 0.001, [], [], 1)
    )  # only identified by MS/MS; no MS1 feature available
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file3", 1, 0.0, 0.001, [], [], 1)
    )

    return peptideIntensityList


@pytest.fixture
def peptideIntensityListMissingValuesNans():
    peptideIntensityList = list()
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file1", 1, 25.0, 0.001, [], [], 1)
    )

    return peptideIntensityList


@pytest.fixture
def peptideIntensityListMultipleCharges():
    peptideIntensityList = list()
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file1", 1, 25.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file2", 1, 10.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file3", 1, 5.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 3, "file1", 1, 15.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 3, "file2", 1, 10.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 3, "file3", 1, 5.0, 0.001, [], [], 1)
    )

    return peptideIntensityList


@pytest.fixture
def peptideIntensityListMultiplePeptides():
    peptideIntensityList = list()
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file1", 1, 25.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file2", 1, 10.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file3", 1, 5.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("BPEPTIDE", 2, "file1", 1, 15.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("BPEPTIDE", 2, "file2", 1, 10.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("BPEPTIDE", 2, "file3", 1, 5.0, 0.001, [], [], 1)
    )

    return peptideIntensityList


@pytest.fixture
def peptideIntensityListMultiplePeptidesJagged():
    peptideIntensityList = list()
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file1", 1, 25.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file2", 1, 10.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file3", 1, 28.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("BPEPTIDE", 2, "file1", 1, 5.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("BPEPTIDE", 2, "file2", 1, 2.0, 0.001, [], [], 1)
    )

    return peptideIntensityList


@pytest.fixture
def peptideIntensityListDuplicatePeptides():
    peptideIntensityList = list()
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file1", 1, 25.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file1", 1, 10.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file2", 1, 10.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file3", 1, 5.0, 0.001, [], [], 1)
    )

    return peptideIntensityList


@pytest.fixture
def peptideIntensityListMultipleFractions():
    peptideIntensityList = list()
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file1", 1, 25.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file1", 2, 10.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file2", 1, 10.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file3", 1, 5.0, 0.001, [], [], 1)
    )

    return peptideIntensityList


@pytest.fixture
def peptideIntensityListMatchBetweenRuns():
    peptideIntensityList = list()
    # peptide, charge, experiment, fraction, intensity, postErrProb, tmtIntensities, silacIntensities
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file1", 1, 25.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file2", 1, 10.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file3", 1, 5.0, np.nan, [], [], 1)
    )

    return peptideIntensityList


@pytest.fixture
def peptideIntensityListMatchBetweenRunsWithUnidentified():
    peptideIntensityList = list()
    # peptide, charge, experiment, fraction, intensity, postErrProb, tmtIntensities, silacIntensities
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file1", 1, 25.0, 0.2, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file2", 1, 10.0, 0.2, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file3", 1, 5.0, np.nan, [], [], 1)
    )

    return peptideIntensityList


@pytest.fixture
def peptideIntensityListTMT():
    peptideIntensityList = list()
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file1", 1, 25.0, 0.001, [1.0, 2.0, 3.0], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file1", 2, 10.0, 0.001, [2.0, 5.0, 1.0], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file2", 1, 10.0, 0.001, [3.0, 4.0, 5.0], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file3", 1, 5.0, 0.001, [5.0, 6.0, 7.0], [], 1)
    )

    return peptideIntensityList


@pytest.fixture
def peptideIntensityListSilac():
    peptideIntensityList = list()
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file1", 1, 25.0, 0.001, [], [11.0, 12.0, 2.0], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file1", 2, 10.0, 0.001, [], [2.0, 5.0, 3.0], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file2", 1, 10.0, 0.001, [], [3.0, 4.0, 3.0], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("APEPTIDE", 2, "file3", 1, 5.0, 0.001, [], [1.5, 2.5, 1.0], 1)
    )

    return peptideIntensityList


@pytest.fixture
def peptideIntensityListLargeRatioStabilizationHighRatio():
    """
    Test fixture for large ratio stabilization where peptide_count_ratio > 5.
    file1: 1 peptide (PEPTIDEA) with intensity 100
    file2: 10+ peptides with mixed intensities summing to 100
    This simulates a situation where file1 has few peptides and file2 has many.
    """
    peptideIntensityList = list()
    # file1: 1 peptide
    peptideIntensityList.append(
        PrecursorQuant("PEPTIDEA", 2, "file1", 1, 100.0, 0.001, [], [], 1)
    )

    # file2: 10 peptides with varying intensities (sum = 100)
    peptides_file2 = [
        "PEPTIDEA",
        "PEPTIDEB",
        "PEPTIDEC",
        "PEPTIDED",
        "PEPTIDEE",
        "PEPTIDEF",
        "PEPTIDEG",
        "PEPTIDEH",
        "PEPTIDEI",
        "PEPTIDEJ",
    ]
    for pep in peptides_file2:
        # Create peptides with intensities that sum to 100
        peptideIntensityList.append(
            PrecursorQuant(pep, 2, "file2", 1, 10.0, 0.001, [], [], 1)
        )

    return peptideIntensityList


@pytest.fixture
def peptideIntensityListLargeRatioStabilizationMediumRatio():
    """
    Test fixture for medium ratio stabilization where 2.5 < peptide_count_ratio <= 5.
    file1: 2 peptides with intensity 100
    file2: 8 peptides with mixed intensities summing to 100
    Ratio = 8/2 = 4.0, which should trigger weighted average
    """
    peptideIntensityList = list()
    # file1: 2 peptides
    peptideIntensityList.append(
        PrecursorQuant("PEPTIDEA", 2, "file1", 1, 50.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("PEPTIDEB", 2, "file1", 1, 50.0, 0.001, [], [], 1)
    )

    # file2: 8 peptides with varying intensities (sum = 100)
    peptides_file2 = [
        "PEPTIDEA",
        "PEPTIDEB",
        "PEPTIDEC",
        "PEPTIDED",
        "PEPTIDEE",
        "PEPTIDEF",
        "PEPTIDEG",
        "PEPTIDEH",
    ]
    for pep in peptides_file2:
        peptideIntensityList.append(
            PrecursorQuant(pep, 2, "file2", 1, 12.5, 0.001, [], [], 1)
        )

    return peptideIntensityList


@pytest.fixture
def peptideIntensityListLargeRatioStabilizationSmallRatio():
    """
    Test fixture for small ratio where peptide_count_ratio <= 2.5.
    file1: 4 peptides with intensity 100
    file2: 8 peptides with mixed intensities summing to 100
    Ratio = 8/4 = 2.0, which should NOT trigger stabilization
    """
    peptideIntensityList = list()
    # file1: 4 peptides
    for i in range(4):
        peptideIntensityList.append(
            PrecursorQuant(f"APEPTIDE{i}", 2, "file1", 1, 25.0, 0.001, [], [], 1)
        )

    # file2: 8 peptides with varying intensities (sum = 100)
    for i in range(8):
        peptideIntensityList.append(
            PrecursorQuant(f"BPEPTIDE{i}", 2, "file2", 1, 12.5, 0.001, [], [], 1)
        )

    return peptideIntensityList


@pytest.fixture
def peptideIntensityListLargeRatioStabilizationZeroCount():
    """
    Test fixture for zero peptide count scenario.
    file1: 1 peptide
    file2: 0 peptides (no PSMs)
    file3: 1 peptide
    """
    peptideIntensityList = list()
    peptideIntensityList.append(
        PrecursorQuant("PEPTIDEA", 2, "file1", 1, 100.0, 0.001, [], [], 1)
    )
    peptideIntensityList.append(
        PrecursorQuant("PEPTIDEC", 2, "file3", 1, 100.0, 0.001, [], [], 1)
    )

    return peptideIntensityList
