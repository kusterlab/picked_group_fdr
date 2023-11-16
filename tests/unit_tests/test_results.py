from picked_group_fdr.results import ProteinGroupResult
from picked_group_fdr.parsers.maxquant import parse_mq_protein_groups_file_row

def test_extend():
    p = ProteinGroupResult()
    p.extend([1, 2, 3])
    assert p.extraColumns == [1, 2, 3]


def test_append():
    p = ProteinGroupResult()
    p.append(1)
    assert p.extraColumns == [1]


def test_from_mq_protein_groups():
    row = {
        "Protein IDs": "123",
        "Majority protein IDs": "456",
        "Peptide counts (unique)": "3",
        "Protein names": "protein1",
        "Gene names": "gene1",
        "Fasta headers": "header1",
        "Number of proteins": "1",
        "Q-value": "0.05",
        "Score": "10",
        "Reverse": "+",
        "Potential contaminant": "",
    }
    cols = {
        "Protein IDs": 0,
        "Majority protein IDs": 1,
        "Peptide counts (unique)": 2,
        "Protein names": 3,
        "Gene names": 4,
        "Fasta headers": 5,
        "Number of proteins": 6,
        "Q-value": 7,
        "Score": 8,
        "Reverse": 9,
        "Potential contaminant": 10,
    }
    p = parse_mq_protein_groups_file_row(list(row.values()), cols)
    assert p.proteinIds == "123"
    assert p.majorityProteinIds == "456"
    assert p.peptideCountsUnique == "3"
    assert p.proteinNames == "protein1"
    assert p.geneNames == "gene1"
    assert p.fastaHeaders == "header1"
    assert p.numberOfProteins == 1
    assert p.qValue == 0.05
    assert p.score == 10
    assert p.reverse == "+"
    assert p.potentialContaminant == ""


def test_from_protein_group():
    proteinGroup = ["protein1", "protein2"]
    peptideScores = [
        (0.05, "peptide1", ["protein1"]),
        (0.03, "peptide2", ["protein1", "protein2"]),
    ]
    reportedFdr = 0.05
    proteinScore = 10
    scoreCutoff = 0.01
    proteinAnnotations = {
        "protein1": ("protein1", "gene1", "header1"),
        "protein2": ("protein2", "gene2", "header2"),
    }
    keep_all_proteins = False
    p = ProteinGroupResult.from_protein_group(
        proteinGroup,
        peptideScores,
        reportedFdr,
        proteinScore,
        scoreCutoff,
        proteinAnnotations,
        keep_all_proteins,
    )
    assert p is None


def test_from_protein_group_keep_all_proteins():
    proteinGroup = ["protein1", "protein2"]
    peptideScores = [
        (0.05, "peptide1", ["protein1"]),
        (0.03, "peptide2", ["protein1", "protein2"]),
    ]
    reportedFdr = 0.05
    proteinScore = 10
    scoreCutoff = 0.01
    proteinAnnotations = {
        "protein1": ("protein1", "gene1", "header1"),
        "protein2": ("protein2", "gene2", "header2"),
    }
    keep_all_proteins = True
    p = ProteinGroupResult.from_protein_group(
        proteinGroup,
        peptideScores,
        reportedFdr,
        proteinScore,
        scoreCutoff,
        proteinAnnotations,
        keep_all_proteins,
    )
    assert p.proteinIds == "protein1;protein2"
    assert p.majorityProteinIds == "protein1;protein2"
    assert p.peptideCountsUnique == "0;0"
    assert p.proteinNames == "protein1;protein2"
    assert p.geneNames == "gene1;gene2"
    assert p.fastaHeaders == "header1;header2"
    assert p.numberOfProteins == 2
    assert p.qValue == 0.05
    assert p.score == 10
    assert p.reverse == ""
    assert p.potentialContaminant == ""


def test_get_peptide_counts():
    scorePeptidePairs = [
        (0.05, "peptide1", ["protein1"]),
        (0.03, "peptide2", ["protein1", "protein2"]),
        (0.02, "peptide3", ["protein2"]),
        (0.01, "peptide4", ["protein3"]),
    ]
    scoreCutoff = 0.02
    p = ProteinGroupResult()
    counts = p._get_peptide_counts(scorePeptidePairs, scoreCutoff)
    assert counts == {"protein2": 1, "protein3": 1}
