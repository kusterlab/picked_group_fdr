from __future__ import annotations

from typing import List
import logging

from . import tsv
from .. import writers
from .. import helpers

# for type hints only
from .. import scoring_strategy
from .. import results

logger = logging.getLogger(__name__)

"""psm.tsv columns:
1 Spectrum
2 Spectrum File
3 Peptide  # example: KEVPKQQAAYR (does not contain modifications)
4 Modified Peptide  # example: APHAM[147]DRDWKESVASK
5 Prev AA
6 Next AA
7 Peptide Length
8 Charge
9 Retention
10 Observed Mass
11 Calibrated Observed Mass
12 Observed M/Z
13 Calibrated Observed M/Z
14 Calculated Peptide Mass
15 Calculated M/Z
16 Delta Mass
17 SpectralSim
18 RTScore
19 Expectation
20 Hyperscore
21 Nextscore
22 PeptideProphet Probability
23 Number of Enzymatic Termini
24 Number of Missed Cleavages
25 Protein Start
26 Protein End
27 Intensity
28 Assigned Modifications
29 Observed Modifications
30 Is Unique
31 Protein
32 Protein ID
33 Entry Name
34 Gene
35 Protein Description
36 Mapped Genes
37 Mapped Proteins
"""


def parse_fragpipe_psm_file(
    reader,
    headers,
    get_proteins,
    score_type: scoring_strategy.ProteinScoringStrategy,
    **kwargs,
):
    pept_col = tsv.get_column_index(headers, "Peptide")
    modpept_col = tsv.get_column_index(headers, "Modified Peptide")
    score_col = tsv.get_column_index(
        headers, "SpectralSim"
    )  # could also use Hyperscore
    post_err_prob_col = tsv.get_column_index(headers, "PeptideProphet Probability")
    protein_col = tsv.get_column_index(headers, "Protein")
    other_proteins_col = tsv.get_column_index(headers, "Mapped Proteins")

    if score_type.get_score_column() == "pep":
        score_col = post_err_prob_col

    logger.info("Parsing FragPipe psm.tsv file")
    for line_idx, row in enumerate(reader):
        if line_idx % 500000 == 0:
            logger.info(f"    Reading line {line_idx}")

        peptide = row[pept_col]
        if len(row[modpept_col]) > 0:
            peptide = row[modpept_col]

        experiment = 1
        score = float(row[score_col])
        if score_type.get_score_column() == "pep":
            score = 1 - score + 1e-16

        proteins = [row[protein_col]]
        if len(row[other_proteins_col]) > 0:
            proteins += row[other_proteins_col].split(", ")

        proteins = get_proteins(peptide, proteins)
        if proteins:
            yield peptide, proteins, experiment, score


def parse_fragpipe_psm_file_for_peptide_remapping(reader, headers):
    protein_col = tsv.get_column_index(headers, "Protein")
    other_proteins_col = tsv.get_column_index(headers, "Mapped Proteins")

    logger.info("Parsing FragPipe psm.tsv file")
    for line_idx, row in enumerate(reader):
        if line_idx % 500000 == 0:
            logger.info(f"    Reading line {line_idx}")

        proteins = [row[protein_col]]
        if len(row[other_proteins_col]) > 0:
            proteins += row[other_proteins_col].split(", ")

        yield row, proteins


def parse_fragpipe_psm_file_for_protein_tsv(reader, headers):
    pept_col = tsv.get_column_index(headers, "Peptide")
    modpept_col = tsv.get_column_index(headers, "Modified Peptide")
    charge_col = tsv.get_column_index(headers, "Charge")
    post_err_prob_col = tsv.get_column_index(headers, "PeptideProphet Probability")
    protein_col = tsv.get_column_index(headers, "Protein")
    other_proteins_col = tsv.get_column_index(headers, "Mapped Proteins")
    assigned_mods_col = tsv.get_column_index(headers, "Assigned Modifications")
    observed_mods_col = tsv.get_column_index(headers, "Observed Modifications")

    logger.info("Parsing FragPipe psm.tsv file")
    for line_idx, row in enumerate(reader):
        if line_idx % 500000 == 0:
            logger.info(f"    Reading line {line_idx}")

        peptide = row[pept_col]
        if len(row[modpept_col]) > 0:
            peptide = row[modpept_col]
        charge = int(row[charge_col])
        post_err_prob = 1 - float(row[post_err_prob_col]) + 1e-16
        assigned_mods = row[assigned_mods_col]
        observed_mods = row[observed_mods_col]

        proteins = [row[protein_col]]
        if len(row[other_proteins_col]) > 0:
            proteins += row[other_proteins_col].split(", ")

        yield peptide, charge, post_err_prob, assigned_mods, observed_mods, proteins


def parse_fragpipe_combined_ion_file(reader, headers):
    pept_col = tsv.get_column_index(headers, "Peptide Sequence")
    modpept_col = tsv.get_column_index(headers, "Modified Sequence")
    charge_col = tsv.get_column_index(headers, "Charge")
    protein_col = tsv.get_column_index(headers, "Protein")
    other_proteins_col = tsv.get_column_index(headers, "Mapped Proteins")
    assigned_mods_col = tsv.get_column_index(headers, "Assigned Modifications")
    intensity_cols = [
        (idx, h.replace(" Intensity", ""))
        for idx, h in enumerate(headers)
        if h.endswith(" Intensity")
    ]

    logger.info("Parsing FragPipe combined_ion.tsv file")
    for line_idx, row in enumerate(reader):
        if line_idx % 100000 == 0:
            logger.info(f"    Reading line {line_idx}")

        peptide = row[pept_col]
        if len(row[modpept_col]) > 0:
            peptide = row[modpept_col]
        charge = int(row[charge_col])
        assigned_mods = row[assigned_mods_col]

        proteins = [row[protein_col]]
        if len(row[other_proteins_col]) > 0:
            proteins += row[other_proteins_col].split(", ")

        intensities = [
            (experiment, float(row[col_idx])) for col_idx, experiment in intensity_cols
        ]

        yield peptide, charge, assigned_mods, proteins, intensities


def parse_fragpipe_combined_protein_file(
    combined_protein_file: str, additional_headers: List[str] = None
) -> results.ProteinGroupResults:
    if additional_headers is None:
        additional_headers = []

    delimiter = "\t"
    if combined_protein_file.endswith(".csv"):
        delimiter = ","

    reader = tsv.get_tsv_reader(combined_protein_file, delimiter)
    headers = next(reader)

    cols = {
        x: headers.index(x)
        for x in writers.PROTEIN_GROUP_HEADERS
        + ["Protein group", "Protein Probability"]
        + additional_headers
        if x in headers
    }

    logger.info("Parsing FragPipe combined_protein_tsv file")
    protein_group_results = []
    for row in reader:
        protein_group_results.append(
            parse_fragpipe_combined_protein_file_row(row, cols, additional_headers)
        )

    protein_group_results = results.ProteinGroupResults(protein_group_results)
    protein_group_results.append_headers(additional_headers)

    return protein_group_results


def parse_fragpipe_combined_protein_file_row(row, cols, additional_headers):
    def _get_field(x, default_value=""):
        return row[cols[x]] if x in cols else default_value

    return results.ProteinGroupResult(
        proteinIds=_get_field("Protein group"),
        majorityProteinIds=_get_field("Majority protein IDs"),
        peptideCountsUnique=_get_field(
            "Peptide counts (unique)",
            ";".join(["1"] * len(_get_field("Protein group").split(";"))),
        ),
        bestPeptide="",
        numberOfProteins=int(_get_field("Number of proteins", -1)),
        qValue=float(_get_field("Q-value", 0.0)),
        score=float(_get_field("Score", _get_field("Protein Probability"))),
        reverse=_get_field(
            "Reverse",
            "+" if helpers.is_decoy(_get_field("Protein group").split(";")) else "",
        ),
        potentialContaminant=_get_field("Potential contaminant"),
        precursorQuants=[],
        extraColumns=[_get_field(x) for x in additional_headers],
    )
