from __future__ import annotations

import logging

import numpy as np

from . import tsv

# for type hints only
from .. import scoring

logger = logging.getLogger(__name__)

"""https://sage-docs.vercel.app/docs/results/search

results.sage.tsv columns:
1 peptide
2 proteins
3 num_proteins
4 filename
5 scannr
6 rank
7 label
8 expmass
9 calcmass
10 charge
11 peptide_len
12 missed_cleavages
13 isotope_error
14 precursor_ppm
15 fragment_ppm
16 hyperscore
17 delta_next
18 delta_best
19 rt
20 aligned_rt
21 predicted_rt
22 delta_rt_model
23 matched_peaks
24 longest_b
25 longest_y
26 longest_y_pct
27 matched_intensity_pct
28 scored_candidates
29 poisson
30 sage_discriminant_score
31 posterior_error
32 spectrum_q
33 peptide_q
34 protein_q
35 ms1_intensity
36 ms2_intensity
"""


def is_sage_results_file(headers):
    return "sage_discriminant_score" in map(str.lower, headers)


def parse_sage_results_file(
    reader,
    headers,
    get_proteins,
    score_type: scoring.ProteinScoringStrategy,
    **kwargs,
):
    pept_col = tsv.get_column_index(headers, "peptide")
    score_col = tsv.get_column_index(
        headers, "sage_discriminant_score"
    )  # could also use Hyperscore
    post_err_prob_col = tsv.get_column_index(headers, "posterior_error")
    protein_col = tsv.get_column_index(headers, "proteins")

    if score_type.get_score_column() == "pep":
        score_col = post_err_prob_col

    logger.info("Parsing Sage results.sage.tsv file")
    for line_idx, row in enumerate(reader):
        if line_idx % 500000 == 0:
            logger.info(f"    Reading line {line_idx}")

        peptide = row[pept_col]
        experiment = 1
        score = float(row[score_col])
        if score_type.get_score_column() == "pep":
            score = np.power(10, score)

        proteins = row[protein_col].split(";")

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
        charge = int(row[charge_col])
        assigned_mods = row[assigned_mods_col]

        proteins = [row[protein_col]]
        if len(row[other_proteins_col]) > 0:
            proteins += row[other_proteins_col].split(", ")

        intensities = [
            (experiment, float(row[col_idx])) for col_idx, experiment in intensity_cols
        ]

        yield peptide, charge, assigned_mods, proteins, intensities
