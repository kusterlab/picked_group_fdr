from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import numpy as np

from . import tsv

# for type hints only
from .. import scoring_strategy

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

lfq.tsv columns:
1 peptide
2 charge
3 proteins
4 q_value
5 score
6 spectral_angle
7 <file_1>
...
"""


def parse_sage_results_file(
    reader,
    headers,
    get_proteins,
    score_type: Optional[scoring_strategy.ProteinScoringStrategy],
    for_quantification: bool = False,
):
    use_post_err_prob = score_type is None or score_type.get_score_column() == "pep"

    pept_col = tsv.get_column_index(headers, "peptide")
    charge_col = tsv.get_column_index(headers, "charge")
    score_col = tsv.get_column_index(headers, "sage_discriminant_score")
    filename_col = tsv.get_column_index(headers, "filename")
    post_err_prob_col = tsv.get_column_index(headers, "posterior_error")
    protein_col = tsv.get_column_index(headers, "proteins")

    if use_post_err_prob:
        score_col = post_err_prob_col

    logger.info("Parsing Sage results.sage.tsv file")
    for line_idx, row in enumerate(reader):
        if line_idx % 500000 == 0:
            logger.info(f"    Reading line {line_idx}")

        modified_peptide = row[pept_col] # example: FC[+57.0215]LPYRMDVEK
        charge = int(row[charge_col])
        filename = Path(row[filename_col]).stem
        score = float(row[score_col])
        if use_post_err_prob:
            # sage's posterior_error column is log10(PEP) transformed
            score = np.power(10, score)

        proteins = row[protein_col].split(";")

        proteins = get_proteins(modified_peptide, proteins)
        if proteins:
            if for_quantification:
                yield modified_peptide, proteins, filename, score, charge
            else:
                yield modified_peptide, proteins, filename, score


def parse_sage_lfq_file(reader, headers):
    pept_col = tsv.get_column_index(headers, "peptide")
    charge_col = tsv.get_column_index(headers, "charge")
    protein_col = tsv.get_column_index(headers, "proteins")

    # here we make the dangerous assumption that the columns after spectral_angle are
    # the intensity columns. For now, there is no better way, since there is no other
    # way to identify the intensity columns.
    spectral_angle_col = tsv.get_column_index(headers, "spectral_angle")
    intensity_cols = [
        (idx, Path(h).stem) for idx, h in enumerate(headers) if idx > spectral_angle_col
    ]

    logger.info("Parsing Sage lfq.tsv file")
    for line_idx, row in enumerate(reader):
        if line_idx % 100000 == 0:
            logger.info(f"    Reading line {line_idx}")

        peptide = row[pept_col]
        charge = int(row[charge_col])

        proteins = row[protein_col].split(";")

        intensities = [
            (filename, float(row[col_idx])) for col_idx, filename in intensity_cols
        ]

        yield peptide, charge, proteins, intensities
