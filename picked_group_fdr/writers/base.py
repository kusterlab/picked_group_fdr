from __future__ import annotations

import collections
import logging
from typing import List, Union

import numpy as np

from .. import fdr
from .. import helpers

# for type hints only
from .. import results
from .. import columns

logger = logging.getLogger(__name__)

PROTEIN_GROUP_HEADERS = [
    "Protein IDs",
    "Majority protein IDs",
    "Peptide counts (unique)",
    "Best peptide",
    "Number of proteins",
    "Q-value",
    "Score",
    "Reverse",
    "Potential contaminant",
]


def retain_only_identified_precursors(
    precursor_list: List[columns.PrecursorQuant], post_err_prob_cutoff
):
    identified_precursors = set()
    for precursor in precursor_list:
        if precursor.post_err_prob <= post_err_prob_cutoff:
            identified_precursors.add((precursor.peptide, precursor.charge))
    return [
        precursor_row
        for precursor_row in precursor_list
        if (precursor_row.peptide, precursor_row.charge) in identified_precursors
    ]


def print_num_peptides_at_fdr(post_err_probs: List, post_err_prob_cutoff: float):
    surviving_mod_peptides = set(
        [x[3] for x in post_err_probs if x[0] <= post_err_prob_cutoff]
    )

    peptides_per_rawfile = collections.defaultdict(list)
    peptides_per_experiment = collections.defaultdict(list)
    peptides_per_rawfile_mbr = collections.defaultdict(list)
    peptides_per_experiment_mbr = collections.defaultdict(list)
    for post_err_prob, rawfile, experiment, peptide in post_err_probs:
        if post_err_prob <= post_err_prob_cutoff:
            peptides_per_rawfile[rawfile].append(peptide)
            peptides_per_experiment[experiment].append(peptide)
        elif helpers.is_mbr(post_err_prob) and peptide in surviving_mod_peptides:
            peptides_per_rawfile_mbr[rawfile].append(peptide)
            peptides_per_experiment_mbr[experiment].append(peptide)

    logger.info("Precursor counts per rawfile (1% PSM-level FDR):")
    for rawfile, peptides in sorted(peptides_per_rawfile.items()):
        num_peptides = len(set(peptides))
        num_peptides_with_mbr = len(set(peptides + peptides_per_rawfile_mbr[rawfile]))
        logger.info(
            f"    {rawfile}: {num_peptides} {'(' + str(num_peptides_with_mbr) + ' with MBR)' if num_peptides_with_mbr > num_peptides else ''}"
        )

    logger.info("Precursor counts per experiment (1% PSM-level FDR):")
    for experiment, peptides in sorted(peptides_per_experiment.items()):
        num_peptides = len(set(peptides))
        num_peptides_with_mbr = len(
            set(peptides + peptides_per_experiment_mbr[experiment])
        )
        logger.info(
            f"    {experiment}: {num_peptides} {'(' + str(num_peptides_with_mbr) + ' with MBR)' if num_peptides_with_mbr > num_peptides else ''}"
        )


def append_quant_columns(
    protein_group_results: results.ProteinGroupResults,
    columns_to_add: List[columns.ProteinGroupColumns],
    post_err_probs: List,
    psm_fdr_cutoff: float,
):
    protein_group_results.remove_protein_groups_without_precursors()

    # (1) technically this is a precursor-level FDR and not a PSM-level FDR
    # (2) in contrast to MaxQuant, we set a global precursor-level FDR
    #         instead of a per raw file PSM-level FDR
    post_err_prob_cutoff = fdr.calc_post_err_prob_cutoff(
        [x[0] for x in post_err_probs if not helpers.is_mbr(x[0])], psm_fdr_cutoff
    )
    logger.info(
        f"PEP-cutoff corresponding to {psm_fdr_cutoff*100:g}% PSM-level FDR: {post_err_prob_cutoff}"
    )

    print_num_peptides_at_fdr(post_err_probs, post_err_prob_cutoff)

    logger.info("Filtering for identified precursors")
    # precursor = (peptide, charge) tuple
    # this filter also ensures that MBR precursors which were matched to
    # unidentified precursors are removed
    for pgr in protein_group_results:
        pgr.precursorQuants = retain_only_identified_precursors(
            pgr.precursorQuants, post_err_prob_cutoff
        )

    for c in columns_to_add:
        c.append(protein_group_results, post_err_prob_cutoff)

    return protein_group_results


def format_extra_columns(x: Union[str, float]) -> str:
    if type(x) == str:
        return x
    if np.isnan(x):
        return ""
    return "%.0f" % (x)
