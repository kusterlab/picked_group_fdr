import logging
from typing import List, Optional
from datetime import datetime

import numpy as np

from . import helpers
from . import entrapment


logger = logging.getLogger(__name__)


def calculate_protein_fdrs(
    protein_groups, protein_scores, protein_group_fdr_threshold: float
):
    logger.info("Calculating protein group-level FDRs")
    num_decoys, num_entrapments, num_targets = 0, 0, 0
    protein_group_info_list = list()
    for protein_group, protein_score in zip(protein_groups, protein_scores):
        if protein_score == -100.0:
            break

        if helpers.is_decoy(protein_group):
            num_decoys += 1
        else:
            num_targets += 1
            if entrapment.is_entrapment(protein_group):
                num_entrapments += 1
        reported_fdr = (num_decoys + 1) / (num_targets + 1)
        observed_fdr = (num_entrapments + 1) / (num_targets + 1)

        skip_for_counting = helpers.is_decoy(protein_group) or helpers.is_obsolete(
            protein_group
        )
        protein_group_info_list.append((reported_fdr, observed_fdr, skip_for_counting))

    logger.info(f"#Protein groups: Target = {num_targets}; Decoys = {num_decoys}")
    if num_entrapments > 1:
        logger.info(
            f"    Entrapments = {num_entrapments}; Targets-Entrapments = {num_targets - num_entrapments}"
        )

    if len(protein_group_info_list) == 0:
        raise Exception(
            "No proteins with scores found, make sure that protein identifiers are consistent in the evidence and fasta files"
        )

    reported_fdrs, observed_fdrs, skip_for_counting = zip(*protein_group_info_list)
    reported_qvals, observed_qvals = fdrs_to_qvals(reported_fdrs), fdrs_to_qvals(
        observed_fdrs
    )
    logger.info(
        f"#Target protein groups at {protein_group_fdr_threshold*100:g}% decoy FDR: {count_below_threshold(reported_qvals, protein_group_fdr_threshold, skip_for_counting)}"
    )
    if num_entrapments > 1:
        logger.info(
            f"#Target protein groups at {protein_group_fdr_threshold*100:g}% entrapment FDR: {count_below_threshold(observed_fdrs, protein_group_fdr_threshold, skip_for_counting)}"
        )
        num_decoys_at_entrapment_fdr = count_below_threshold(
            observed_fdrs, protein_group_fdr_threshold
        )
        reported_qvals_index = min(
            len(reported_qvals) - 1, num_decoys_at_entrapment_fdr
        )
        logger.info(
            f"Decoy FDR at {protein_group_fdr_threshold*100:g}% entrapment FDR: {'%.2g' % (reported_qvals[reported_qvals_index])}"
        )
        num_entrapments_at_decoy_fdr = count_below_threshold(
            reported_qvals, protein_group_fdr_threshold
        )
        observed_fdrs_index = min(len(observed_fdrs) - 1, num_entrapments_at_decoy_fdr)
        logger.info(
            f"Entrapment FDR at {protein_group_fdr_threshold*100:g}% decoy FDR: {'%.2g' % (observed_fdrs[observed_fdrs_index])}"
        )

        # write_reported_and_entrapment_fdrs(reported_qvals, observed_qvals)

    return reported_qvals, observed_qvals


def calculate_peptide_fdrs(peptide_scores, score_type, peptide_fdr_threshold: float):
    decoy_scores, entrapment_scores, pool_scores = list(), list(), list()
    fdrs = list()
    reported_fdr, observed_fdr = 0, 0
    sum_post_err_prob = np.nextafter(0, 1)
    for score, _, proteins in peptide_scores:
        if helpers.is_contaminant(proteins):
            continue

        if helpers.is_decoy(proteins):
            decoy_scores.append(score)
        else:
            if entrapment.is_entrapment(proteins):
                entrapment_scores.append(score)
            else:
                pool_scores.append(score)
            if "PEPavg" in score_type.long_description():
                sum_post_err_prob += score
                reported_fdr = sum_post_err_prob / (
                    len(entrapment_scores) + len(pool_scores)
                )
            else:
                reported_fdr = (len(decoy_scores) + 1) / (
                    len(pool_scores) + len(entrapment_scores) + 1
                )
            observed_fdr = (len(entrapment_scores) + 1) / (
                len(pool_scores) + len(entrapment_scores) + 1
            )
            fdrs.append((reported_fdr, observed_fdr))
        # print(score, reported_fdr, observed_fdr, helpers.is_decoy(proteins), entrapment.is_entrapment(proteins), sep='\t')

    num_decoys = len(decoy_scores)
    num_entrapments = len(entrapment_scores)
    num_targets = len(pool_scores) + len(entrapment_scores)
    logger.info(f"#Peptides: Target = {num_targets}; Decoys = {num_decoys}")
    if num_entrapments > 1:
        logger.info(
            f"    Entrapments = {num_entrapments}; Targets-Entrapments = {num_targets - num_entrapments}"
        )

    reported_fdrs, observed_fdrs = zip(*fdrs)
    reported_qvals, observed_qvals = fdrs_to_qvals(reported_fdrs), fdrs_to_qvals(
        observed_fdrs
    )
    logger.info(
        f"#Target peptides at {peptide_fdr_threshold*100:g}% decoy FDR: {count_below_threshold(reported_qvals, peptide_fdr_threshold)}"
    )
    if num_entrapments > 1:
        logger.info(
            f"#Target peptides at {peptide_fdr_threshold*100:g}% entrapment FDR: {count_below_threshold(observed_fdrs, peptide_fdr_threshold)}"
        )
        num_decoys_at_entrapment_fdr = count_below_threshold(
            observed_fdrs, peptide_fdr_threshold
        )
        reported_qvals_index = min(
            len(reported_qvals) - 1, num_decoys_at_entrapment_fdr
        )
        logger.info(
            f"Decoy FDR at {peptide_fdr_threshold*100:g}% entrapment FDR: {'%.2g' % (reported_qvals[reported_qvals_index])}"
        )
        num_entrapments_at_decoy_fdr = count_below_threshold(
            reported_qvals, peptide_fdr_threshold
        )
        observed_fdrs_index = min(len(observed_fdrs) - 1, num_entrapments_at_decoy_fdr)
        logger.info(
            f"Entrapment FDR at {peptide_fdr_threshold*100:g}% decoy FDR: {'%.2g' % (observed_fdrs[observed_fdrs_index])}"
        )

    return reported_qvals, observed_qvals


def write_reported_and_entrapment_fdrs(reported_qvals, observed_qvals):
    import csv

    writer = csv.writer(
        open(
            f'protein_fdr_calibration_{datetime.now().strftime("%d%m%Y_%H%M%S")}.txt',
            "w",
        ),
        delimiter="\t",
    )
    for reported_qval, observed_qval in zip(reported_qvals, observed_qvals):
        writer.writerow([reported_qval, observed_qval])


def fdrs_to_qvals(fdrs: List[float]) -> np.array:
    """
    Makes a list of FDRs monotonically increasing (sometimes referred to as q-values after monotonization)
    """
    return np.minimum.accumulate(fdrs[::-1])[::-1]


def count_below_threshold(
    qvals: List[float],
    qval_threshold: float,
    skip_for_counting: Optional[List[bool]] = None,
):
    """
    Counts number of q-values below a threshold, if skip_for_counting are provided, only the targets are counted
    """
    if skip_for_counting is None:
        return len([1 for x in qvals if x < qval_threshold])
    else:
        return len(
            [
                1
                for x, skip in zip(qvals, skip_for_counting)
                if x < qval_threshold and not skip
            ]
        )


def calc_post_err_prob_cutoff(post_err_probs, psm_qval_cutoff):
    post_err_prob_cutoff = 1.0
    sum_post_err_probs = 0.0
    num_psms = 0
    for post_err_prob in sorted(post_err_probs):
        if not np.isfinite(post_err_prob):
            continue
        sum_post_err_probs += post_err_prob
        num_psms += 1
        if sum_post_err_probs / num_psms > psm_qval_cutoff:
            post_err_prob_cutoff = post_err_prob
            break
    return post_err_prob_cutoff
