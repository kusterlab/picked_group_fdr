from __future__ import annotations

import collections
import logging
from typing import Dict, List, Union

import numpy as np

from .. import fdr
from .. import helpers

# for type hints only
from .. import results
from .. import columns
from ..precursor_quant import PrecursorQuant

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


class ProteinGroupsWriter:
    def get_header_dict(
        self, protein_group_results: results.ProteinGroupResults
    ) -> Dict[str, str]:
        return {x: x for x in protein_group_results.headers}

    def get_columns(self) -> List[columns.ProteinGroupColumns]:
        pass

    def get_extra_columns_formatter(self):
        return format_extra_columns

    def append_quant_columns(
        self,
        protein_group_results: results.ProteinGroupResults,
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

        _print_num_peptides_at_fdr(post_err_probs, post_err_prob_cutoff)

        logger.info("Filtering for identified precursors")
        # precursor = (peptide, charge) tuple
        # this filter also ensures that MBR precursors which were matched to
        # unidentified precursors are removed
        for pgr in protein_group_results:
            pgr.precursorQuants = _retain_only_identified_precursors(
                pgr.precursorQuants, post_err_prob_cutoff
            )

        for c in self.get_columns():
            c.append(protein_group_results, post_err_prob_cutoff)

        return protein_group_results

    def write(
        self,
        protein_group_results: results.ProteinGroupResults,
        protein_groups_out_file: str,
    ):
        protein_group_results.write(
            protein_groups_out_file,
            header_dict=self.get_header_dict(protein_group_results),
            format_extra_columns=self.get_extra_columns_formatter(),
        )


def _retain_only_identified_precursors(
    precursor_list: List[PrecursorQuant], post_err_prob_cutoff
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


def _print_num_peptides_at_fdr(post_err_probs: List, post_err_prob_cutoff: float):
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


def format_extra_columns(x: Union[str, float]) -> str:
    if type(x) == str:
        return x
    if np.isnan(x):
        return ""
    return "%.0f" % (x)


