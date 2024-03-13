from __future__ import annotations

import logging
from typing import List, Optional

import numpy as np
import pandas as pd

from .. import helpers
from ..parsers import psm
from ..parsers import parsers
from ..precursor_quant import PrecursorQuant

# for type hints only
from .. import digest
from .. import results
from .. import protein_groups as pg
from .. import scoring_strategy

logger = logging.getLogger(__name__)


def add_precursor_quants(
    mq_evidence_files: List[str],
    mq_quantification_files: List[str],
    protein_groups: pg.ProteinGroups,
    protein_group_results: results.ProteinGroupResults,
    peptide_to_protein_maps: List[digest.PeptideToProteinMap],
    experimental_design: Optional[pd.DataFrame],
    discard_shared_peptides: bool,
    score_type: scoring_strategy.ProteinScoringStrategy,
    suppress_missing_peptide_warning: bool,
):
    file_mapping = None
    if experimental_design is not None:
        protein_group_results.experiments = experimental_design["Experiment"].unique().tolist()
        file_mapping = parsers.get_file_mapping(experimental_design)

    post_err_probs = list()
    shared_peptide_precursors, unique_peptide_precursors = 0, 0
    parsed_experiments = set()
    missing_peptides_in_protein_groups = 0

    for (
        peptide,
        proteins,
        charge,
        raw_file,
        experiment,
        fraction,
        intensity,
        post_err_prob,
        tmt_cols,
        silac_cols,
        evidence_id,
    ) in psm.parse_evidence_file_multiple(
        mq_evidence_files,
        peptide_to_protein_maps=peptide_to_protein_maps,
        score_type=score_type,
        for_quantification=True,
    ):
        if protein_group_results.num_tmt_channels == -1:
            # There are 3 columns per TMT channel:
            #     Reporter intensity corrected,
            #     Reporter intensity
            #     Reporter intensity count
            protein_group_results.num_tmt_channels = int(len(tmt_cols) / 3)
        if protein_group_results.num_silac_channels == -1:
            protein_group_results.num_silac_channels = len(silac_cols)

        # override the parsed experiment and fraction if --file_list_file option is used
        if file_mapping:
            experiment, fraction = file_mapping[raw_file]
        elif experiment not in parsed_experiments:
            parsed_experiments.add(experiment)

        protein_group_idxs = protein_groups.get_protein_group_idxs(proteins)

        # removes peptides not present in the proteinGroups.txt file
        if helpers.is_missing_in_protein_groups(protein_group_idxs):
            if not suppress_missing_peptide_warning:
                logger.debug(
                    f"Could not find any of the proteins {proteins} in proteinGroups.txt"
                )
            missing_peptides_in_protein_groups += 1
            continue

        if discard_shared_peptides and helpers.is_shared_peptide(protein_group_idxs):
            shared_peptide_precursors += 1
            continue

        unique_peptide_precursors += 1

        if not helpers.is_decoy(proteins):
            post_err_probs.append((post_err_prob, raw_file, experiment, peptide))

        if len(tmt_cols) > 0:
            tmt_cols = np.array(tmt_cols, dtype="float64")
        else:
            tmt_cols = None

        if len(silac_cols) > 0:
            silac_cols = np.array(silac_cols, dtype="float64")
        else:
            silac_cols = None

        for protein_group_idx in protein_group_idxs:
            precursor_quant = PrecursorQuant(
                peptide,
                charge,
                experiment,
                fraction,
                intensity,
                post_err_prob,
                tmt_cols,
                silac_cols,
                evidence_id,
            )
            protein_group_results[protein_group_idx].precursorQuants.append(
                precursor_quant
            )

    if missing_peptides_in_protein_groups > 0:
        logger.debug(
            f"Skipped {missing_peptides_in_protein_groups} precursors from proteins not present in proteinGroups.txt file"
        )

    logger.info(
        f"Found {unique_peptide_precursors} precursors from unique and {shared_peptide_precursors} precursors from shared peptides"
    )

    if len(parsed_experiments) > 0:
        protein_group_results.experiments = sorted(list(parsed_experiments))

    return protein_group_results, post_err_probs
