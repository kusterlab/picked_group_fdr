from __future__ import annotations

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from .. import helpers
from ..parsers import sage, tsv, parsers
from ..precursor_quant import PrecursorQuant

# for type hints only
from .. import digest
from .. import results
from .. import protein_groups as pg
from .. import scoring_strategy

logger = logging.getLogger(__name__)


def add_precursor_quants(
    fragpipe_psm_file: str,
    protein_group_results: results.ProteinGroupResults,
    protein_groups: pg.ProteinGroups,
    experimental_design: Optional[pd.DataFrame],
    discard_shared_peptides: bool,
    suppress_missing_peptide_warning: bool,
):
    file_mapping = None
    if experimental_design is not None:
        protein_group_results.experiments = (
            experimental_design["Experiment"].unique().tolist()
        )
        file_mapping = parsers.get_file_mapping(experimental_design)

    delimiter = tsv.get_delimiter(fragpipe_psm_file)
    reader = tsv.get_tsv_reader(fragpipe_psm_file, delimiter)
    headers = next(reader)

    get_proteins = lambda peptide, proteins: proteins

    post_err_probs = []
    parsed_experiments = set()
    for (
        peptide,
        proteins,
        filename,
        post_err_prob,
        charge,
    ) in sage.parse_sage_results_file(
        reader, headers, get_proteins, score_type=None, for_quantification=True
    ):
        protein_group_idxs = protein_groups.get_protein_group_idxs(proteins)

        if helpers.is_missing_in_protein_groups(protein_group_idxs):
            if not suppress_missing_peptide_warning:
                logger.debug(
                    f"Could not find any of the proteins {proteins} in proteinGroups.txt"
                )
            continue

        if discard_shared_peptides and helpers.is_shared_peptide(protein_group_idxs):
            continue

        experiment, fraction = filename, -1
        if file_mapping:
            experiment, fraction = file_mapping[filename]
        elif experiment not in parsed_experiments:
            parsed_experiments.add(experiment)

        if not helpers.is_decoy(proteins):
            post_err_probs.append((post_err_prob, filename, experiment, peptide))

        for protein_group_idx in protein_group_idxs:
            precursorQuant = PrecursorQuant(
                peptide=peptide,
                charge=charge,
                experiment=experiment,
                fraction=fraction,
                intensity=np.nan,
                post_err_prob=post_err_prob,
                tmt_intensities=None,  # TODO: add TMT support
                silac_intensities=None,  # TODO: add SILAC support
                evidence_id=-1,
            )
            protein_group_results[protein_group_idx].precursorQuants.append(
                precursorQuant
            )

    if len(parsed_experiments) > 0:
        protein_group_results.experiments.extend(sorted(list(parsed_experiments)))

    return protein_group_results, post_err_probs


def update_precursor_quants(
    protein_group_results: results.ProteinGroupResults,
    protein_groups: pg.ProteinGroups,
    sage_lfq_tsvs: List[str],
    experimental_design: Optional[pd.DataFrame],
    discard_shared_peptides: bool,
    suppress_missing_peptide_warning: bool,
):
    for sage_lfq_tsv in sage_lfq_tsvs:
        protein_group_results = update_precursor_quants_single(
            protein_group_results,
            protein_groups,
            sage_lfq_tsv,
            experimental_design,
            discard_shared_peptides,
            suppress_missing_peptide_warning,
        )
    return protein_group_results


def update_precursor_quants_single(
    protein_group_results: results.ProteinGroupResults,
    protein_groups: pg.ProteinGroups,
    sage_lfq_tsv: str,
    experimental_design: Optional[pd.DataFrame],
    discard_shared_peptides: bool,
    suppress_missing_peptide_warning: bool,
):
    file_mapping = None
    if experimental_design is not None:
        file_mapping = parsers.get_file_mapping(experimental_design)

    delimiter = tsv.get_delimiter(sage_lfq_tsv)
    reader = tsv.get_tsv_reader(sage_lfq_tsv, delimiter)
    headers = next(reader)

    for (
        peptide,
        charge,
        proteins,
        intensities,
    ) in sage.parse_sage_lfq_file(reader, headers):
        protein_group_idxs = protein_groups.get_protein_group_idxs(proteins)

        if helpers.is_missing_in_protein_groups(protein_group_idxs):
            if not suppress_missing_peptide_warning:
                logger.debug(
                    f"Could not find any of the proteins {proteins} in proteinGroups.txt"
                )
            continue

        if discard_shared_peptides and helpers.is_shared_peptide(protein_group_idxs):
            continue

        for protein_group_idx in protein_group_idxs:
            precursors_to_update: Dict[str, Tuple[float, int]] = {}
            for pq_idx, pq in enumerate(
                protein_group_results[protein_group_idx].precursorQuants
            ):
                # if quant.lfq_settings.combine_charge_state is set to true, charge is
                # always -1
                if pq.peptide == peptide and (pq.charge == charge or charge == -1):
                    if (
                        pq.post_err_prob
                        < precursors_to_update.get(pq.experiment, (1.01, np.nan))[0]
                    ):
                        precursors_to_update[pq.experiment] = (pq.post_err_prob, pq_idx)

            for filename, intensity in intensities:
                experiment, fraction = filename, -1
                if file_mapping:
                    experiment, fraction = file_mapping[filename]

                if intensity == 0.0:
                    continue

                if experiment in precursors_to_update:
                    pq_idx = precursors_to_update[experiment][1]
                    protein_group_results[protein_group_idx].precursorQuants[
                        pq_idx
                    ].intensity = intensity
                    if charge == -1:
                        protein_group_results[protein_group_idx].precursorQuants[
                            pq_idx
                        ].charge = charge
                else:
                    # match-between-runs hit
                    precursorQuant = PrecursorQuant(
                        peptide=peptide,
                        charge=charge,
                        experiment=experiment,
                        fraction=fraction,
                        intensity=intensity,
                        post_err_prob=np.nan,
                        tmt_intensities=None,  # TODO: add TMT support
                        silac_intensities=None,  # TODO: add SILAC support
                        evidence_id=-1,
                    )
                    protein_group_results[protein_group_idx].precursorQuants.append(
                        precursorQuant
                    )
    return protein_group_results


def add_precursor_quants_multiple(
    sage_results_files: List[str],
    sage_lfq_tsv: List[str],
    protein_groups: pg.ProteinGroups,
    protein_group_results: results.ProteinGroupResults,
    peptide_to_protein_maps: List[digest.PeptideToProteinMap],
    experimental_design: Optional[pd.DataFrame],
    discard_shared_peptides: bool,
    score_type: scoring_strategy.ProteinScoringStrategy,
    suppress_missing_peptide_warning: bool,
):
    post_err_probs_combined = []
    for sage_results_file in sage_results_files:
        protein_group_results, post_err_probs = add_precursor_quants(
            sage_results_file,
            protein_group_results,
            protein_groups,
            experimental_design,
            discard_shared_peptides,
            suppress_missing_peptide_warning,
        )
        post_err_probs_combined.extend(post_err_probs)

    protein_group_results = update_precursor_quants(
        protein_group_results,
        protein_groups,
        sage_lfq_tsv,
        experimental_design,
        discard_shared_peptides,
        suppress_missing_peptide_warning,
    )
    return protein_group_results, post_err_probs_combined
