from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd

from .. import helpers
from ..parsers import fragpipe, tsv
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
    experiment: str,
    discard_shared_peptides: bool,
    suppress_missing_peptide_warning: bool,
):
    protein_group_results.experiments.append(experiment)

    delimiter = tsv.get_delimiter(fragpipe_psm_file)
    reader = tsv.get_tsv_reader(fragpipe_psm_file, delimiter)
    headers = next(reader)

    post_err_probs = []
    for (
        peptide,
        charge,
        post_err_prob,
        assigned_mods,
        observed_mods,
        proteins,
    ) in fragpipe.parse_fragpipe_psm_file_for_protein_tsv(reader, headers):
        protein_group_idxs = protein_groups.get_protein_group_idxs(proteins)

        if helpers.is_missing_in_protein_groups(protein_group_idxs):
            if not suppress_missing_peptide_warning:
                logger.debug(
                    f"Could not find any of the proteins {proteins} in proteinGroups.txt"
                )
            continue

        if discard_shared_peptides and helpers.is_shared_peptide(protein_group_idxs):
            continue

        if not helpers.is_decoy(proteins):
            post_err_probs.append((post_err_prob, experiment, experiment, peptide))

        for protein_group_idx in protein_group_idxs:
            precursor_quant = PrecursorQuant(
                peptide=peptide,
                charge=charge,
                experiment=experiment,
                fraction=-1,
                intensity=np.nan,
                post_err_prob=post_err_prob,
                tmt_intensities=None,  # TODO: add TMT support
                silac_intensities=None,  # TODO: add SILAC support
                evidence_id=-1,
                assigned_mods=assigned_mods,
                observed_mods=observed_mods,
            )
            protein_group_results[protein_group_idx].precursorQuants.append(
                precursor_quant
            )
    return protein_group_results, post_err_probs


def update_precursor_quants(
    protein_group_results: results.ProteinGroupResults,
    protein_groups: pg.ProteinGroups,
    combined_ion_files: List[str],
    discard_shared_peptides: bool,
    suppress_missing_peptide_warning: bool,
):
    for combined_ion_file in combined_ion_files:
        protein_group_results = update_precursor_quants_single(
            protein_group_results,
            protein_groups,
            combined_ion_file,
            discard_shared_peptides,
            suppress_missing_peptide_warning,
        )
    return protein_group_results


def update_precursor_quants_single(
    protein_group_results: results.ProteinGroupResults,
    protein_groups: pg.ProteinGroups,
    combined_ion_file: str,
    discard_shared_peptides: bool,
    suppress_missing_peptide_warning: bool,
):
    delimiter = tsv.get_delimiter(combined_ion_file)
    reader = tsv.get_tsv_reader(combined_ion_file, delimiter)
    headers = next(reader)

    for (
        peptide,
        charge,
        assigned_mods,
        proteins,
        intensities,
    ) in fragpipe.parse_fragpipe_combined_ion_file(reader, headers):
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
                if (
                    pq.peptide == peptide
                    and pq.charge == charge
                    and pq.assigned_mods == assigned_mods
                ):
                    if (
                        pq.post_err_prob
                        < precursors_to_update.get(pq.experiment, (1.01, np.nan))[0]
                    ):
                        precursors_to_update[pq.experiment] = (pq.post_err_prob, pq_idx)

            for experiment, intensity in intensities:
                if intensity == 0.0:
                    continue

                if experiment in precursors_to_update:
                    pq_idx = precursors_to_update[experiment][1]
                    protein_group_results[protein_group_idx].precursorQuants[
                        pq_idx
                    ].intensity = intensity
                else:
                    # match-between-runs hit
                    precursor_quant = PrecursorQuant(
                        peptide=peptide,
                        charge=charge,
                        experiment=experiment,
                        fraction=-1,
                        intensity=intensity,
                        post_err_prob=np.nan,
                        tmt_intensities=None,  # TODO: add TMT support
                        silac_intensities=None,  # TODO: add SILAC support
                        evidence_id=-1,
                        assigned_mods=assigned_mods,
                        observed_mods="",
                    )
                    protein_group_results[protein_group_idx].precursorQuants.append(
                        precursor_quant
                    )
    return protein_group_results


def add_precursor_quants_multiple(
    fragpipe_psm_files: List[str],
    combined_ion_files: List[str],
    protein_groups: pg.ProteinGroups,
    protein_group_results: results.ProteinGroupResults,
    peptide_to_protein_maps: List[digest.PeptideToProteinMap],
    experimental_design: Optional[pd.DataFrame],
    discard_shared_peptides: bool,
    score_type: scoring_strategy.ProteinScoringStrategy,
    suppress_missing_peptide_warning: bool,
):
    post_err_probs_combined = []
    for fragpipe_psm_file in fragpipe_psm_files:
        experiment = Path(fragpipe_psm_file).parent.name
        protein_group_results, post_err_probs = add_precursor_quants(
            fragpipe_psm_file,
            protein_group_results,
            protein_groups,
            experiment,
            discard_shared_peptides,
            suppress_missing_peptide_warning,
        )
        post_err_probs_combined.extend(post_err_probs)

    if combined_ion_files is not None:
        protein_group_results = update_precursor_quants(
            protein_group_results,
            protein_groups,
            combined_ion_files,
            discard_shared_peptides,
            suppress_missing_peptide_warning,
        )
    return protein_group_results, post_err_probs_combined
