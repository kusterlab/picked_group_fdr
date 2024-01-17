import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd

from .. import helpers
from ..parsers import fragpipe, tsv
from ..precursor_quant import PrecursorQuant
from ..protein_groups import ProteinGroups
from ..results import ProteinGroupResults

# for type hints only
from .. import digest

logger = logging.getLogger(__name__)


def add_precursor_quants(
    fragpipe_psm_file: str,
    protein_group_results: ProteinGroupResults,
    protein_groups: ProteinGroups,
    experiment: str,
    discard_shared_peptides: bool,
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

        if len(protein_group_idxs) == 0:
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
    protein_group_results: ProteinGroupResults,
    protein_groups: ProteinGroups,
    combined_ion_file: str,
    discard_shared_peptides: bool,
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

        if len(protein_group_idxs) == 0:
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
    combined_ion_file: List[str],
    protein_groups: ProteinGroups,
    protein_group_results: ProteinGroupResults,
    peptide_to_protein_maps: List[digest.PeptideToProteinMap],
    experimental_design: Optional[pd.DataFrame],
    discard_shared_peptides: bool,
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
        )
        post_err_probs_combined.extend(post_err_probs)

    if combined_ion_file is not None:
        protein_group_results = update_precursor_quants(
            protein_group_results,
            protein_groups,
            combined_ion_file,
            discard_shared_peptides,
        )
    return protein_group_results, post_err_probs_combined
