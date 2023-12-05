import logging
from typing import Dict, List, Tuple

import numpy as np

from .. import helpers
from ..parsers import sage, tsv
from ..precursor_quant import PrecursorQuant
from ..protein_groups import ProteinGroups
from ..results import ProteinGroupResults
from ..scoring import ProteinScoringStrategy

logger = logging.getLogger(__name__)


def add_precursor_quants(
    fragpipe_psm_file: str,
    protein_group_results: ProteinGroupResults,
    protein_groups: ProteinGroups,
    discard_shared_peptides: bool,
):
    delimiter = tsv.get_delimiter(fragpipe_psm_file)
    reader = tsv.get_tsv_reader(fragpipe_psm_file, delimiter)
    headers = next(reader)

    get_proteins = lambda peptide, proteins: proteins
    score_type = ProteinScoringStrategy("Sage bestPEP")

    post_err_probs = []
    for (
        peptide,
        proteins,
        experiment,
        post_err_prob,
        charge,
    ) in sage.parse_sage_results_file(
        reader, headers, get_proteins, score_type, for_quantification=True
    ):
        protein_group_idxs = protein_groups.get_protein_group_idxs(proteins)

        if len(protein_group_idxs) == 0:
            logger.debug(
                f"Could not find any of the proteins {proteins} in proteinGroups.txt"
            )
            continue

        if discard_shared_peptides and helpers.is_shared_peptide(protein_group_idxs):
            continue

        if not helpers.is_decoy(proteins):
            post_err_probs.append((post_err_prob, "", experiment, peptide))

        for protein_group_idx in protein_group_idxs:
            precursorQuant = PrecursorQuant(
                peptide=peptide,
                charge=charge,
                experiment=experiment,
                fraction=-1,
                intensity=np.nan,
                post_err_prob=post_err_prob,
                tmt_intensities=None,  # TODO: add TMT support
                silac_intensities=None,  # TODO: add SILAC support
                evidence_id=-1,
            )
            protein_group_results[protein_group_idx].precursorQuants.append(
                precursorQuant
            )
    return protein_group_results, post_err_probs


def update_precursor_quants(
    protein_group_results: ProteinGroupResults,
    protein_groups: ProteinGroups,
    sage_lfq_tsv: str,
    discard_shared_peptides: bool,
):
    delimiter = tsv.get_delimiter(sage_lfq_tsv)
    reader = tsv.get_tsv_reader(sage_lfq_tsv, delimiter)
    headers = next(reader)

    protein_group_results.experiments = sage.get_experiments_from_sage_lfq_headers(
        headers
    )

    for (
        peptide,
        charge,
        proteins,
        intensities,
    ) in sage.parse_sage_lfq_file(reader, headers):
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
                # if quant.lfq_settings.combine_charge_state is set to true, charge is
                # always -1
                if pq.peptide == peptide and (pq.charge == charge or charge == -1):
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
                        fraction=-1,
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
    combined_ion_file: str,
    protein_groups: ProteinGroups,
    protein_group_results: ProteinGroupResults,
    discard_shared_peptides: bool = True,
):
    post_err_probs_combined = []
    for sage_results_file in sage_results_files:
        protein_group_results, post_err_probs = add_precursor_quants(
            sage_results_file,
            protein_group_results,
            protein_groups,
            discard_shared_peptides,
        )
        post_err_probs_combined.extend(post_err_probs)

    protein_group_results = update_precursor_quants(
        protein_group_results,
        protein_groups,
        combined_ion_file,
        discard_shared_peptides,
    )
    return protein_group_results, post_err_probs_combined