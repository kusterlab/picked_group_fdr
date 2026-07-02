from __future__ import annotations

from typing import List, Dict
import logging

import numpy as np

from .. import helpers
from .base import ProteinGroupColumns

# for type hints only
from .. import precursor_quant
from .. import results

logger = logging.getLogger(__name__)


def get_silac_channels(num_silac_channels: int):
    if num_silac_channels == 3:
        return ["L", "M", "H"]
    elif num_silac_channels == 2:
        return ["L", "H"]
    elif num_silac_channels > 0:
        raise ValueError("ERROR: Found a number of SILAC channels not equal to 2 or 3")
    return list()


class SummedIntensityAndIbaqColumns(ProteinGroupColumns):
    """Summed intensity & iBAQ"""

    num_ibaq_peptides_per_protein: Dict[str, int]

    def __init__(
        self,
        num_ibaq_peptides_per_protein: Dict[str, int],
        protein_group_fdr_threshold: float = 0.01,
    ):
        self.num_ibaq_peptides_per_protein = num_ibaq_peptides_per_protein

        # only used for reporting number of protein groups at the given threshold
        self.protein_group_fdr_threshold = protein_group_fdr_threshold

    def append_headers(
        self,
        protein_group_results: results.ProteinGroupResults,
    ) -> None:
        silac_channels = get_silac_channels(protein_group_results.num_silac_channels)
        protein_group_results.append_header("Intensity")
        for experiment in protein_group_results.experiments:
            protein_group_results.append_header("Intensity " + experiment)
            for silacChannel in silac_channels:
                protein_group_results.append_header(
                    "Intensity " + silacChannel + " " + experiment
                )

        protein_group_results.append_header("Number of theoretical peptides iBAQ")
        protein_group_results.append_header("iBAQ")
        for experiment in protein_group_results.experiments:
            protein_group_results.append_header("iBAQ " + experiment)
            for silacChannel in silac_channels:
                protein_group_results.append_header(
                    "iBAQ " + silacChannel + " " + experiment
                )

    def append_columns(
        self,
        protein_group_results: results.ProteinGroupResults,
        post_err_prob_cutoff: float,
    ) -> None:
        logger.info("Doing quantification: summed peptide intensity")

        experiment_to_idx_map = protein_group_results.get_experiment_to_idx_map()
        silac_channels = get_silac_channels(protein_group_results.num_silac_channels)
        num_silac_channels = len(silac_channels)

        protein_group_counts = np.zeros(
            len(experiment_to_idx_map) * (1 + num_silac_channels), dtype="int"
        )
        for pgr in protein_group_results:
            intensities = _get_intensities(
                pgr.precursorQuants,
                experiment_to_idx_map,
                post_err_prob_cutoff,
                num_silac_channels,
            )

            total_intensity = sum(intensities[:: num_silac_channels + 1])
            pgr.append(total_intensity)
            pgr.extend(intensities)

            if pgr.qValue < self.protein_group_fdr_threshold:
                protein_group_counts += np.array(
                    [1 if intensity > 0 else 0 for intensity in intensities]
                )

            # iBAQ: divide intensity by number of (fully tryptic) theoretical peptides
            num_theoretical_peptides = [
                self.num_ibaq_peptides_per_protein[p] for p in pgr.proteinIds.split(";")
            ]
            leading_protein_num_peptides = max([1, num_theoretical_peptides[0]])
            pgr.append(";".join(map(str, num_theoretical_peptides)))
            pgr.append(total_intensity / leading_protein_num_peptides)
            pgr.extend([i / leading_protein_num_peptides for i in intensities])

        logger.info(
            f"#Protein groups quantified ({self.protein_group_fdr_threshold*100:g}% protein group-level FDR, summed intensity / iBAQ):"
        )
        for experiment, num_protein_groups in zip(
            experiment_to_idx_map.keys(),
            helpers.chunks(protein_group_counts, num_silac_channels + 1),
        ):
            logger.info(f"    {experiment}: {num_protein_groups[0]}")
            for silac_idx, silac_channel in enumerate(silac_channels):
                logger.info(
                    f"    {experiment} {silac_channel}: {num_protein_groups[silac_idx+1]}"
                )


def _get_intensities(
    peptide_intensity_list: List[precursor_quant.PrecursorQuant],
    experiment_to_idx_map,
    post_err_prob_cutoff,
    num_silac_channels,
    remove_silac_summed_intensity_columns: bool = False,
) -> List[float]:
    intensities = [0.0] * (len(experiment_to_idx_map) * (1 + num_silac_channels))
    for precursor in peptide_intensity_list:
        if np.isnan(precursor.intensity):
            continue

        if (
            helpers.is_mbr(precursor.post_err_prob)
            or precursor.post_err_prob <= post_err_prob_cutoff
        ):
            intensities[
                experiment_to_idx_map[precursor.experiment] * (1 + num_silac_channels)
            ] += precursor.intensity
            if precursor.silac_intensities is not None:
                for silac_idx, silac_intensity in enumerate(
                    precursor.silac_intensities
                ):
                    intensities[
                        experiment_to_idx_map[precursor.experiment]
                        * (1 + num_silac_channels)
                        + silac_idx
                        + 1
                    ] += silac_intensity
    if remove_silac_summed_intensity_columns and num_silac_channels > 0:
        # removes the column intensities per experiment with the SILAC
        # channels summed up
        del intensities[:: (num_silac_channels + 1)]
    return intensities
