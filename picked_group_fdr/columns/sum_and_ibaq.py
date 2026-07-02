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
        self, num_ibaq_peptides_per_protein: Dict[str, int],
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

        proteinGroupCounts = np.zeros(
            len(experiment_to_idx_map) * (1 + num_silac_channels), dtype="int"
        )
        for pgr in protein_group_results:
            intensities = _get_intensities(
                pgr.precursorQuants,
                experiment_to_idx_map,
                post_err_prob_cutoff,
                num_silac_channels,
            )

            totalIntensity = sum(intensities[:: num_silac_channels + 1])
            pgr.append(totalIntensity)
            pgr.extend(intensities)

            if pgr.qValue < self.protein_group_fdr_threshold:
                proteinGroupCounts += np.array(
                    [1 if intensity > 0 else 0 for intensity in intensities]
                )

            # iBAQ: divide intensity by number of (fully tryptic) theoretical peptides
            numTheoreticalPeptides = [
                self.num_ibaq_peptides_per_protein[p] for p in pgr.proteinIds.split(";")
            ]
            leadingProteinNumPeptides = max([1, numTheoreticalPeptides[0]])
            pgr.append(";".join(map(str, numTheoreticalPeptides)))
            pgr.append(totalIntensity / leadingProteinNumPeptides)
            pgr.extend([i / leadingProteinNumPeptides for i in intensities])

        logger.info(
            f"#Protein groups quantified ({self.protein_group_fdr_threshold*100:g}% protein group-level FDR, summed intensity / iBAQ):"
        )
        for experiment, numProteinGroups in zip(
            experiment_to_idx_map.keys(),
            helpers.chunks(proteinGroupCounts, num_silac_channels + 1),
        ):
            logger.info(f"    {experiment}: {numProteinGroups[0]}")
            for silacIdx, silacChannel in enumerate(silac_channels):
                logger.info(
                    f"    {experiment} {silacChannel}: {numProteinGroups[silacIdx+1]}"
                )


def _get_intensities(
    peptideIntensityList: List[precursor_quant.PrecursorQuant],
    experimentToIdxMap,
    postErrProbCutoff,
    numSilacChannels,
) -> List[float]:
    intensities = [0.0] * (len(experimentToIdxMap) * (1 + numSilacChannels))
    for precursor in peptideIntensityList:
        if np.isnan(precursor.intensity):
            continue

        if (
            helpers.is_mbr(precursor.post_err_prob)
            or precursor.post_err_prob <= postErrProbCutoff
        ):
            intensities[
                experimentToIdxMap[precursor.experiment] * (1 + numSilacChannels)
            ] += precursor.intensity
            if precursor.silac_intensities is not None:
                for silacIdx, silacIntensity in enumerate(precursor.silac_intensities):
                    intensities[
                        experimentToIdxMap[precursor.experiment]
                        * (1 + numSilacChannels)
                        + silacIdx
                        + 1
                    ] += silacIntensity
    return intensities
