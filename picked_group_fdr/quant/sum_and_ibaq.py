from typing import List, Dict
import logging

import numpy as np

from .. import helpers
from .base import ProteinGroupColumns

# for type hints only
from .precursor_quant import PrecursorQuant
from ..results import ProteinGroupResults

logger = logging.getLogger(__name__)


class SummedIntensityAndIbaqColumns(ProteinGroupColumns):
    """Summed intensity & iBAQ"""

    silacChannels: List[str]
    numIbaqPeptidesPerProtein: Dict[str, int]

    def __init__(
        self, silacChannels: List[str], numIbaqPeptidesPerProtein: Dict[str, int]
    ):
        self.silacChannels = silacChannels
        self.numIbaqPeptidesPerProtein = numIbaqPeptidesPerProtein

    def append_headers(
        self,
        proteinGroupResults: ProteinGroupResults,
        experiments: List[str],
    ) -> None:
        proteinGroupResults.append_header("Intensity")
        for experiment in experiments:
            proteinGroupResults.append_header("Intensity " + experiment)
            for silacChannel in self.silacChannels:
                proteinGroupResults.append_header(
                    "Intensity " + silacChannel + " " + experiment
                )

        proteinGroupResults.append_header("Number of theoretical peptides iBAQ")
        proteinGroupResults.append_header("iBAQ")
        for experiment in experiments:
            proteinGroupResults.append_header("iBAQ " + experiment)
            for silacChannel in self.silacChannels:
                proteinGroupResults.append_header(
                    "iBAQ " + silacChannel + " " + experiment
                )

    def append_columns(
        self,
        proteinGroupResults: ProteinGroupResults,
        experimentToIdxMap: Dict[str, int],
        postErrProbCutoff: float,
    ) -> None:
        logger.info("Doing quantification: summed peptide intensity")
        numSilacChannels = len(self.silacChannels)

        proteinGroupCounts = np.zeros(
            len(experimentToIdxMap) * (1 + numSilacChannels), dtype="int"
        )
        for pgr in proteinGroupResults:
            intensities = _get_intensities(
                pgr.precursorQuants,
                experimentToIdxMap,
                postErrProbCutoff,
                numSilacChannels,
            )

            totalIntensity = sum(intensities[:: numSilacChannels + 1])
            pgr.append(totalIntensity)
            pgr.extend(intensities)

            if pgr.qValue < 0.01:
                proteinGroupCounts += np.array(
                    [1 if intensity > 0 else 0 for intensity in intensities]
                )

            # iBAQ: divide intensity by number of (fully tryptic) theoretical peptides
            numTheoreticalPeptides = [
                self.numIbaqPeptidesPerProtein[p] for p in pgr.proteinIds.split(";")
            ]
            leadingProteinNumPeptides = max([1, numTheoreticalPeptides[0]])
            pgr.append(";".join(map(str, numTheoreticalPeptides)))
            pgr.append(totalIntensity / leadingProteinNumPeptides)
            pgr.extend([i / leadingProteinNumPeptides for i in intensities])

        logger.info(
            "#Protein groups quantified (1% protein group-level FDR, summed intensity / iBAQ):"
        )
        for experiment, numProteinGroups in zip(
            experimentToIdxMap.keys(),
            helpers.chunks(proteinGroupCounts, numSilacChannels + 1),
        ):
            logger.info(f"    {experiment}: {numProteinGroups[0]}")
            for silacIdx, silacChannel in enumerate(self.silacChannels):
                logger.info(
                    f"    {experiment} {silacChannel}: {numProteinGroups[silacIdx+1]}"
                )


def _get_intensities(
    peptideIntensityList: List[PrecursorQuant],
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
