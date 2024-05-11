from __future__ import annotations

from typing import List

import numpy as np

from . import psm
from .. import helpers

# for type hints only
from .. import digest
from .. import scoring_strategy
from .. import peptide_info


def parse_evidence_files(
    evidence_files: List[str],
    peptide_to_protein_maps: List[digest.PeptideToProteinMap],
    score_type: scoring_strategy.ProteinScoringStrategy,
    suppress_missing_peptide_warning: bool,
) -> peptide_info.PeptideInfoList:
    """Returns best score per peptide"""
    peptide_info_list = dict()
    for modified_peptide, proteins, _, score in psm.parse_evidence_file_multiple(
        evidence_files,
        peptide_to_protein_maps=peptide_to_protein_maps,
        score_type=score_type,
        suppress_missing_peptide_warning=suppress_missing_peptide_warning,
    ):
        peptide = helpers.remove_modifications(modified_peptide)
        if np.isnan(score) or score >= peptide_info_list.get(peptide, [np.inf])[0]:
            continue

        peptide_info_list[peptide] = [score, proteins]

    return peptide_info_list
