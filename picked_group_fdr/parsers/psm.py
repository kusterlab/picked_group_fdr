from __future__ import annotations

from typing import Callable, Dict, List, Optional
import logging

from .. import digest
from .. import helpers
from .. import scoring_strategy

from . import tsv


logger = logging.getLogger(__name__)


def get_peptide_to_protein_mapper(
    peptide_to_protein_map: digest.PeptideToProteinMap,
    score_type: scoring_strategy.ProteinScoringStrategy,
    suppress_missing_peptide_warning: bool,
) -> Callable[[str, List[str]], List[str]]:
    """Returns a function that updates the list of proteins for a peptide.

    In most cases, the returned function simply returns the same list of proteins as 
    were used as input. In some cases, the mapping from peptide to protein is not
    accurate (e.g. MaxQuant) and a remapping is applied. Note that this is also the 
    place where razor peptide decisions are made.

    Args:
        peptide_to_protein_map (Dict): _description_
        score_type (scoring.ProteinScoringStrategy): _description_
        suppress_missing_peptide_warning (bool): _description_

    Returns:
        Callable[[str, List[str]], List[str]]: a function that maps a peptide and a list
            of proteins to an updated list of proteins.
    """    
    def get_proteins(peptide, tmp_proteins):
        if score_type.remaps_peptides_to_proteins():
            proteins = digest.get_proteins(peptide_to_protein_map, peptide)
            if len(proteins) == 0:
                if (
                    not helpers.is_contaminant(tmp_proteins)
                    and not suppress_missing_peptide_warning
                ):
                    logger.warning(f"Missing peptide: {peptide} {tmp_proteins}")
                return None
        else:
            proteins = tmp_proteins

        # filtering for razor peptide approach
        proteins = score_type.filter_proteins(proteins)

        return helpers.remove_decoy_proteins_from_target_peptides(proteins)

    return get_proteins


def parse_evidence_file_single(
    evidence_file: str,
    peptide_to_protein_map: Dict,
    score_type: Optional[scoring_strategy.ProteinScoringStrategy],
    for_quantification: bool = False,
    suppress_missing_peptide_warning: bool = False,
):
    if score_type is None:
        score_type = scoring_strategy.ProteinScoringStrategy("bestPEP")

    delimiter = tsv.get_delimiter(evidence_file)
    reader = tsv.get_tsv_reader(evidence_file, delimiter)
    headers = next(reader)

    get_proteins = get_peptide_to_protein_mapper(
        peptide_to_protein_map, score_type, suppress_missing_peptide_warning
    )

    yield from score_type.get_evidence_parser()(
        reader, headers, get_proteins, score_type, for_quantification=for_quantification
    )


def parse_evidence_file_multiple(
    evidence_files: List[str],
    peptide_to_protein_maps: List[Dict],
    score_type: scoring_strategy.ProteinScoringStrategy,
    for_quantification: bool = False,
    suppress_missing_peptide_warning: bool = False,
):
    for evidence_file, peptide_to_protein_map in zip(
        evidence_files, peptide_to_protein_maps
    ):
        yield from parse_evidence_file_single(
            evidence_file,
            peptide_to_protein_map,
            score_type,
            for_quantification,
            suppress_missing_peptide_warning,
        )
