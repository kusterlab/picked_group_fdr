from __future__ import annotations

from typing import Dict, List
import logging

from .. import digest
from .. import helpers
from .. import scoring

from . import tsv
from . import fragpipe
from . import maxquant
from . import percolator


logger = logging.getLogger(__name__)


def get_peptide_to_protein_mapper(
    peptide_to_protein_map: Dict,
    score_type: scoring.ProteinScoringStrategy,
    suppress_missing_peptide_warning: bool,
):
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
    score_type: scoring.ProteinScoringStrategy,
    for_quantification: bool = False,
    suppress_missing_peptide_warning: bool = False,
):
    delimiter = tsv.get_delimiter(evidence_file)
    reader = tsv.get_tsv_reader(evidence_file, delimiter)
    headers = next(reader)

    get_proteins = get_peptide_to_protein_mapper(
        peptide_to_protein_map, score_type, suppress_missing_peptide_warning
    )

    if percolator.is_percolator_file(headers):
        yield from percolator.parse_percolator_out_file(
            reader, headers, get_proteins, score_type
        )
    elif fragpipe.is_fragpipe_psm_file(headers):
        yield from fragpipe.parse_fragpipe_psm_file(
            reader, headers, get_proteins, score_type
        )
    else:
        # convert headers to lowercase since MQ changes the capitalization frequently
        headers = list(map(str.lower, headers))
        if evidence_file.endswith(".csv"):
            headers = [x.replace(".", " ") for x in headers]
        yield from maxquant.parse_mq_evidence_file(
            reader, headers, get_proteins, score_type, for_quantification
        )


def parse_evidence_file_multiple(
    evidence_files: List[str],
    peptide_to_protein_maps: List[Dict],
    score_type: scoring.ProteinScoringStrategy,
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