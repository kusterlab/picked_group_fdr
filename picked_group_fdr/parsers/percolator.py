from __future__ import annotations

import os
import logging
from typing import Dict, Tuple, List

from .modifications import FIXED_MODS_DICTS, FIXED_MODS_UNIMOD
from . import modifications
from . import tsv

# for type hints only
from .. import scoring_strategy

logger = logging.getLogger(__name__)


# tuple of (search engine score, posterior error probability)
ScorePEPPair = Tuple[float, float]
# tuple of (scan number, modified sequence)
ScanPeptidePair = Tuple[int, str]
# dictionary of raw_file => dict(scan peptide pair => percolator result)
ResultsDict = Dict[str, Dict[ScanPeptidePair, ScorePEPPair]]


def is_native_percolator_file(headers):
    return "psmid" in headers


def is_mokapot_file(headers):
    return "specid" in headers


def is_ms2rescore_file(headers):
    return "spectrum_id" in headers


def get_percolator_column_idxs(headers: List[str]):
    headers = list(map(str.lower, headers))
    get_header_col = tsv.get_header_col_func(headers)
    if is_native_percolator_file(headers):
        id_col = get_header_col("psmid", required=True)
        filename_col = get_header_col("filename")
        pept_col = get_header_col("peptide", required=True)
        score_col = get_header_col("score", required=True)
        qval_col = get_header_col("q-value", required=True)
        post_err_prob_col = get_header_col("posterior_error_prob", required=True)
        protein_col = get_header_col("proteinids", required=True)
    elif is_mokapot_file(headers):
        id_col = get_header_col("specid", required=True)
        filename_col = get_header_col("filename")
        pept_col = get_header_col("peptide", required=True)
        score_col = get_header_col("mokapot score", required=True)
        qval_col = get_header_col("mokapot q-value", required=True)
        post_err_prob_col = get_header_col("mokapot pep", required=True)
        protein_col = get_header_col("proteins", required=True)
    elif is_ms2rescore_file(headers):
        id_col = get_header_col("spectrum_id", required=True)
        filename_col = get_header_col("run")
        if "peptidoform" in headers: # peptide-level output
            pept_col = get_header_col("peptidoform")
            score_col = get_header_col("score", required=True)
            qval_col = get_header_col("qvalue", required=True)
            post_err_prob_col = get_header_col("pep", required=True)
        else:
            pept_col = get_header_col("peptide", required=True)
            score_col = get_header_col("mokapot score", required=True)
            qval_col = get_header_col("mokapot q-value", required=True)
            post_err_prob_col = get_header_col("mokapot pep", required=True)    
        protein_col = get_header_col("protein_list", required=True)
    else:
        raise ValueError(
            "Could not determine percolator input file format. The file should either "
            "contain a column named PSMId (native percolator), SpecId (mokapot) or spectrum_id (ms2rescore)."
        )
    return (
        id_col,
        filename_col,
        pept_col,
        score_col,
        qval_col,
        post_err_prob_col,
        protein_col,
    )


def parse_percolator_out_file(
    reader,
    headers,
    get_proteins,
    score_type: scoring_strategy.ProteinScoringStrategy,
    **kwargs,
):
    (
        _,
        _,
        pept_col,
        score_col,
        _,
        post_err_prob_col,
        protein_col,
    ) = get_percolator_column_idxs(headers)
    headers = list(map(str.lower, headers))

    if score_type.get_score_column() == "posterior_error_prob":
        score_col = post_err_prob_col

    logger.info("Parsing Percolator output file")
    for line_idx, row in enumerate(reader):
        if line_idx % 500000 == 0:
            logger.info(f"    Reading line {line_idx}")

        modified_peptide = row[pept_col]
        if line_idx == 0:
            peptide_contains_flanks = modified_peptide.startswith(
                "-."
            ) and modified_peptide.endswith(".-")

        if peptide_contains_flanks:
            modified_peptide = modified_peptide[2:-2]

        experiment = 1
        score = float(row[score_col])

        if is_native_percolator_file(headers):
            proteins = row[protein_col:]
        elif is_mokapot_file(headers):
            proteins = row[protein_col].split("\t")
        elif is_ms2rescore_file(headers):
            proteins = eval(row[protein_col])
        else:
            raise ValueError("Could not deduce percolator output format from headers.")

        proteins = get_proteins(modified_peptide, proteins)
        if proteins:
            yield modified_peptide, proteins, experiment, score


def parse_percolator_out_file_to_dict(
    perc_out_file: str, results_dict: dict, input_type: str = ""
) -> Tuple[Dict[str, str], ResultsDict]:
    if not os.path.isfile(perc_out_file):
        raise FileNotFoundError(
            f"Could not find percolator output file {perc_out_file}. Please check if this file exists."
        )

    delimiter = tsv.get_delimiter(perc_out_file)
    reader = tsv.get_tsv_reader(perc_out_file, delimiter)
    headers = next(reader)  # save the header

    id_col, filename_col, pept_col, score_col, _, post_err_prob_col, _ = (
        get_percolator_column_idxs(headers)
    )

    logger.info("Parsing Percolator output file")
    convert_to_proforma = modifications.prosit_mod_to_proforma()
    fixed_mod_idx = -1
    first = True
    for line_idx, row in enumerate(reader):
        if line_idx % 500000 == 0:
            logger.info(f"    Reading line {line_idx}")

        modified_sequence = row[pept_col][2:-2]
        score = float(row[score_col])
        post_err_prob = float(row[post_err_prob_col])
        filename = ""
        if filename_col >= 0:
            filename = row[filename_col]

        psm_id = row[id_col]
        if input_type == "prosit":
            raw_file, scan_number, modified_sequence = parse_prosit_psmid_and_peptide(
                psm_id, modified_sequence, filename, convert_to_proforma
            )

            if first:
                for i, fixed_mod in enumerate(FIXED_MODS_UNIMOD):
                    if fixed_mod in modified_sequence:
                        fixed_mod_idx = i
                first = False
            elif fixed_mod_idx >= 0:
                if FIXED_MODS_UNIMOD[fixed_mod_idx] not in modified_sequence:
                    fixed_mod_idx = -1
        else:
            raw_file, scan_number, modified_sequence = (
                parse_andromeda_psmid_and_peptide(psm_id, modified_sequence)
            )

        results_dict[raw_file][(scan_number, modified_sequence)] = (
            score,
            post_err_prob,
        )

    return FIXED_MODS_DICTS[fixed_mod_idx + 1], results_dict


def parse_prosit_psmid_and_peptide(
    psm_id: str, modified_sequence: str, filename: str, convert_to_proforma
) -> tuple[str, int, str]:
    """Parse the Prosit PSMId and peptide.

    Args:
        psm_id (str): The Prosit PSMId is dash-separated string containing
            (raw_file, scan_number, modified_sequence, charge and optionally
            scan_event_number). The raw_file and modified_sequence can contain
            dashes themselves.
        peptide (str): The peptide sequence.
        convert_to_proforma (function): A function to convert the peptide to
            ProForma notation.

    Returns:
        tuple[str, int, str]: A tuple containing the parsed raw file, scan
            number, and peptide in ProForma notation.
    """
    # The original Prosit output files did not have a filename column, but
    # always included a scan event number, i.e. num_fields = 5.
    # The new Oktoberfest output has a filename column, but the scan event
    # number is optional, so instead we count the number of fields using the
    # number of dashes in the filename column.
    num_fields = 5
    if len(filename) > 0:
        num_fields = (
            psm_id.count("-") - filename.count("-") - modified_sequence.count("-") + 1
        )

    scan_number_idx = -1 * (num_fields - 1) - modified_sequence.count("-")
    scan_number = int(float(psm_id.split("-")[scan_number_idx]))

    if len(filename) == 0:
        filename = "-".join(psm_id.split("-")[:scan_number_idx])

    modified_sequence = convert_to_proforma(modified_sequence)
    return filename, scan_number, modified_sequence


def parse_andromeda_psmid_and_peptide(
    psm_id: str, modified_sequence: str
) -> tuple[str, int, str]:
    modified_sequence = modified_sequence.replace("[42]", "(ac)").replace(
        "M[16]", "M(ox)"
    )
    filename = "_".join(psm_id.split("_")[:-3])
    scan_number = int(psm_id.split("_")[-3])
    return filename, scan_number, modified_sequence


PERCOLATOR_NATIVE_HEADERS = [
    "PSMId",
    "score",
    "q-value",
    "posterior_error_prob",
    "peptide",
    "proteinIds",
]
