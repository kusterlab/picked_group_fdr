import os
import csv
import re
from typing import List, Dict, Optional
import logging

import numpy as np

from picked_group_fdr import digest, helpers

logger = logging.getLogger(__name__)


# csv.field_size_limit(sys.maxsize)
csv.field_size_limit(2147483647)


TMT_UNIMOD = "[UNIMOD:737]"  # TMT6/10/11
TMTPRO_UNIMOD = "[UNIMOD:2016]"
ITRAQ4_UNIMOD = "[UNIMOD:214]"
ITRAQ8_UNIMOD = "[UNIMOD:730]"


# copied from fundamentals/constants.py
MAXQUANT_VAR_MODS = {
    "(ox)": "[UNIMOD:35]",
    "(Oxidation (M))": "[UNIMOD:35]",
    "(tm)": "[UNIMOD:737]",
    "_(tm)": f"_{TMT_UNIMOD}",
    "K(tm)": f"K{TMT_UNIMOD}",
    "_(TMTPro (N-term))": f"_{TMTPRO_UNIMOD}",
    "K(TMTPro (K))": f"K{TMTPRO_UNIMOD}",
    "_(iTRAQ4plex (N-term))": f"_{ITRAQ4_UNIMOD}",
    "K(iTRAQ4plex (K))": f"K{ITRAQ4_UNIMOD}",
    "_(iTRAQ8plex (N-term))": f"_{ITRAQ8_UNIMOD}",
    "K(iTRAQ8plex (K))": f"K{ITRAQ8_UNIMOD}",
    "(ph)": "[UNIMOD:21]",
    "K(Lys8)": "K[UNIMOD:259]",
    "R(Arg10)": "R[UNIMOD:267]",
    "C(Carbamidomethyl (C))": "C[UNIMOD:4]",
}

DEFAULT_FIXED_MODS = {"C": "C[UNIMOD:4]"}

SILAC_HEAVY_FIXED_MODS = {
    "C": "C[UNIMOD:4]",
    "K": "K[UNIMOD:259]",
    "R": "R[UNIMOD:267]",
}

TMT_FIXED_MODS = {"C": "C[UNIMOD:4]", "^_": f"_{TMT_UNIMOD}", "K": f"K{TMT_UNIMOD}"}

TMTPRO_FIXED_MODS = {
    "C": "C[UNIMOD:4]",
    "^_": f"_{TMTPRO_UNIMOD}",
    "K": f"K{TMTPRO_UNIMOD}",
}

ITRAQ4_FIXED_MODS = {
    "C": "C[UNIMOD:4]",
    "^_": f"_{ITRAQ4_UNIMOD}",
    "K": f"K{ITRAQ4_UNIMOD}",
}

ITRAQ8_FIXED_MODS = {
    "C": "C[UNIMOD:4]",
    "^_": f"_{ITRAQ8_UNIMOD}",
    "K": f"K{ITRAQ8_UNIMOD}",
}

FIXED_MODS_UNIMOD = [TMT_UNIMOD, TMTPRO_UNIMOD, ITRAQ4_UNIMOD, ITRAQ8_UNIMOD]
FIXED_MODS_DICTS = [
    DEFAULT_FIXED_MODS,
    TMT_FIXED_MODS,
    TMTPRO_FIXED_MODS,
    ITRAQ4_FIXED_MODS,
    ITRAQ8_FIXED_MODS,
]

PERCOLATOR_NATIVE_HEADERS = [
    "PSMId",
    "score",
    "q-value",
    "posterior_error_prob",
    "peptide",
    "proteinIds",
]


def parse_protein_groups_file_multiple(
    protein_groups_files: List[str], are_decoy_file: List[bool], **kwargs
):
    for protein_groups_file, is_decoy_file in zip(protein_groups_files, are_decoy_file):
        yield from parse_protein_groups_file_single(
            protein_groups_file, is_decoy_file=is_decoy_file, **kwargs
        )


def parse_protein_groups_file_single(
    protein_groups_file: str,
    protein_column: str = "Protein IDs",
    score_column: str = "Score",
    is_decoy_file: bool = False,
):
    """Parse protein groups file from MaxQuant or ProteomeDiscoverer.

    For PD with Chimerys, the column names are:
    - protein_column="Accession"
    - score_column="Score CHIMERY CHIMERYS"

    Args:
        protein_groups_file (_type_): _description_
        protein_column (str, optional): _description_. Defaults to 'Protein IDs'.
        score_column (str, optional): _description_. Defaults to 'Score'.
        is_decoy_file (bool, optional): _description_. Defaults to False.

    Yields:
        _type_: _description_
    """
    delimiter = get_delimiter(protein_groups_file)

    reader = get_tsv_reader(protein_groups_file, delimiter)
    headers = next(reader)  # save the header

    score_col = get_column_index(headers, score_column)
    protein_col = get_column_index(headers, protein_column)

    logger.info(f"Parsing proteinGroups file: {protein_groups_file}")
    for row in reader:
        proteins = list(map(str.strip, row[protein_col].split(";")))
        if is_decoy_file:
            proteins = [f"REV__{p}" for p in proteins]
        score = -100.0
        if len(row[score_col]) > 0:
            score = float(row[score_col])
        yield proteins, score


def parse_peptides_files_multiple(
    peptides_files: List[str], are_decoy_file: List[bool], **kwargs
):
    for peptide_file, is_decoy_file in zip(peptides_files, are_decoy_file):
        yield from parse_peptides_file_single(
            peptide_file, is_decoy_file=is_decoy_file, **kwargs
        )


def parse_peptides_file_single(
    peptides_file: str,
    peptide_column: str = "Modified sequence",
    protein_column: str = "Protein IDs",
    score_column: str = "Score",
    is_decoy_file: bool = False,
):
    """Parse peptide-level file from MaxQuant (evidence/msms.txt) or ProteomeDiscoverer.

    For PD with Chimerys, the column names are:
    - protein_column="Accession"
    - score_column="Score CHIMERY CHIMERYS"

    Args:
        peptides_file (_type_): evidence.txt or msms.txt
        peptide_column (str, optional): column name for peptide sequence. Defaults to 'Protein IDs'.
        protein_column (str, optional): column name for protein identifiers. Defaults to 'Protein IDs'.
        score_column (str, optional): column name for score. Defaults to 'Score'.
        is_decoy_file (bool, optional): if this file only contains decoy peptides. Defaults to False.

    Yields:
        (str, str, str, float): Peptide, Proteins, Experiment, Score
    """
    delimiter = get_delimiter(peptides_file)

    reader = get_tsv_reader(peptides_file, delimiter)
    headers = next(reader)  # save the header

    peptide_col = get_column_index(headers, peptide_column)
    score_col = get_column_index(headers, score_column)
    if not is_decoy_file:
        protein_col = get_column_index(headers, protein_column)

    experiment = "Experiment1"

    logger.info(f"Parsing peptides file: {peptides_file}")
    for row in reader:
        proteins = []
        if is_decoy_file:
            proteins = ["REV__protein"]
        elif len(row[protein_col]) > 0:
            proteins = list(map(str.strip, row[protein_col].split(";")))

        score = -100.0
        if len(row[score_col]) > 0:
            score = float(row[score_col])
            if np.isnan(score):
                continue
        yield row[peptide_col], proteins, experiment, score


def parse_evidence_file_multiple(
    evidence_files,
    peptide_to_protein_maps,
    score_type,
    for_quantification=False,
    suppress_missing_peptide_warning=False,
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


def parse_evidence_file_single(
    evidence_file,
    peptide_to_protein_map,
    score_type,
    for_quantification=False,
    suppress_missing_peptide_warning=False,
):
    delimiter = get_delimiter(evidence_file)
    reader = get_tsv_reader(evidence_file, delimiter)
    headers = next(reader)  # save the header

    get_proteins = get_peptide_to_protein_mapper(
        peptide_to_protein_map, score_type, suppress_missing_peptide_warning
    )

    if is_percolator_file(headers):
        yield from parse_percolator_out_file(reader, headers, get_proteins, score_type)
    else:
        headers = list(
            map(str.lower, headers)
        )  # convert headers to lowercase since MQ changes the capitalization frequently
        if evidence_file.endswith(".csv"):
            headers = [x.replace(".", " ") for x in headers]
        yield from parse_mq_evidence_file(
            reader, headers, get_proteins, score_type, for_quantification
        )


def get_peptide_to_protein_mapper(
    peptide_to_protein_map, score_type, suppress_missing_peptide_warning
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


def parse_mq_evidence_file(
    reader, headers, get_proteins, score_type, for_quantification=False
):
    """
    Reads in approximately 100,000 lines per second with for_quantification=False
    and 50,000 lines per second with for_quantification=True

    Columns needed for identification:
    - Modified sequence
    - Leading proteins (not used for all methods, but should still be in the input file)
    - Score/PEP
    - Experiment

    Extra columns needed for quantification:
    - Charge
    - Intensity
    - Raw file
    - Fraction (optional)
    - Id
    - Reporter intensity corrected (for TMT)
    - Reporter intensity (for TMT)
    - Reporter intensity count (for TMT)
    - Intensity L (for SILAC)
    - Intensity H (for SILAC)
    - Intensity M (optional, for SILAC)
    """
    get_header_col = get_header_col_func(headers)
    get_header_cols_starting_with = get_header_cols_starting_with_func(headers)

    pept_col = get_header_col("modified sequence", required=True)

    protein_col = get_header_col(
        "leading proteins", required=True
    )  # all protein groups, each represented by the first protein in the group
    if score_type.use_razor:
        protein_col = get_header_col(
            "leading razor protein", required=True
        )  # best scoring protein group, represented by the first protein in the group

    score_col = get_header_col(score_type.get_score_column(), required=True)

    experiment_col = get_header_col("experiment")
    charge_col = get_header_col("charge", required=for_quantification)

    intensity_col = get_header_col("intensity", required=for_quantification)
    tmt_cols = get_header_cols_starting_with("reporter intensity ")
    silac_cols = list()
    if "intensity l" in headers:  # SILAC
        silac_cols.append(get_header_col("intensity l", required=for_quantification))
        if "intensity m" in headers:
            silac_cols.append(
                get_header_col("intensity m", required=for_quantification)
            )
        if "intensity h" in headers:
            silac_cols.append(
                get_header_col("intensity h", required=for_quantification)
            )

    raw_file_col = get_header_col("raw file", required=for_quantification)
    fraction_col = get_header_col("fraction")
    evidence_id_col = get_header_col("id", required=for_quantification)

    logger.info("Parsing MaxQuant evidence.txt file")
    for line_idx, row in enumerate(reader):
        if line_idx % 500000 == 0:
            logger.info(f"    Reading line {line_idx}")

        peptide = row[pept_col]
        proteins = get_proteins(
            helpers.clean_peptide(peptide), row[protein_col].split(";")
        )
        if not proteins:
            continue

        experiment = "Experiment1"
        if experiment_col >= 0:
            experiment = row[experiment_col]

        score = float(row[score_col]) if len(row[score_col]) > 0 else float("nan")

        if for_quantification:
            charge = int(row[charge_col])
            intensity = (
                float(row[intensity_col]) if len(row[intensity_col]) > 0 else 0.0
            )
            if fraction_col >= 0:
                fraction = row[fraction_col]
            else:
                fraction = -1
            raw_file = row[raw_file_col]
            tmt_intensities = [row[tmt_col] for tmt_col in tmt_cols]
            silac_intensities = [
                row[silac_col] if len(row[silac_col]) > 0 else 0
                for silac_col in silac_cols
            ]
            evidence_id = int(row[evidence_id_col])
            yield peptide, proteins, charge, raw_file, experiment, fraction, intensity, score, tmt_intensities, silac_intensities, evidence_id
        else:
            yield peptide, proteins, experiment, score


def get_header_col_func(headers):
    def get_header_col(name, required=False):
        if required:
            return get_column_index(headers, name)
        else:
            if name in headers:
                return get_column_index(headers, name)
            return -1

    return get_header_col


def get_header_cols_starting_with_func(headers):
    def get_header_cols_starting_with(name):
        cols = [idx for idx, h in enumerate(headers) if h.startswith(name)]
        return cols

    return get_header_cols_starting_with


def parse_percolator_out_file(reader, headers, get_proteins, score_type="PEP"):
    (
        _,
        pept_col,
        score_col,
        _,
        post_err_prob_col,
        protein_col,
    ) = get_percolator_column_idxs(headers)

    if score_type.get_score_column() == "posterior_error_prob":
        score_col = post_err_prob_col

    logger.info("Parsing Percolator output file")
    for line_idx, row in enumerate(reader):
        if line_idx % 500000 == 0:
            logger.info(f"    Reading line {line_idx}")

        peptide = row[pept_col][1:-1]
        experiment = 1
        score = float(row[score_col])

        if is_native_percolator_file(headers):
            proteins = row[protein_col:]
        elif is_mokapot_file(headers):
            proteins = row[protein_col].split("\t")

        proteins = get_proteins(peptide, proteins)
        if proteins:
            yield peptide, proteins, experiment, score


def is_percolator_file(headers):
    return is_native_percolator_file(headers) or is_mokapot_file(headers)


def is_native_percolator_file(headers):
    return "psmid" in map(str.lower, headers)


def is_mokapot_file(headers):
    return "specid" in map(str.lower, headers)


def get_delimiter(filename: str):
    if filename.endswith(".csv"):
        return ","
    else:
        return "\t"


def get_column_index(headers: List[str], column_name: str):
    if column_name not in headers:
        raise ValueError(
            f"Column {column_name} is missing. Please check your input file."
        )
    return headers.index(column_name)


def get_percolator_column_idxs(headers):
    if is_native_percolator_file(headers):
        id_col = get_column_index(headers, "PSMId")
        pept_col = get_column_index(headers, "peptide")
        score_col = get_column_index(headers, "score")
        qval_col = get_column_index(headers, "q-value")
        post_err_prob_col = get_column_index(headers, "posterior_error_prob")
        protein_col = get_column_index(headers, "proteinIds")
    elif is_mokapot_file(headers):
        id_col = get_column_index(headers, "SpecId")
        pept_col = get_column_index(headers, "Peptide")
        score_col = get_column_index(headers, "mokapot score")
        qval_col = get_column_index(headers, "mokapot q-value")
        post_err_prob_col = get_column_index(headers, "mokapot PEP")
        protein_col = get_column_index(headers, "Proteins")
    else:
        raise ValueError(
            "Could not determine percolator input file format. The file should either "
            "contain a column named PSMId (native percolator) or SpecId (mokapot)."
        )
    return id_col, pept_col, score_col, qval_col, post_err_prob_col, protein_col


def parse_percolator_out_file_to_dict(
    perc_out_file: str, results_dict: dict, input_type: str = ""
):
    if not os.path.isfile(perc_out_file):
        raise FileNotFoundError(
            f"Could not find percolator output file {perc_out_file}. Please check if this file exists."
        )

    delimiter = get_delimiter(perc_out_file)
    reader = get_tsv_reader(perc_out_file, delimiter)
    headers = next(reader)  # save the header

    id_col, pept_col, score_col, _, post_err_prob_col, _ = get_percolator_column_idxs(
        headers
    )

    logger.info("Parsing Percolator output file")
    fixed_mod_idx = -1
    first = True
    for line_idx, row in enumerate(reader):
        if line_idx % 500000 == 0:
            logger.info(f"    Reading line {line_idx}")

        peptide = row[pept_col][2:-2]
        score = float(row[score_col])
        post_err_prob = float(row[post_err_prob_col])

        psm_id = row[id_col]
        if input_type == "prosit":
            peptide = peptide.replace("m", "M(ox)")
            if first:
                for i, fixed_mod in enumerate(FIXED_MODS_UNIMOD):
                    if fixed_mod in peptide:
                        fixed_mod_idx = i
                first = False
            elif fixed_mod_idx >= 0:
                if FIXED_MODS_UNIMOD[fixed_mod_idx] not in peptide:
                    fixed_mod_idx = -1
            raw_file = "-".join(psm_id.split("-")[:-4])
            scan_number = int(float(psm_id.split("-")[-4]))
        else:
            peptide = peptide.replace("[42]", "(ac)").replace("M[16]", "M(ox)")
            raw_file = "_".join(psm_id.split("_")[:-3])
            scan_number = int(psm_id.split("_")[-3])

        results_dict[raw_file][(scan_number, peptide)] = (score, post_err_prob)

    return FIXED_MODS_DICTS[fixed_mod_idx + 1], results_dict


# copied from fundamentals/mod_string.py
def maxquant_to_internal(
    sequences: List[str], fixed_mods: Optional[Dict[str, str]] = {"C": "C[UNIMOD:4]"}
) -> List[str]:
    """
    Function to translate a MaxQuant modstring to the Prosit format
    :param sequences: List[str] of sequences
    :param fixed_mods: Optional dictionary of modifications with key aa and value mod, e.g. 'M': 'M(UNIMOD:35)'.
    Fixed modifications must be included in the variable modificatons dictionary throws Assertion error otherwise.
    :return: List[str] of modified sequences.
    """
    err_msg = f"Provided illegal fixed mod, supported modifications are {set(MAXQUANT_VAR_MODS.values())}."
    assert all(x in MAXQUANT_VAR_MODS.values() for x in fixed_mods.values()), err_msg

    replacements = {**MAXQUANT_VAR_MODS, **fixed_mods}

    def custom_regex_escape(key: str) -> str:
        """
        Subfunction to escape only normal brackets in the modstring
        :param key: The match to escape.
        :return match with escaped special characters.
        """
        for k, v in {"(": "\(", ")": "\)"}.items():
            key = key.replace(k, v)
        return key

    regex = re.compile("|".join(map(custom_regex_escape, replacements.keys())))

    def find_replacement(match: re) -> str:
        """
        Subfunction to find the corresponding substitution for a match.
        :param match: an re.Match object found by re.sub
        :return substitution string for the given match
        """
        key = match.string[match.start() : match.end()]
        if "_" in key:  # If _ is in the match we need to differentiate n and c term
            if match.start() == 0:
                key = f"^{key}"
            else:
                key = f"{key}$"

        return replacements[key]

    return [regex.sub(find_replacement, seq)[1:-1] for seq in sequences]


def get_tsv_reader(filename: str, delimiter: str = "\t"):
    return csv.reader(
        open(filename, "r", newline="", encoding="utf-8-sig"), delimiter=delimiter
    )


def get_tsv_writer(filename: str, delimiter: str = "\t"):
    return csv.writer(open(filename, "w", newline=""), delimiter=delimiter)
