import logging
import re
from typing import Dict, List, Optional


from .. import helpers
from .modifications import ITRAQ4_UNIMOD, ITRAQ8_UNIMOD, TMT_UNIMOD, TMTPRO_UNIMOD
from .tsv import (
    get_header_col_func,
    get_header_cols_starting_with_func,
)

logger = logging.getLogger(__name__)

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
