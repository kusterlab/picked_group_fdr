import re
from typing import Dict, List, Optional


TMT_UNIMOD = "[UNIMOD:737]"  # TMT6/10/11
TMTPRO_UNIMOD = "[UNIMOD:2016]"
ITRAQ4_UNIMOD = "[UNIMOD:214]"
ITRAQ8_UNIMOD = "[UNIMOD:730]"


DEFAULT_FIXED_MODS = {"C": "C[UNIMOD:4]"}
SILAC_HEAVY_FIXED_MODS = {
    "C": "C[UNIMOD:4]",
    "K": "K[UNIMOD:259]",
    "R": "R[UNIMOD:267]",
}

TMT_FIXED_MODS = {
    "C": "C[UNIMOD:4]",
    "^_": f"_{TMT_UNIMOD}-",
    "K": f"K{TMT_UNIMOD}",
}
TMTPRO_FIXED_MODS = {
    "C": "C[UNIMOD:4]",
    "^_": f"_{TMTPRO_UNIMOD}-",
    "K": f"K{TMTPRO_UNIMOD}",
}
ITRAQ4_FIXED_MODS = {
    "C": "C[UNIMOD:4]",
    "^_": f"_{ITRAQ4_UNIMOD}-",
    "K": f"K{ITRAQ4_UNIMOD}",
}
ITRAQ8_FIXED_MODS = {
    "C": "C[UNIMOD:4]",
    "^_": f"_{ITRAQ8_UNIMOD}-",
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
# copied from fundamentals/constants.py
MAXQUANT_VAR_MODS = {
    "(ox)": "[UNIMOD:35]",
    "(Oxidation (M))": "[UNIMOD:35]",
    "(tm)": "[UNIMOD:737]",
    "^_(tm)": f"_{TMT_UNIMOD}-",
    "K(tm)": f"K{TMT_UNIMOD}",
    "^_(TMTPro (N-term))": f"_{TMTPRO_UNIMOD}-",
    "K(TMTPro (K))": f"K{TMTPRO_UNIMOD}",
    "^_(iTRAQ4plex (N-term))": f"_{ITRAQ4_UNIMOD}-",
    "K(iTRAQ4plex (K))": f"K{ITRAQ4_UNIMOD}",
    "^_(iTRAQ8plex (N-term))": f"_{ITRAQ8_UNIMOD}-",
    "K(iTRAQ8plex (K))": f"K{ITRAQ8_UNIMOD}",
    "(ph)": "[UNIMOD:21]",
    "K(Lys8)": "K[UNIMOD:259]",
    "R(Arg10)": "R[UNIMOD:267]",
    "C(Carbamidomethyl (C))": "C[UNIMOD:4]",
}


# adapted from fundamentals/mod_string.py
def maxquant_mod_to_proforma(
    fixed_mods: Optional[Dict[str, str]] = {"C": "C[UNIMOD:4]"}
) -> List[str]:
    """
    Function to translate a MaxQuant modstring to the ProForma format
    :param sequence: modified sequence in MaxQuant format, e.g. _(ac)ACM(ox)PEPTIDE_
    :param fixed_mods: Optional dictionary of modifications with key aa and value mod, e.g. 'C': 'C[UNIMOD:4]'.
    Fixed modifications must be included in the variable modificatons dictionary throws Assertion error otherwise.
    :return: modified sequence in proforma notation, e.g. [UNIMOD:1]-AC[UNIMOD:4]M[UNIMOD:35]PEPTIDE
    """
    err_msg = f"Provided illegal fixed mod, supported modifications are {set(MAXQUANT_VAR_MODS.values())}."
    assert all(x in MAXQUANT_VAR_MODS.values() for x in fixed_mods.values()), err_msg

    replacements = {**MAXQUANT_VAR_MODS, **fixed_mods}
    replacement_values = list(replacements.values())

    def custom_regex_escape(key: str) -> str:
        """
        Subfunction to escape only normal brackets in the modstring
        :param key: The match to escape.
        :return match with escaped special characters.
        """
        for k, v in {"(": "\\(", ")": "\\)"}.items():
            key = key.replace(k, v)
        return f"({key})"

    regex = re.compile("|".join(map(custom_regex_escape, replacements.keys())))

    def first_not_none_index(tup):
        return next((i for i, x in enumerate(tup) if x is not None), -1)

    def find_replacement(match: re) -> str:
        """
        Subfunction to find the corresponding substitution for a match.
        :param match: an re.Match object found by re.sub
        :return substitution string for the given match
        """
        value_idx = first_not_none_index(match.groups())
        return replacement_values[value_idx]

    def replace(sequence: str):
        return regex.sub(find_replacement, sequence)[1:-1]
    
    return replace


PROSIT_VAR_MODS = {
    "m": "M[UNIMOD:35]",
    f"^{TMT_UNIMOD}-?": f"{TMT_UNIMOD}-",
    f"^{TMTPRO_UNIMOD}-?": f"{TMTPRO_UNIMOD}-",
    f"^{ITRAQ4_UNIMOD}-?": f"{ITRAQ4_UNIMOD}-",
    f"^{ITRAQ8_UNIMOD}-?": f"{ITRAQ8_UNIMOD}-",
}


def prosit_mod_to_proforma(
    fixed_mods: Optional[Dict[str, str]] = None,
    variable_mods: Optional[Dict[str, str]] = None,
) -> List[str]:
    """
    Function to translate a Prosit modstring to the ProForma format
    :param sequence: modified sequence in (old) Prosit format, e.g. [UNIMOD:1]AC[UNIMOD:4]mPEPTIDE
    :param fixed_mods: Optional dictionary of modifications with key aa and value mod, e.g. 'C': 'C[UNIMOD:4]'.
    Fixed modifications must be included in the variable modificatons dictionary throws Assertion error otherwise.
    :return: modified sequence in proforma notation, e.g. [UNIMOD:1]-AC[UNIMOD:4]M[UNIMOD:35]PEPTIDE
    """
    if fixed_mods is None:
        fixed_mods = {}

    if variable_mods is None:
        variable_mods = {}

    err_msg = f"Provided illegal fixed mod, supported modifications are {set(PROSIT_VAR_MODS.values())}."
    assert all(x in PROSIT_VAR_MODS.values() for x in fixed_mods.values()), err_msg

    replacements = {**PROSIT_VAR_MODS, **fixed_mods, **variable_mods}
    replacement_values = list(replacements.values())

    def custom_regex_escape(key: str) -> str:
        """
        Subfunction to escape only square brackets in the modstring
        :param key: The match to escape.
        :return match with escaped special characters.
        """
        for k, v in {"[": "\\[", "]": "\\]"}.items():
            key = key.replace(k, v)
        return f"({key})"

    regex = re.compile("|".join(map(custom_regex_escape, replacements.keys())))

    def first_not_none_index(tup):
        return next((i for i, x in enumerate(tup) if x is not None), -1)

    def find_replacement(match: re) -> str:
        """
        Subfunction to find the corresponding substitution for a match.
        :param match: an re.Match object found by re.sub
        :return substitution string for the given match
        """
        value_idx = first_not_none_index(match.groups())
        return replacement_values[value_idx]

    def replace(sequence: str):
        return regex.sub(find_replacement, sequence)
    
    return replace
