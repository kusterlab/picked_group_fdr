import os
import sys
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


def parseProteinGroupsFiles(
    protein_groups_files: List[str], are_decoy_file: List[bool], **kwargs
):
    for mqProteinGroupsFile, is_decoy_file in zip(protein_groups_files, are_decoy_file):
        yield from parseProteinGroupsFile(
            mqProteinGroupsFile, is_decoy_file=is_decoy_file, **kwargs
        )


def parseProteinGroupsFile(
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
        mqProteinGroupsFile (_type_): _description_
        protein_column (str, optional): _description_. Defaults to 'Protein IDs'.
        score_column (str, optional): _description_. Defaults to 'Score'.
        is_decoy_file (bool, optional): _description_. Defaults to False.

    Yields:
        _type_: _description_
    """
    delimiter = getDelimiter(protein_groups_file)

    reader = getTsvReader(protein_groups_file, delimiter)
    headers = next(reader)  # save the header

    scoreCol = get_column_index(headers, score_column)
    proteinCol = get_column_index(headers, protein_column)

    logger.info(f"Parsing proteinGroups file: {protein_groups_file}")
    for row in reader:
        proteins = list(map(str.strip, row[proteinCol].split(";")))
        if is_decoy_file:
            proteins = [f"REV__{p}" for p in proteins]
        score = -100.0
        if len(row[scoreCol]) > 0:
            score = float(row[scoreCol])
        yield proteins, score


def parsePeptidesFiles(peptides_files: List[str], are_decoy_file: List[bool], **kwargs):
    for peptide_file, is_decoy_file in zip(peptides_files, are_decoy_file):
        yield from parsePeptidesFile(
            peptide_file, is_decoy_file=is_decoy_file, **kwargs
        )


def parsePeptidesFile(
    peptides_file: str,
    peptide_column: str = "Modified sequence",
    protein_column: str = "Protein IDs",
    score_column: str = "Score",
    is_decoy_file: bool = False,
):
    """Parse protein groups file from MaxQuant or ProteomeDiscoverer.

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
    delimiter = getDelimiter(peptides_file)

    reader = getTsvReader(peptides_file, delimiter)
    headers = next(reader)  # save the header

    peptideCol = get_column_index(headers, peptide_column)
    scoreCol = get_column_index(headers, score_column)
    if not is_decoy_file:
        proteinCol = get_column_index(headers, protein_column)

    experiment = "Experiment1"

    logger.info(f"Parsing peptides file: {peptides_file}")
    for row in reader:
        proteins = []
        if is_decoy_file:
            proteins = ["REV__protein"]
        elif len(row[proteinCol]) > 0:
            proteins = list(map(str.strip, row[proteinCol].split(";")))

        score = -100.0
        if len(row[scoreCol]) > 0:
            score = float(row[scoreCol])
            if np.isnan(score):
                continue
        yield row[peptideCol], proteins, experiment, score


def parseEvidenceFiles(
    evidenceFiles,
    peptideToProteinMaps,
    scoreType,
    forQuantification=False,
    suppressMissingPeptideWarning=False,
):
    for evidenceFile, peptideToProteinMap in zip(evidenceFiles, peptideToProteinMaps):
        yield from parseEvidenceFile(
            evidenceFile,
            peptideToProteinMap,
            scoreType,
            forQuantification,
            suppressMissingPeptideWarning,
        )


def parseEvidenceFile(
    evidenceFile,
    peptideToProteinMap,
    scoreType,
    forQuantification=False,
    suppressMissingPeptideWarning=False,
):
    delimiter = getDelimiter(evidenceFile)
    reader = getTsvReader(evidenceFile, delimiter)
    headers = next(reader)  # save the header

    get_proteins = get_peptide_to_protein_mapper(
        peptideToProteinMap, scoreType, suppressMissingPeptideWarning
    )

    if isPercolatorFile(headers):
        yield from parsePercolatorOutFile(reader, headers, get_proteins, scoreType)
    else:
        headers = list(
            map(str.lower, headers)
        )  # convert headers to lowercase since MQ changes the capitalization frequently
        if evidenceFile.endswith(".csv"):
            headers = [x.replace(".", " ") for x in headers]
        yield from parseMqEvidenceFile(
            reader, headers, get_proteins, scoreType, forQuantification
        )


def get_peptide_to_protein_mapper(
    peptideToProteinMap, scoreType, suppressMissingPeptideWarning
):
    def get_proteins(peptide, tmp_proteins):
        if scoreType.remaps_peptides_to_proteins():
            proteins = digest.getProteins(peptideToProteinMap, peptide)
            if len(proteins) == 0:
                if (
                    not helpers.isContaminant(tmp_proteins)
                    and not suppressMissingPeptideWarning
                ):
                    logger.warning(f"Missing peptide: {peptide} {tmp_proteins}")
                return None
        else:
            proteins = tmp_proteins
        proteins = scoreType.filter_proteins(
            proteins
        )  # filtering for razor peptide approach
        return helpers.removeDecoyProteinsFromTargetPeptides(proteins)

    return get_proteins


def parseMqEvidenceFile(
    reader, headers, get_proteins, scoreType, forQuantification=False
):
    """
    Reads in approximately 100,000 lines per second with forQuantification=False
    and 50,000 lines per second with forQuantification=True

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
    getHeaderCol = getHeaderColFunc(headers)
    getHeaderColsStartingWith = getHeaderColsStartingWithFunc(headers)

    peptCol = getHeaderCol("modified sequence", required=True)

    proteinCol = getHeaderCol(
        "leading proteins", required=True
    )  # all protein groups, each represented by the first protein in the group
    if scoreType.use_razor:
        proteinCol = getHeaderCol(
            "leading razor protein", required=True
        )  # best scoring protein group, represented by the first protein in the group

    scoreCol = getHeaderCol(scoreType.get_score_column(), required=True)

    experimentCol = getHeaderCol("experiment")
    chargeCol = getHeaderCol("charge", required=forQuantification)

    intensityCol = getHeaderCol("intensity", required=forQuantification)
    tmtCols = getHeaderColsStartingWith("reporter intensity ")
    silacCols = list()
    if "intensity l" in headers:  # SILAC
        silacCols.append(getHeaderCol("intensity l", required=forQuantification))
        if "intensity m" in headers:
            silacCols.append(getHeaderCol("intensity m", required=forQuantification))
        if "intensity h" in headers:
            silacCols.append(getHeaderCol("intensity h", required=forQuantification))

    rawFileCol = getHeaderCol("raw file", required=forQuantification)
    fractionCol = getHeaderCol("fraction")
    evidenceIdCol = getHeaderCol("id", required=forQuantification)

    logger.info("Parsing MaxQuant evidence.txt file")
    for lineIdx, row in enumerate(reader):
        if lineIdx % 500000 == 0:
            logger.info(f"    Reading line {lineIdx}")

        peptide = row[peptCol]
        proteins = get_proteins(
            helpers.cleanPeptide(peptide), row[proteinCol].split(";")
        )
        if not proteins:
            continue

        experiment = "Experiment1"
        if experimentCol >= 0:
            experiment = row[experimentCol]

        score = float(row[scoreCol]) if len(row[scoreCol]) > 0 else float("nan")

        if forQuantification:
            charge = int(row[chargeCol])
            intensity = float(row[intensityCol]) if len(row[intensityCol]) > 0 else 0.0
            if fractionCol >= 0:
                fraction = row[fractionCol]
            else:
                fraction = -1
            rawFile = row[rawFileCol]
            tmtIntensities = [row[tmtCol] for tmtCol in tmtCols]
            silacIntensities = [
                row[silacCol] if len(row[silacCol]) > 0 else 0 for silacCol in silacCols
            ]
            evidenceId = int(row[evidenceIdCol])
            yield peptide, proteins, charge, rawFile, experiment, fraction, intensity, score, tmtIntensities, silacIntensities, evidenceId
        else:
            yield peptide, proteins, experiment, score


def getHeaderColFunc(headers):
    def getHeaderCol(name, required=False):
        if required:
            return get_column_index(headers, name)
        else:
            if name in headers:
                return get_column_index(headers, name)
            return -1

    return getHeaderCol


def getHeaderColsStartingWithFunc(headers):
    def getHeaderColsStartingWith(name):
        cols = [idx for idx, h in enumerate(headers) if h.startswith(name)]
        return cols

    return getHeaderColsStartingWith


def parsePercolatorOutFile(reader, headers, get_proteins, scoreType="PEP"):
    _, peptCol, scoreCol, _, postErrProbCol, proteinCol = getPercolatorColumnIdxs(
        headers
    )

    if scoreType.get_score_column() == "posterior_error_prob":
        scoreCol = postErrProbCol

    logger.info("Parsing Percolator output file")
    for lineIdx, row in enumerate(reader):
        if lineIdx % 500000 == 0:
            logger.info(f"    Reading line {lineIdx}")

        peptide = row[peptCol][1:-1]
        experiment = 1
        score = float(row[scoreCol])

        if isNativePercolatorFile(headers):
            proteins = row[proteinCol:]
        elif isMokapotFile(headers):
            proteins = row[proteinCol].split("\t")

        proteins = get_proteins(peptide, proteins)
        if proteins:
            yield peptide, proteins, experiment, score


def getPercolatorNativeHeaders():
    return [
        "PSMId",
        "score",
        "q-value",
        "posterior_error_prob",
        "peptide",
        "proteinIds",
    ]


def isPercolatorFile(headers):
    return isNativePercolatorFile(headers) or isMokapotFile(headers)


def isNativePercolatorFile(headers):
    return "psmid" in map(str.lower, headers)


def isMokapotFile(headers):
    return "specid" in map(str.lower, headers)


def getDelimiter(filename: str):
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


def getPercolatorColumnIdxs(headers):
    if isNativePercolatorFile(headers):
        idCol = get_column_index(headers, "PSMId")
        peptCol = get_column_index(headers, "peptide")
        scoreCol = get_column_index(headers, "score")
        qvalCol = get_column_index(headers, "q-value")
        postErrProbCol = get_column_index(headers, "posterior_error_prob")
        proteinCol = get_column_index(headers, "proteinIds")
    elif isMokapotFile(headers):
        idCol = get_column_index(headers, "SpecId")
        peptCol = get_column_index(headers, "Peptide")
        scoreCol = get_column_index(headers, "mokapot score")
        qvalCol = get_column_index(headers, "mokapot q-value")
        postErrProbCol = get_column_index(headers, "mokapot PEP")
        proteinCol = get_column_index(headers, "Proteins")
    else:
        raise ValueError(
            "Could not determine percolator input file format. The file should either contain a column named PSMId (native percolator) or SpecId (mokapot)."
        )
    return idCol, peptCol, scoreCol, qvalCol, postErrProbCol, proteinCol


def parsePercolatorOutFileToDict(percOutFile, resultsDict, inputType=""):
    if not os.path.isfile(percOutFile):
        raise FileNotFoundError(
            f"Could not find percolator output file {percOutFile}. Please check if this file exists."
        )

    delimiter = getDelimiter(percOutFile)
    reader = getTsvReader(percOutFile, delimiter)
    headers = next(reader)  # save the header

    idCol, peptCol, scoreCol, _, postErrProbCol, _ = getPercolatorColumnIdxs(headers)

    logger.info("Parsing Percolator output file")
    fixed_mod_idx = -1
    first = True
    for lineIdx, row in enumerate(reader):
        if lineIdx % 500000 == 0:
            logger.info(f"    Reading line {lineIdx}")

        peptide = row[peptCol][2:-2]
        score = float(row[scoreCol])
        postErrProb = float(row[postErrProbCol])

        psmId = row[idCol]
        if inputType == "prosit":
            peptide = peptide.replace("m", "M(ox)")
            if first:
                for i, fixed_mod in enumerate(FIXED_MODS_UNIMOD):
                    if fixed_mod in peptide:
                        fixed_mod_idx = i
                first = False
            elif fixed_mod_idx >= 0:
                if FIXED_MODS_UNIMOD[fixed_mod_idx] not in peptide:
                    fixed_mod_idx = -1
            rawFile = "-".join(psmId.split("-")[:-4])
            scanNumber = int(float(psmId.split("-")[-4]))
        else:
            peptide = peptide.replace("[42]", "(ac)").replace("M[16]", "M(ox)")
            rawFile = "_".join(psmId.split("_")[:-3])
            scanNumber = int(psmId.split("_")[-3])

        resultsDict[rawFile][(scanNumber, peptide)] = (score, postErrProb)

    return FIXED_MODS_DICTS[fixed_mod_idx + 1], resultsDict


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


def getTsvReader(filename, delimiter="\t"):
    # Python 3
    if sys.version_info[0] >= 3:
        return csv.reader(open(filename, "r", newline=""), delimiter=delimiter)
    # Python 2
    else:
        return csv.reader(open(filename, "rb"), delimiter=delimiter)


def getTsvWriter(filename, delimiter="\t"):
    # Python 3
    if sys.version_info[0] >= 3:
        return csv.writer(open(filename, "w", newline=""), delimiter=delimiter)
    # Python 2
    else:
        return csv.writer(open(filename, "wb"), delimiter=delimiter)
