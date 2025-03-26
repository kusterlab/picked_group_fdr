from __future__ import annotations

import _csv
from dataclasses import dataclass
from enum import Enum
import logging
from typing import List, Iterator

from . import tsv
from .. import results  # TODO: get rid of this import
from .. import writers  # TODO: get rid of this import

# for type hints only
from .. import scoring_strategy


logger = logging.getLogger(__name__)


def parse_mq_evidence_file(
    reader,
    headers,
    get_proteins,
    score_type: scoring_strategy.ProteinScoringStrategy,
    for_quantification: bool = False,
    **kwargs,
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
    # convert headers to lowercase since MQ changes the capitalization frequently
    headers = list(map(str.lower, headers))

    get_header_col = tsv.get_header_col_func(headers)
    get_header_cols_starting_with = tsv.get_header_cols_starting_with_func(headers)

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
    silac_cols = get_silac_cols(headers, get_header_col, for_quantification)

    raw_file_col = get_header_col("raw file", required=for_quantification)
    fraction_col = get_header_col("fraction")
    evidence_id_col = get_header_col("id", required=for_quantification)

    logger.info("Parsing MaxQuant evidence.txt file")
    for line_idx, row in enumerate(reader):
        if line_idx % 500000 == 0:
            logger.info(f"    Reading line {line_idx}")

        modified_peptide = row[pept_col][1:-1]
        proteins = get_proteins(modified_peptide, row[protein_col].split(";"))
        if not proteins:
            continue

        experiment = "Experiment1"
        if experiment_col >= 0:
            experiment = row[experiment_col]

        score = float(row[score_col]) if len(row[score_col]) > 0 else float("nan")

        if not for_quantification:
            yield modified_peptide, proteins, experiment, score
        else:
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
            yield modified_peptide, proteins, charge, raw_file, experiment, fraction, intensity, score, tmt_intensities, silac_intensities, evidence_id


def get_silac_cols(headers, get_header_col, for_quantification):
    silac_cols = list()
    if "intensity l" in headers:
        silac_cols.append(get_header_col("intensity l", required=for_quantification))
        if "intensity m" in headers:
            silac_cols.append(
                get_header_col("intensity m", required=for_quantification)
            )
        if "intensity h" in headers:
            silac_cols.append(
                get_header_col("intensity h", required=for_quantification)
            )
    return silac_cols


def parse_mq_protein_groups_file(
    mq_protein_groups_file: str, additional_headers: List[str] = None
) -> results.ProteinGroupResults:
    if additional_headers is None:
        additional_headers = []

    delimiter = "\t"
    if mq_protein_groups_file.endswith(".csv"):
        delimiter = ","

    reader = tsv.get_tsv_reader(mq_protein_groups_file, delimiter)
    headers = next(reader)

    cols = {
        x: headers.index(x)
        for x in writers.PROTEIN_GROUP_HEADERS + additional_headers
        if x in headers
    }

    logger.info("Parsing MaxQuant proteinGroups.txt file")
    protein_group_results = []
    for row in reader:
        protein_group_results.append(
            parse_mq_protein_groups_file_row(row, cols, additional_headers)
        )

    protein_group_results = results.ProteinGroupResults(protein_group_results)
    protein_group_results.append_headers(additional_headers)

    return protein_group_results


def parse_mq_protein_groups_file_row(row, cols, additional_headers):
    _get_field = lambda x: row[cols[x]] if x in cols else ""

    return results.ProteinGroupResult(
        proteinIds=_get_field("Protein IDs"),
        majorityProteinIds=_get_field("Majority protein IDs"),
        peptideCountsUnique=_get_field("Peptide counts (unique)"),
        bestPeptide="",
        numberOfProteins=int(_get_field("Number of proteins")),
        qValue=float(_get_field("Q-value")),
        score=float(_get_field("Score")),
        reverse=_get_field("Reverse"),
        potentialContaminant=_get_field("Potential contaminant"),
        precursorQuants=[],
        extraColumns=[_get_field(x) for x in additional_headers],
    )


class LabelingState(Enum):
    NOT_SILAC = -3  # missing labeling state column
    UNKNOWN = -2  # empty value
    LIGHT_MAYBE = (
        -1
    )  # cannot figure out what this actually means, but it seems to be light as well...
    LIGHT = 0
    HEAVY = 1


@dataclass
class EvidenceRow:
    raw_file: str
    scannr: int
    peptide: str
    score: float
    post_err_prob: float
    is_decoy: bool
    is_contaminant: bool
    id_type: str
    labeling_state: str


def parse_evidence_file_for_percolator_matching(
    reader: "_csv._reader", headers: List[str]
) -> Iterator[List[str], EvidenceRow]:
    """Parses MaxQuant evidence.txt or msms.txt file such that the score and 
    PEP can be updated using percolator output files.

    Args:
        reader (_csv._reader): CSV reader object for an evidence file
        headers (List[str]): headers of the evidence file, all in lower case

    Yields:
        Iterator[List[str], EvidenceRow]: tuples of (csv reader row, Evidence row object)
    
    Columns needed for percolator matching:
    - Raw file
    - Scan number
    - Score
    - PEP
    - Modified sequence
    - Reverse
    - Potential contaminant
    - Labeling state (for SILAC)
    """    
    score_col = tsv.get_column_index(headers, "score")
    post_err_prob_col = tsv.get_column_index(headers, "pep")

    # these columns are needed to retrieve the PSM
    raw_file_col = tsv.get_column_index(headers, "raw file")
    if "ms/ms scan number" in headers:
        scannr_col = tsv.get_column_index(headers, "ms/ms scan number")  # evidence.txt
    else:
        scannr_col = tsv.get_column_index(headers, "scan number")  # msms.txt
    pept_col = tsv.get_column_index(headers, "modified sequence")
    id_type_col = tsv.get_column_index(
        headers, "type"
    )  # MULTI-MSMS MULTI-MATCH MSMS MULTI-SECPEP MULTI-MATCH-MSMS
    reverse_col = tsv.get_column_index(headers, "reverse")
    contaminant_col = tsv.get_column_index(headers, "potential contaminant")
    labeling_state_col = None
    if "labeling state" in headers:
        labeling_state_col = tsv.get_column_index(headers, "labeling state")

    for row in reader:
        # if scanNr is empty, it is an MBR evidence row which we encode with scanNr = -1
        scanNr = -1
        if len(row[scannr_col]) > 0:
            scanNr = int(row[scannr_col])

        labeling_state = LabelingState.NOT_SILAC
        if labeling_state_col is not None:
            labeling_state = LabelingState.UNKNOWN
            if len(row[labeling_state_col]) > 0:
                labeling_state = int(row[labeling_state_col])

        yield row, EvidenceRow(
            raw_file=row[raw_file_col],
            scannr=scanNr,
            score=float(row[score_col]),
            post_err_prob=float(row[post_err_prob_col]),
            peptide=row[pept_col],
            is_decoy=(row[reverse_col] == "+"),
            is_contaminant=(row[contaminant_col] == "+"),
            id_type=row[id_type_col],
            labeling_state=labeling_state,
        )


def parse_peptides_file_for_percolator_matching(
    reader: "_csv._reader", headers: List[str]
) -> Iterator[List[str], EvidenceRow]:
    score_col = tsv.get_column_index(headers, "score")
    post_err_prob_col = tsv.get_column_index(headers, "pep")

    # these columns are needed to retrieve the PSM
    pept_col = tsv.get_column_index(headers, "sequence")
    reverse_col = tsv.get_column_index(headers, "reverse")
    contaminant_col = tsv.get_column_index(headers, "potential contaminant")

    for row in reader:
        yield row, EvidenceRow(
            raw_file="",
            scannr=-1,
            score=float(row[score_col]),
            post_err_prob=float(row[post_err_prob_col]),
            peptide="_" + row[pept_col] + "_",
            is_decoy=(row[reverse_col] == "+"),
            is_contaminant=(row[contaminant_col] == "+"),
            id_type="Unknown",
            labeling_state=LabelingState.NOT_SILAC,
        )


def has_unknown_silac_label(labeling_state: LabelingState):
    return (
        labeling_state == LabelingState.UNKNOWN
        or labeling_state == LabelingState.LIGHT_MAYBE
    )


def is_heavy_labeled(labeling_state: LabelingState):
    return labeling_state == LabelingState.HEAVY
