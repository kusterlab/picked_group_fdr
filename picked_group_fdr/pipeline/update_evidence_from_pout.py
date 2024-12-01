import _csv
import sys
import os
import collections
import logging
from typing import Dict, List, Tuple, Set, Callable, Iterator
from dataclasses import dataclass, field

from .. import __version__, __copyright__
from .. import helpers
from ..picked_group_fdr import ArgumentParserWithLogger
from ..parsers import tsv
from ..parsers import modifications
from ..parsers import percolator
from ..parsers import maxquant

# hacky way to get package logger when running as module
logger = logging.getLogger(__package__ + "." + __file__)


def parse_args(argv):
    import argparse

    apars = ArgumentParserWithLogger(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    apars.add_argument(
        "--mq_evidence",
        default=None,
        metavar="EV",
        nargs="+",
        required=True,
        help="""MaxQuant evidence.txt or msms.txt file(s). If you 
                want to combine multiple evidence files, use 
                spaces to separate the file paths.""",
    )

    apars.add_argument(
        "--mq_evidence_out",
        default=None,
        metavar="EV",
        required=True,
        help="""MaxQuant evidence.txt or msms.txt combined output 
                file with updated percolator results.""",
    )

    apars.add_argument(
        "--perc_results",
        default=list(),
        metavar="POUT",
        nargs="+",
        help="""Percolator output file(s) with PSMs separated by 
                spaces. Also include the decoy results. If this 
                flag is not set, the evidence files are simply 
                concatenated without updating the PSMs.""",
    )

    apars.add_argument(
        "--pout_input_type",
        default="andromeda",
        metavar="SE",
        help="""Input type of percolator output file. Can be 
                "andromeda" or "prosit".""",
    )

    apars.add_argument(
        "--mq_input_type",
        default="auto",
        metavar="IT",
        help="""Input type of MaxQuant output file. Can be 
                "peptides", "evidence", "msms" or "auto" 
                (automatic detection).""",
    )

    apars.add_argument(
        "--suppress_missing_peptide_warning",
        help="Suppress missing peptide warning when mapping peptides to proteins.",
        action="store_true",
    )

    # ------------------------------------------------
    args = apars.parse_args(argv)

    return args


def main(argv) -> None:
    logger.info(f"UpdateEvidence version {__version__}\n{__copyright__}")
    logger.info(
        f'Issued command: {os.path.basename(__file__)} {" ".join(map(str, argv))}'
    )

    args = parse_args(argv)

    if os.path.isfile(args.mq_evidence_out):
        logger.info(
            f"Found updated evidence file {args.mq_evidence_out}, remove this file to rerun update_evidence."
        )
        return

    update_evidence_files(
        args.mq_evidence,
        args.perc_results,
        args.mq_evidence_out,
        args.mq_input_type,
        args.pout_input_type,
        args.suppress_missing_peptide_warning,
    )


def update_evidence_files(
    evidence_files: List[str],
    pout_files: List[str],
    out_evidence_file: str,
    mq_input_type: str,
    pout_input_type: str,
    suppress_missing_peptide_warning: bool,
) -> None:
    fixed_mods, results_dict = get_percolator_results(pout_files, pout_input_type)
    parser_func = maxquant.parse_evidence_file_for_percolator_matching
    if mq_input_type == "peptides":
        fixed_mods = dict()
        results_dict = convert_PSM_dict_to_peptide_dict(results_dict)
        parser_func = maxquant.parse_peptides_file_for_percolator_matching

    logger.info("Writing updated combined evidence file")
    writer = tsv.get_tsv_writer(out_evidence_file + ".tmp")
    first_headers = []
    for evidence_file in evidence_files:
        first_headers = update_evidence_single(
            evidence_file,
            writer,
            first_headers,
            fixed_mods,
            results_dict,
            parser_func,
            pout_input_type,
            suppress_missing_peptide_warning,
        )

    os.rename(out_evidence_file + ".tmp", out_evidence_file)

    logger.info(f"Results written to {out_evidence_file}")


@dataclass
class EvidenceUpdateStats:
    mq_peps: List[Tuple[float, str]] = field(default_factory=list)
    prosit_peps: List[Tuple[float, str]] = field(default_factory=list)
    unexplained_missing_psms: int = 0
    unexplained_peptides: List[str] = field(default_factory=list)
    rows_written: int = 0
    last_written_row: int = -1

    def log_progress(self):
        if (
            self.rows_written % 500000 == 0
            and self.rows_written != self.last_written_row
        ):
            logger.info(f"    Writing line {self.rows_written}")
            self.last_written_row = self.rows_written

    def add_missing_peptide(self, peptide: str) -> None:
        if self.unexplained_missing_psms < 10:
            self.unexplained_peptides.append(peptide)
        self.unexplained_missing_psms += 1

    def calculate_unexplained_percentage(self) -> int:
        if self.rows_written > 0:
            return int(self.unexplained_missing_psms / self.rows_written * 100)
        return 0

    def log_unexplained_peptides_summary(self) -> None:
        unexplained_percentage = self.calculate_unexplained_percentage()
        logger.info(
            f"Unexplained missing PSMs in Percolator results: {self.unexplained_missing_psms} out of {self.rows_written} ({unexplained_percentage}%)"
        )
        if self.unexplained_missing_psms > 0:
            logger.info(
                "If this percentage is low (<5%), it is probably the result of second peptide hits."
            )
            logger.info("\tFirst 10 missing peptides:")
            for peptide in self.unexplained_peptides:
                logger.info("\t" + peptide)

    def log_summary(self, pout_input_type: str) -> None:
        self.log_unexplained_peptides_summary()

        logger.info("#MQ identifications:")
        self.count_below_FDR(self.mq_peps)

        logger.info(f"#{pout_input_type} identifications:")
        self.count_below_FDR(self.prosit_peps)

    @staticmethod
    def count_below_FDR(
        post_err_probs: List[Tuple[float, str]], fdr_threshold: float = 0.01
    ):
        post_err_probs = sorted(post_err_probs)
        counts = collections.defaultdict(int)
        summed_PEP = 0.0
        for i, (pep, id_type) in enumerate(post_err_probs):
            summed_PEP += pep
            if summed_PEP / (i + 1) > fdr_threshold:
                for k, v in sorted(counts.items()):
                    logger.info(f"  {k}: {v}")
                break
            counts[id_type] += 1


def update_evidence_single(
    evidence_file: str,
    writer: "_csv._writer",
    first_headers: List[str],
    fixed_mods: Dict[str, str],
    results_dict: percolator.ResultsDict,
    parser_func: Callable[
        ["_csv._reader", List[str]], Iterator[Tuple[List[str], maxquant.EvidenceRow]]
    ],
    pout_input_type: str,
    suppress_missing_peptide_warning: bool,
) -> List[str]:
    logger.info(f"Processing {evidence_file}")

    reader = tsv.get_tsv_reader(evidence_file)
    headers = initialize_headers(reader, first_headers, writer)

    # these columns can be updated
    score_col = tsv.get_column_index(headers, "score")
    post_err_prob_col = tsv.get_column_index(headers, "pep")

    stats = EvidenceUpdateStats()
    missing_raw_files = set()

    for row, psm in parser_func(reader, headers):
        stats.log_progress()

        if skip_row_due_to_mbr_or_empty_results_dict(psm, results_dict):
            write_row(writer, row, stats)
            continue

        if raw_file_missing(psm, results_dict, missing_raw_files):
            continue

        if not psm.is_decoy:
            stats.mq_peps.append((psm.post_err_prob, psm.id_type))

        perc_result, peptide = find_percolator_psm(
            psm, fixed_mods, results_dict, pout_input_type
        )

        if not perc_result:
            if is_unexplainable_missing_psm(psm, peptide, pout_input_type):
                stats.add_missing_peptide(psm.peptide)
                if not suppress_missing_peptide_warning:
                    logger.debug(
                        f"Missing PSM: {psm.raw_file}, {peptide}, {psm.scannr}"
                    )
            continue

        process_percolator_result(
            psm, perc_result, score_col, post_err_prob_col, row, stats
        )
        write_row(writer, row, stats)

    stats.log_summary(pout_input_type)
    return headers


def initialize_headers(
    reader: "_csv._reader", first_headers: List[str], writer: "_csv._writer"
) -> List[str]:
    headers_original = next(reader)
    headers = [header.lower() for header in headers_original]

    if len(first_headers) == 0:
        writer.writerow(headers_original)
    elif headers != first_headers:
        warn_for_header_difference(first_headers, headers)

    return headers


def skip_row_due_to_mbr_or_empty_results_dict(
    psm: maxquant.EvidenceRow, results_dict: percolator.ResultsDict
) -> bool:
    return len(results_dict) == 0 or is_mbr_evidence_row(psm)


def raw_file_missing(
    psm: maxquant.EvidenceRow,
    results_dict: percolator.ResultsDict,
    missing_raw_files: Set[str],
) -> bool:
    if len(results_dict.get(psm.raw_file, {})) == 0:
        if psm.raw_file not in missing_raw_files:
            logger.warning(
                f"Found no PSMs for {psm.raw_file} in percolator result files"
            )
            missing_raw_files.add(psm.raw_file)
        return True
    return False


def process_percolator_result(
    psm: maxquant.EvidenceRow,
    perc_result: percolator.ScorePEPPair,
    score_col: int,
    post_err_prob_col: int,
    row: List[str],
    stats: EvidenceUpdateStats,
) -> None:
    perc_score, perc_post_err_prob = perc_result

    if not psm.is_decoy:
        stats.prosit_peps.append((perc_post_err_prob, psm.id_type))

    row[score_col] = perc_score
    row[post_err_prob_col] = perc_post_err_prob


def write_row(
    writer: "_csv._writer", row: List[str], stats: EvidenceUpdateStats
) -> None:
    writer.writerow(row)
    stats.rows_written += 1


def warn_for_header_difference(first_headers: List[str], headers: List[str]):
    logger.warning("Current column names are different from the first evidence file")
    logger.info("Column\tFirstHeaders\tCurrentHeaders")
    logger.info(
        "\n".join(
            [
                str(i + 1) + "\t" + x + "\t" + y
                for i, (x, y) in enumerate(zip(first_headers, headers))
                if x != y
            ]
        )
    )


def is_unexplainable_missing_psm(
    psm: maxquant.EvidenceRow, peptide: str, pout_input_type: str
):
    return (
        not psm.is_contaminant
        and psm.score > 0.0
        and not (pout_input_type == "prosit" and not is_valid_prosit_peptide(peptide))
    )


def find_percolator_psm(
    psm: maxquant.EvidenceRow,
    fixed_mods: Dict[str, str],
    results_dict: percolator.ResultsDict,
    pout_input_type: str,
) -> Tuple[percolator.ScorePEPPair, str]:
    """Finds percolator PSM corresponding to a row in the evidence.txt file.

    Special care is needed for SILAC experiments. Prosit explicitly annotates
    these modifications in the peptide sequence, but MaxQuant does not.
    Moreover, MaxQuant does not annotate which peptides are heavy labeled if
    there are multiple MS/MS scans mapping to the precursor. We therefore
    test for both heavy and light labeled versions of the peptide.
    """
    peptide = psm.peptide[1:-1]
    if pout_input_type == "prosit":
        fixed_mods_tmp = fixed_mods
        if maxquant.is_heavy_labeled(psm.labeling_state):
            fixed_mods_tmp = modifications.SILAC_HEAVY_FIXED_MODS
        peptide = modifications.maxquant_mod_to_proforma(fixed_mods=fixed_mods_tmp)(
            psm.peptide
        )

    perc_result = results_dict[psm.raw_file].get((psm.scannr, peptide), None)
    if (
        pout_input_type == "prosit"
        and not perc_result
        and maxquant.has_unknown_silac_label(psm.labeling_state)
    ):
        fixed_mods_tmp = modifications.SILAC_HEAVY_FIXED_MODS
        peptide = modifications.maxquant_mod_to_proforma(fixed_mods=fixed_mods_tmp)(
            psm.peptide
        )
        perc_result = results_dict[psm.raw_file].get((psm.scannr, peptide), None)
    return perc_result, peptide


def get_percolator_results(
    pout_files: List[str], pout_input_type: str
) -> Tuple[Dict[str, str], percolator.ResultsDict]:
    results_dict = collections.defaultdict(dict)
    fixed_mods = None
    for pout_file in pout_files:
        logger.info(f"Processing {pout_file}")
        fixed_mods, results_dict = percolator.parse_percolator_out_file_to_dict(
            pout_file, results_dict, pout_input_type
        )
        logger.debug(f"Found fixed modifications: {fixed_mods}")

    logger.info("Finished parsing percolator output files")
    logger.info("#PSMs per raw file:")
    for raw_file, psms in results_dict.items():
        logger.info(f"  {raw_file}: {len(psms)} PSMs")
    return fixed_mods, results_dict


def convert_PSM_dict_to_peptide_dict(
    results_dict: percolator.ResultsDict,
) -> percolator.ResultsDict:
    peptide_results_dict = dict()
    for _, results in results_dict.items():
        for (_, modified_peptide), (score, PEP) in results.items():
            peptide = helpers.remove_modifications(modified_peptide)
            curr_score, curr_PEP = peptide_results_dict.get(peptide, (-1e10, 1e10))
            # raw file = "" and scannr = -1 are hard coded values when parsing the peptides.txt file
            peptide_results_dict[""][(-1, peptide)] = (
                max(curr_score, score),
                min(curr_PEP, PEP),
            )
    return peptide_results_dict


def is_mbr_evidence_row(psm: maxquant.EvidenceRow) -> bool:
    return psm.scannr == -1


def is_valid_prosit_peptide(peptide: str) -> bool:
    return len(helpers.remove_modifications(peptide)) <= 30 and not "(ac)" in peptide


if __name__ == "__main__":
    main(sys.argv[1:])
