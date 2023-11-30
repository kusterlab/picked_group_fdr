import _csv
import sys
import os
import collections
import logging
from typing import Dict, List, Tuple

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
        "--mq_msms",
        default=None,
        metavar="M",
        nargs="+",
        help="""MaxQuant msms.txt file(s) (optional). Can help 
                resolve MS1 features with multiple scans. If you
                want to combine multiple evidence files, use 
                spaces to separate the file paths. Not 
                implemented yet...""",
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

    update_func = update_evidence_file
    if args.mq_input_type == "peptides":
        update_func = update_peptides_file

    update_func(
        args.mq_evidence,
        args.perc_results,
        args.mq_evidence_out,
        args.mq_msms,
        args.pout_input_type,
    )

    os.rename(args.mq_evidence_out + ".tmp", args.mq_evidence_out)


def update_evidence_file(
    evidence_files: List[str],
    pout_files: List[str],
    out_evidence_file: str,
    msms_files: List[str],
    pout_input_type: str,
) -> None:
    fixed_mods, results_dict = get_percolator_results(pout_files, pout_input_type)

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
            pout_input_type,
        )

    logger.info(f"Results written to {out_evidence_file}")


def update_evidence_single(
    evidence_file: str,
    writer: "_csv._writer",
    first_headers: List[str],
    fixed_mods: Dict[str, str],
    results_dict: Dict[str, Dict[Tuple[int, str], Tuple[float, float]]],
    pout_input_type: str,
) -> List[str]:
    logger.info(f"Processing {evidence_file}")
    reader = tsv.get_tsv_reader(evidence_file)

    headers_original = next(reader)  # save the header
    headers = list(map(lambda x: x.lower(), headers_original))
    if len(first_headers) == 0:
        writer.writerow(headers_original)
    elif headers != first_headers:
        warn_for_header_difference(first_headers, headers)

    # these columns will be updated
    score_col = tsv.get_column_index(headers, "score")
    post_err_prob_col = tsv.get_column_index(headers, "pep")

    mq_PEPs, prosit_PEPs = list(), list()
    unexplained_missing_PSMs = 0
    unexplained_peptides = list()
    rows_written = 0
    missing_raw_files = set()
    for row, psm in maxquant.parse_evidence_file_for_percolator_matching(
        reader, headers
    ):
        if rows_written % 500000 == 0:
            logger.info(f"    Writing line {rows_written}")

        if len(results_dict) == 0 or is_mbr_evidence_row(psm):
            rows_written += 1
            writer.writerow(row)
            continue

        if len(results_dict[psm.raw_file]) == 0:
            if psm.raw_file not in missing_raw_files:
                logger.warning(
                    f"Found no PSMs for {psm.raw_file} in percolator result files"
                )
                missing_raw_files.add(psm.raw_file)
            continue

        if not psm.is_decoy:
            mq_PEPs.append((psm.post_err_prob, psm.id_type))

        perc_result, peptide = find_percolator_psm(
            psm, fixed_mods, results_dict, pout_input_type
        )
        if not perc_result:
            if is_unexplainable_missing_psm(psm, peptide, pout_input_type):
                unexplained_missing_PSMs += 1
                logger.debug("Unexplained missing PSM:")
                if unexplained_missing_PSMs <= 10:
                    unexplained_peptides.append(psm.peptide)
            logger.debug(
                f"Missing PSM in percolator output: {psm.raw_file}, {peptide}, {psm.scannr}"
            )
            continue

        perc_score, perc_post_err_prob = perc_result
        if not psm.is_decoy:
            prosit_PEPs.append((perc_post_err_prob, psm.id_type))

        row[score_col] = perc_score
        row[post_err_prob_col] = perc_post_err_prob

        rows_written += 1
        writer.writerow(row)

    unexplained_percentage = int(unexplained_missing_PSMs / rows_written * 100)
    logger.info(
        f"Unexplained missing PSMs in Percolator results: {unexplained_missing_PSMs} out of {rows_written} ({unexplained_percentage}%)"
    )
    if unexplained_missing_PSMs > 0:
        logger.info(
            "If this percentage is low (<5%), it is probably the result of second peptide hits."
        )
        logger.info("\tFirst 10 missing peptides:")
        for peptide in unexplained_peptides:
            logger.info("\t" + peptide)

    logger.info("#MQ identifications:")
    count_below_FDR(mq_PEPs)

    logger.info(f"#{pout_input_type} identifications:")
    count_below_FDR(prosit_PEPs)

    return headers


def warn_for_header_difference(first_headers, headers):
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


def is_unexplainable_missing_psm(psm, peptide, pout_input_type):
    return (
        not psm.is_contaminant
        and psm.score > 0.0
        and not (pout_input_type == "prosit" and not is_valid_prosit_peptide(peptide))
    )


def find_percolator_psm(
    psm: maxquant.EvidenceRow,
    fixed_mods: Dict[str, str],
    results_dict: Dict[str, Dict[Tuple[int, str], Tuple[float, float]]],
    pout_input_type: str,
) -> Tuple[Tuple[float, float], str]:
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
        peptide = modifications.maxquant_mod_to_unimod_single(
            psm.peptide, fixed_mods=fixed_mods_tmp
        )

    perc_result = results_dict[psm.raw_file].get((psm.scannr, peptide), None)
    if (
        pout_input_type == "prosit"
        and not perc_result
        and maxquant.has_unknown_silac_label(psm.labeling_state)
    ):
        fixed_mods_tmp = modifications.SILAC_HEAVY_FIXED_MODS
        peptide = modifications.maxquant_mod_to_unimod_single(
            psm.peptide, fixed_mods=fixed_mods_tmp
        )
        perc_result = results_dict[psm.raw_file].get((psm.scannr, peptide), None)
    return perc_result, peptide


def update_peptides_file(
    peptide_files: List[str],
    pout_files: List[str],
    out_peptide_file: str,
    msms_files: List[str],
    pout_input_type: str,
) -> None:
    _, results_dict = get_percolator_results(pout_files, pout_input_type)

    peptide_results_dict = convert_PSM_dict_to_peptide_dict(results_dict)

    logger.info("Writing updated peptides file")
    writer = tsv.get_tsv_writer(out_peptide_file + ".tmp")
    first = True
    mq_PEPs, prosit_PEPs = list(), list()
    unexplained_missing_peptides = 0
    unexplained_peptides = list()
    for peptide_file in peptide_files:
        logger.info(f"Processing {peptide_file}")
        reader = tsv.get_tsv_reader(peptide_file)
        headers_original = next(reader)  # save the header
        headers = list(map(lambda x: x.lower(), headers_original))

        if first:
            first = False
            writer.writerow(headers_original)
            first_headers = headers
        elif headers != first_headers:
            logger.warning(
                "Current column names are different from the first evidence file"
            )
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

        score_col = tsv.get_column_index(headers, "score")
        post_err_prob_col = tsv.get_column_index(headers, "pep")

        peptide_col = tsv.get_column_index(headers, "sequence")
        reverse_col = tsv.get_column_index(headers, "reverse")

        for row in reader:
            if len(pout_files) > 0:
                peptide = row[peptide_col]
                is_decoy = row[reverse_col] == "+"

                perc_result = peptide_results_dict.get(peptide, None)

                if not is_decoy:
                    mq_PEPs.append((float(row[post_err_prob_col]), "Unknown"))
                    if not perc_result and (
                        pout_input_type != "prosit" or is_valid_prosit_peptide(peptide)
                    ):  # mysterious missing PSMs
                        unexplained_missing_peptides += 1
                        if unexplained_missing_peptides <= 10:
                            unexplained_peptides.append(row[peptide_col])

                if perc_result:
                    if not is_decoy:
                        prosit_PEPs.append((perc_result[1], "Unknown"))
                    row[score_col] = perc_result[0]
                    row[post_err_prob_col] = perc_result[1]

            writer.writerow(row)

    logger.info(
        f"Unexplained missing peptides in Percolator results: {unexplained_missing_peptides}"
    )
    if unexplained_missing_peptides > 0:
        logger.info("\tFirst 10 missing peptides:")
        for peptide in unexplained_peptides:
            logger.info("\t" + peptide)
    logger.info(f"Results written to {out_peptide_file}")

    logger.info("MQ")
    count_below_FDR(mq_PEPs)

    logger.info(pout_input_type)
    count_below_FDR(prosit_PEPs)


def get_percolator_results(
    pout_files: List[str], pout_input_type: str
) -> Tuple[Dict[str, str], Dict[str, Dict[Tuple[int, str], Tuple[float, float]]]]:
    results_dict = collections.defaultdict(dict)
    fixed_mods = None
    for pout_file in pout_files:
        logger.info(f"Processing {pout_file}")
        fixed_mods, results_dict = percolator.parse_percolator_out_file_to_dict(
            pout_file, results_dict, pout_input_type
        )

    logger.info("Finished parsing percolator output files")
    logger.info("#PSMs per raw file:")
    for raw_file, psms in results_dict.items():
        logger.info(f"  {raw_file}: {len(psms)} PSMs")
    return fixed_mods, results_dict


def convert_PSM_dict_to_peptide_dict(
    results_dict: Dict[str, Dict[Tuple[int, str], Tuple[float, float]]],
) -> Dict[str, Tuple[float, float]]:
    peptide_results_dict = dict()
    for _, results in results_dict.items():
        for (_, peptide), (score, PEP) in results.items():
            peptide = helpers.clean_peptide(peptide, remove_flanks=False)
            curr_score, curr_PEP = peptide_results_dict.get(peptide, (-1e10, 1e10))
            peptide_results_dict[peptide] = (max(curr_score, score), min(curr_PEP, PEP))
    return peptide_results_dict


def is_mbr_evidence_row(psm: maxquant.EvidenceRow) -> bool:
    return psm.scannr == -1


def is_valid_prosit_peptide(peptide: str) -> bool:
    return (
        len(helpers.clean_peptide(peptide, remove_flanks=False)) <= 30
        and not "(ac)" in peptide
    )


def count_below_FDR(post_err_probs: List[float], fdr_threshold: float = 0.01):
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


if __name__ == "__main__":
    main(sys.argv[1:])
