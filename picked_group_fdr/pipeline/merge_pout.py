import collections
import sys
import os
import logging

import numpy as np
import triqler.qvality as qvality

from .. import __version__, __copyright__
from .. import helpers
from .. import digest
from ..picked_group_fdr import ArgumentParserWithLogger
from ..parsers import tsv, percolator
from ..digestion_params import add_digestion_arguments

# hacky way to get package logger when running as module
logger = logging.getLogger(__package__ + "." + __file__)


def parseArgs(argv):
    import argparse

    apars = ArgumentParserWithLogger(
        description="""Merges percolator output files to a list of ranked 
                        peptides, again in percolator output format.
                        If a fasta file is provided, peptides are remapped 
                        to their protein sequences based on an in-silico 
                        digestion of the fasta database.
                        Peptides are then ranked based on their q-values.
                        For target peptides, q-values are re-estimated from
                        the posterior error probabilities to obtain a
                        higher resolution in the high-confident region.""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    apars.add_argument(
        "--perc_results",
        default=list(),
        metavar="POUT",
        nargs="+",
        help="""Percolator output file(s), separated by 
                spaces. Also include the decoy results!""",
    )

    apars.add_argument(
        "--fasta",
        default=None,
        metavar="F",
        required=True,
        nargs="+",
        help="""Fasta file(s) with protein sequences""",
    )

    apars.add_argument(
        "--perc_merged",
        default="merged_pout.txt",
        metavar="M",
        help="""Output file with merged pout results.""",
    )

    add_digestion_arguments(apars)

    # ------------------------------------------------
    args = apars.parse_args(argv)

    return args


def main(argv):
    logger.info(f"MergePout version {__version__}\n{__copyright__}")
    logger.info(
        f'Issued command: {os.path.basename(__file__)} {" ".join(map(str, argv))}'
    )

    args = parseArgs(argv)

    if os.path.isfile(args.perc_merged):
        logger.info(
            f"Found merged pout file {args.perc_merged}, remove this file to rerun merge_pout."
        )
        return

    peptide_to_protein_map = collections.defaultdict(list)
    if len(args.fasta) > 0:
        digestion_params_list = digest.get_digestion_params_list(args)
        peptide_to_protein_map = digest.get_peptide_to_protein_map_from_params(
            args.fasta, digestion_params_list
        )

    merge_pout(
        args.perc_results,
        peptide_to_protein_map,
        args.perc_merged,
        args.suppress_missing_peptide_warning,
    )


def merge_pout(
    perc_results, peptideToProteinMap, perc_merged, suppress_missing_peptide_warning
):
    seenPeptides = dict()
    missingPeptides, matchedPeptides = 0, 0
    for poutFile in perc_results:
        poutReader = tsv.get_tsv_reader(poutFile)
        headers = next(poutReader)

        (
            idCol,
            _,
            peptCol,
            scoreCol,
            qvalCol,
            postErrProbCol,
            proteinCol,
        ) = percolator.get_percolator_column_idxs(headers)

        sumPEP = 0.0
        matchedBefore = matchedPeptides
        logger.info(f"Parsing {poutFile}")
        # TODO: merge with parsers.parsePercolatorOutFile()
        for i, row in enumerate(poutReader):
            if i % 1000000 == 0:
                logger.info(
                    f"Processing row {i}: #Missing peptides: {missingPeptides}, #Matched peptides: {matchedPeptides}"
                )

            # convert peptide string to upper case, since prosit converts modified amino acids to lower case
            peptide = helpers.remove_modifications(row[peptCol][2:-2].upper())

            if len(peptideToProteinMap) > 0:
                proteins = digest.get_proteins(peptideToProteinMap, peptide)
            elif percolator.is_native_percolator_file(headers):
                proteins = row[proteinCol:]
            elif percolator.is_mokapot_file(headers):
                proteins = row[proteinCol].split("\t")

            isDecoy = helpers.is_decoy(proteins)

            if isDecoy:
                qValue = float(row[qvalCol])
            else:
                # for targets, estimate "high resolution" q-value based on PEPs
                sumPEP += float(row[postErrProbCol])
                qValue = sumPEP / (i + 1)

            if qValue < seenPeptides.get(peptide, (1.0, []))[0]:
                if len(proteins) > 0:
                    row = [
                        row[idCol],
                        row[scoreCol],
                        row[qvalCol],
                        row[postErrProbCol],
                        "-." + peptide + ".-",
                    ] + proteins
                    matchedPeptides += 1
                    seenPeptides[peptide] = (qValue, row, isDecoy)
                else:
                    if (
                        not helpers.is_contaminant(proteins)
                        and not suppress_missing_peptide_warning
                    ):
                        logger.debug(
                            f"Could not find peptide {peptide} in fasta file, check your database and if the correct digestion parameters were specified"
                        )
                    missingPeptides += 1

        logger.info(
            f"Processing row {i}: #Missing peptides: {missingPeptides}, #Matched peptides: {matchedPeptides}"
        )

        if matchedBefore == matchedPeptides:
            logger.warning(f"No new peptides added by {poutFile}")

    psm_infos = sorted(seenPeptides.values())
    peps = get_peptide_PEPs(psm_infos)

    write_updated_PSMs(perc_merged, psm_infos, peps, update_qvals=False)


def get_peptide_PEPs(psm_infos):
    targetScores, decoyScores = list(), list()
    for qval, _, isDecoy in psm_infos:
        if isDecoy:
            decoyScores.append(-1 * np.log10(qval))
        else:
            targetScores.append(-1 * np.log10(qval))

    _, peps = qvality.getQvaluesFromScores(
        targetScores, decoyScores, includeDecoys=True, plotRegressionCurve=False
    )
    return peps


def write_updated_PSMs(perc_merged, psm_infos, peps, update_qvals=False):
    writer = tsv.get_tsv_writer(perc_merged + ".tmp")

    headers = percolator.PERCOLATOR_NATIVE_HEADERS
    writer.writerow(headers)
    _, _, _, _, qvalCol, _, _ = percolator.get_percolator_column_idxs(headers)

    sumPEP = 0.0
    counts = 0
    for (_, row, isDecoy), pep in zip(psm_infos, peps):
        if update_qvals:
            if not isDecoy:
                sumPEP += pep
                counts += 1
            row[qvalCol] = sumPEP / counts
        writer.writerow(row)

    os.rename(perc_merged + ".tmp", perc_merged)


if __name__ == "__main__":
    main(sys.argv[1:])
