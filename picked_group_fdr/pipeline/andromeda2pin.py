#!/usr/bin/python

"""
Converts Andromeda evidence.txt output to tab delimited percolator input file
"""

import sys
import os
import logging
from typing import List

import numpy as np

from ..parsers import tsv
from .. import __version__, __copyright__
from .. import digest
from .. import helpers
from ..picked_group_fdr import ArgumentParserWithLogger
from ..digestion_params import get_digestion_params_list, add_digestion_arguments


# hacky way to get package logger when running as module
logger = logging.getLogger(__package__ + "." + __file__)


# TODO allow mqpar.xml as input
def main(argv):
    logger.info(f"Andromeda2Pin version {__version__}\n{__copyright__}")
    logger.info(
        f'Issued command: {os.path.basename(__file__)} {" ".join(map(str, argv))}'
    )

    args = parseArgs(argv)

    andromedaTargetOutFNs = getMqEvidenceFiles(args.mq_evidence_file)

    percInFN = args.outputTab
    decoyPattern = args.pattern
    # numHits = args.matches
    numHits = 1

    if len(percInFN) > 0:
        if os.path.isfile(percInFN):
            logger.info(
                f"Found output file {percInFN}, remove this file to re-run andromeda2pin."
            )
            return
        logger.info(f"Writing results to: {percInFN}")
        writer = tsv.get_tsv_writer(percInFN + ".tmp")
    else:
        logger.info("Writing results to stdout")
        writer = tsv.get_tsv_writer(sys.stdout)

    digestion_params_list = get_digestion_params_list(args)

    charges = list(range(2, 7))
    writeHeaders(writer, charges)

    for andromedaTargetOutFN, digestion_params in zip(
        andromedaTargetOutFNs, digestion_params_list
    ):
        peptideToProteinMap = digest.get_peptide_to_protein_map_from_params(
            args.databases, [digestion_params]
        )
        convertAndromedaOutToPin(
            andromedaTargetOutFN,
            writer,
            charges,
            numHits,
            peptideToProteinMap,
            args.suppress_missing_peptide_warning,
            decoyPattern=decoyPattern,
        )

    if len(percInFN) > 0:
        os.rename(percInFN + ".tmp", percInFN)
    logger.info("Finished writing percolator input")


def parseArgs(argv):
    import argparse

    apars = ArgumentParserWithLogger(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    apars.add_argument(
        "mq_evidence_file",
        default=None,
        nargs="+",
        metavar="evidence.txt",
        help="""MaxQuant evidence file(s), or a meta file. If you want to combine 
             multiple evidence files, use spaces to separate the file paths or
             use a meta file. Meta files are text files containing the paths of 
             evidence files, one path per line.
             """,
    )

    apars.add_argument(
        "-o",
        "--outputTab",
        default=None,
        metavar="pin.tab",
        help="""Save output in a tab delimited file""",
    )

    # apars.add_argument(
    #     "-m",
    #     "--matches",
    #     default=1,
    #     metavar="M",
    #     type=int,
    #     help="""Maximal number of matches to take in consideration per spectrum.""",
    # )

    apars.add_argument(
        "-P",
        "--pattern",
        default="REV__",
        metavar="P",
        help="""Pattern used to identify the decoy PSMs.""",
    )

    apars.add_argument(
        "-F",
        "--databases",
        default=None,
        metavar="F",
        nargs="+",
        help="""Fasta database used in the search against the spectra file.""",
    )

    apars.add_argument(
        "--suppress_missing_peptide_warning",
        help="Suppress missing peptide warning when mapping peptides to proteins.",
        action="store_true",
    )

    add_digestion_arguments(apars)

    # ------------------------------------------------
    args = apars.parse_args(argv)

    return args


def isMqEvidenceHeader(header: str, delimiter: str):
    """Check if header looks like it comes from an MaxQuant evidence file.

    If there are at least 5 columns and contains the word "sequence", we can safely
    assume this is a MQ evidence file and not a metafile.
    """
    return header.count(delimiter) > 4 and "sequence" in header.lower()


def getMqEvidenceFiles(mq_evidence_files: List[str]):
    andromedaTargetOutFNs = []
    delimiter = tsv.get_delimiter(mq_evidence_files[0])
    with open(mq_evidence_files[0], "r") as f:
        firstLine = True
        for line in f:
            if firstLine:
                if isMqEvidenceHeader(line, delimiter):
                    andromedaTargetOutFNs = mq_evidence_files
                    break
                else:
                    logger.info(
                        "Meta file detected, interpreting each line as a path to an evidence file"
                    )
            andromedaTargetOutFNs.append(line.rstrip())
            firstLine = False
    return andromedaTargetOutFNs


def writeHeaders(writer, charges):
    writer.writerow(
        [
            "SpecId",
            "Label",
            "FileName",
            "ScanNr",
            "ExpMass",
            "AndromedaScore",
            "DeltaScore",
            "PepLen",
        ]
        + ["Charge" + str(i) for i in charges]
        + [
            "Mass",
            "enzN",
            "enzC",
            "enzInt",
            "numMods",
            "dM",
            "absdM",
            "Peptide",
            "Proteins",
        ]
    )
    # writer.writerow(["DefaultDirection", "-", "-", "-", "-", 1, 0.5, 0] + [0 for i in charges] + [0, 0, 0, -1.5, -2, 0, -1])


def parseMqEvidenceFile(mqEvidenceFile, razor=False):
    """
    Columns needed (all headers are converted to lower case for the comparison, so no need to fix that):
    - Modified sequence
    - MS/MS Scan number
    - Raw file
    - Charge
    - Mass
    - Mass error [ppm] (optional)
    - One out of: Leading proteins / Proteins
    - Score
    - Delta score
    - Experiment (optional)
    """
    delimiter = tsv.get_delimiter(mqEvidenceFile)
    reader = tsv.get_tsv_reader(mqEvidenceFile, delimiter)
    headers = next(reader)  # save the header
    headers = list(map(lambda x: x.lower(), headers))

    # People often provide custom evidence files generated by R which replaces spaces
    # with periods, undo that here
    if mqEvidenceFile.endswith(".csv"):
        headers = [x.replace(".", " ") for x in headers]

    peptCol = headers.index("modified sequence")
    idCol = headers.index("ms/ms scan number")
    fileCol = headers.index("raw file")
    chargeCol = headers.index("charge")
    massCol = headers.index("mass")

    deltaMassCol = -1
    if "mass error [ppm]" in headers:
        deltaMassCol = headers.index("mass error [ppm]")

    if "leading proteins" in headers:
        proteinCol = headers.index(
            "leading proteins"
        )  # all proteins the peptide matches to, does not return REV__ proteins
    else:
        proteinCol = headers.index("proteins")
    scoreCol = headers.index("score")
    deltaScoreCol = headers.index("delta score")

    experimentCol = -1
    if "experiment" in headers:
        experimentCol = headers.index("experiment")

    logger.info("Parsing MaxQuant evidence.txt file")
    for lineIdx, row in enumerate(reader):
        if lineIdx % 500000 == 0:
            logger.info(f"    Reading line {lineIdx}")

        if len(row[idCol]) == 0:
            continue

        scanNr = int(row[idCol])
        charge = int(row[chargeCol])
        fileName = row[fileCol]
        peptide = (
            "-."
            + row[peptCol][1:-1]
            .replace("pS", "S[80]")
            .replace("pT", "T[80]")
            .replace("pY", "Y[80]")
            .replace("(ph)", "[80]")
            .replace("(ox)", "[16]")
            .replace("(ac)", "[42]")
            .replace("(Acetyl (Protein N-term))", "[42]")
            .replace("(Oxidation (M))", "[16]")
            + ".-"
        )
        proteins = row[proteinCol].split(";")
        experiment = "Experiment1"
        if experimentCol >= 0:
            experiment = row[experimentCol]            

        score = float(row[scoreCol])
        deltaScore = float(row[deltaScoreCol])
        mass = float(row[massCol])
        deltaMass = 0.0
        if deltaMassCol >= 0:
            deltaMass = float(row[deltaMassCol])
            if np.isnan(deltaMass):
                deltaMass = 0.0

        if not np.isnan(score) and score > 0.0:
            yield scanNr, charge, fileName, peptide, proteins, experiment, score, deltaScore, mass, deltaMass


def convertAndromedaOutToPin(
    andromedaOutFN,
    writer,
    charges,
    numHits,
    peptideToProteinMap,
    suppress_missing_peptide_warning,
    decoyPattern="",
):
    logger.info(f"Reading {andromedaOutFN}")

    for (
        scanNr,
        charge,
        fileName,
        peptide,
        tmp_proteins,
        experiment,
        andromedaScore,
        deltaScore,
        expMass,
        deltaMass,
    ) in parseMqEvidenceFile(andromedaOutFN):
        rank = 1
        psmId = fileName + "_" + str(scanNr) + "_" + str(charge) + "_" + str(rank)
        modPeptide, cleanPeptide, pepLen, enzN, enzC, enzInt, numMods = getPeptideStats(
            peptide, deltaMass
        )
        absDeltaMass = abs(deltaMass)

        if pepLen < 6:
            continue

        if len(peptideToProteinMap) > 0:
            proteins = digest.get_proteins(peptideToProteinMap, cleanPeptide[2:-2])
            if len(proteins) == 0:
                if (
                    not helpers.is_contaminant(tmp_proteins)
                    and not suppress_missing_peptide_warning
                ):
                    logger.warning(
                        f"Could not find peptide {peptide} ({str(tmp_proteins)}) in fasta database, skipping PSM"
                    )
                continue

        if len(decoyPattern) > 0:
            if sum(1 for p in proteins if p.startswith(decoyPattern)) == len(proteins):
                label = -1
            else:
                label = 1

        r = (
            [
                psmId,
                label,
                fileName,
                scanNr,
                expMass,
                andromedaScore,
                deltaScore,
                pepLen,
            ]
            + [1 if charge == i else 0 for i in charges]
            + [
                expMass,
                enzN,
                enzC,
                enzInt,
                numMods,
                deltaMass,
                absDeltaMass,
                modPeptide,
            ]
            + proteins
        )
        writer.writerow(r)


def getPeptideStats(peptide, deltaMass, accurateModMasses=False):
    cleanPeptide = peptide[:2]
    modPeptide = peptide[:2]
    numMods, enzInt = 0, 0
    aaIdx = 2
    while aaIdx < len(peptide):
        if peptide[aaIdx] == "[":
            modStart = aaIdx
            while peptide[aaIdx + 1].isdigit():
                aaIdx += 1
            mod = "["

            multFactor = 1
            if peptide[modStart] == "-":
                multFactor = -1
            modMass = multFactor * float(peptide[modStart + 1 : aaIdx + 1])

            # modification masses are typically not integers, add the deltamass to the first mod to ensure calcmass = expmass
            if accurateModMasses and numMods == 0:
                modMass += deltaMass
                mod += str(round(modMass, 4))
            else:
                mod += str(int(modMass))

            mod += "]"
            aaIdx += 1
            modPeptide += mod
            numMods += 1
        elif peptide[aaIdx] == ".":
            modPeptide += peptide[aaIdx:]
            cleanPeptide += peptide[aaIdx:]
            break
        else:
            modPeptide += peptide[aaIdx]
            cleanPeptide += peptide[aaIdx]
            if digest.is_enzymatic_advanced(
                cleanPeptide[-2],
                cleanPeptide[-1],
                not_post=[],
                methionine_cleavage=False,
            ):
                enzInt += 1
        aaIdx += 1
    enzN = int(
        digest.is_enzymatic_advanced(cleanPeptide[0], cleanPeptide[2], not_post=[])
    )  # Andromeda uses Trypsin/P
    enzC = int(
        digest.is_enzymatic_advanced(cleanPeptide[-3], cleanPeptide[-1], not_post=[])
    )
    pepLen = len(cleanPeptide) - 4
    return modPeptide, cleanPeptide, pepLen, enzN, enzC, enzInt, numMods


if __name__ == "__main__":
    main(sys.argv[1:])
