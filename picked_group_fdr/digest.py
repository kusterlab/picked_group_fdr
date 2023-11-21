from __future__ import print_function, annotations
from pathlib import Path

import sys
import csv
import itertools
import collections
import logging
from typing import Dict, Iterator, List

from .digestion_params import (
    DigestionParams,
    add_digestion_arguments,
    get_digestion_params_list,
)


logger = logging.getLogger(__name__)

ENZYME_CLEAVAGE_RULES = {
    "trypsin": {"pre": ["K", "R"], "not_post": ["P"], "post": []},
    "trypsinp": {"pre": ["K", "R"], "not_post": [], "post": []},
    "no_enzyme": {"pre": [], "not_post": [], "post": []},
    "chymotrypsin": {"pre": ["F", "W", "Y", "L"], "not_post": ["P"], "post": []},
    "chymotrypsin+": {"pre": ["F", "W", "Y", "L", "M"], "not_post": [], "post": []},
    "proteinasek": {
        "pre": ["A", "E", "F", "I", "L", "T", "V", "W", "Y"],
        "not_post": [],
        "post": [],
    },
    "elastase": {"pre": ["L", "V", "A", "G"], "not_post": ["P"], "post": []},
    "clostripain": {"pre": ["R"], "not_post": [""], "post": []},
    "cyanogen-bromide": {"pre": ["M"], "not_post": [""], "post": []},
    "iodosobenzoate": {"pre": ["W"], "not_post": [""], "post": []},
    "proline-endopeptidase": {"pre": ["P"], "not_post": [""], "post": []},
    "staph-protease": {"pre": ["E"], "not_post": [""], "post": []},
    "asp-n": {"pre": [""], "not_post": [""], "post": ["D"]},
    "lys-c": {"pre": ["K"], "not_post": ["P"], "post": []},
    "lys-cp": {"pre": ["K"], "not_post": [""], "post": []},
    "lys-n": {"pre": [""], "not_post": [""], "post": ["K"]},
    "arg-c": {"pre": ["R"], "not_post": ["P"], "post": []},
    "glu-c": {"pre": ["E"], "not_post": ["P"], "post": []},
    "pepsin-a": {"pre": ["F", "L"], "not_post": ["P"], "post": []},
    "elastase-trypsin-chymotrypsin": {
        "pre": ["A", "L", "I", "V", "F", "K", "R", "W", "F", "Y"],
        "not_post": ["P"],
        "post": [],
    },
    "lysarginase": {"pre": [""], "not_post": [""], "post": ["K", "R"]},
    "v8-de": {"pre": ["N", "D", "E", "Q"], "not_post": ["P"], "post": []},
}


def main(argv):
    args = parseArgs()

    digestion_params_list = get_digestion_params_list(args)

    if args.prosit_input:
        writer = getTsvWriter(args.prosit_input, delimiter=",")
        writer.writerow(
            "modified_sequence,collision_energy,precursor_charge".split(",")
        )

        prositInputFileWithProteins = args.prosit_input.replace(
            ".csv", "_with_proteins.csv"
        )
        writerWithProteins = getTsvWriter(prositInputFileWithProteins, delimiter=",")
        writerWithProteins.writerow(
            "modified_sequence,collision_energy,precursor_charge,protein".split(",")
        )

        for peptide, proteins in get_peptide_to_protein_map_from_params(
            args.fasta, digestion_params_list
        ).items():
            if not validPrositPeptide(peptide):
                continue

            for charge in [2, 3, 4]:
                writer.writerow([peptide, 30, charge])
                writerWithProteins.writerow([peptide, 30, charge, proteins[0]])

    if args.peptide_protein_map:
        with open(args.peptide_protein_map + ".params.txt", "w") as f:
            f.write(" ".join(sys.argv))

        writer = getTsvWriter(args.peptide_protein_map, delimiter="\t")
        for peptide, proteins in get_peptide_to_protein_map_from_params(
            args.fasta, digestion_params_list
        ).items():
            writer.writerow([peptide, ";".join(proteins)])

    if args.ibaq_map:
        writer = getTsvWriter(args.ibaq_map, delimiter="\t")

        numPeptidesPerProtein = getNumIbaqPeptidesPerProtein(
            args.fasta, digestion_params_list
        )
        for protein, numPeptides in numPeptidesPerProtein.items():
            writer.writerow([protein, numPeptides])

    # writeProteinToGeneMap(fastaFile, outputFile)


def validPrositPeptide(peptide):
    return len(peptide) <= 30 and "U" not in peptide and "X" not in peptide


def parseArgs():
    import argparse

    apars = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    apars.add_argument(
        "--fasta",
        default=None,
        metavar="F",
        required=True,
        nargs="+",
        help="""Fasta file used as input
                                                    """,
    )

    apars.add_argument(
        "--prosit_input",
        default=None,
        metavar="M",
        required=False,
        help="""Path to file where to write the prosit input file.
                                                    """,
    )

    apars.add_argument(
        "--peptide_protein_map",
        default=None,
        metavar="M",
        required=False,
        help="""Write mapping from peptides to all its proteins to 
                                                         the specified file.
                                                    """,
    )

    apars.add_argument(
        "--ibaq_map",
        default=None,
        metavar="M",
        required=False,
        help="""Write number of peptides per protein to the specified 
                                                         file that meet the iBAQ criteria
                                                         (6 <= pepLen <= 30, no miscleavages).
                                                    """,
    )

    add_digestion_arguments(apars)

    # ------------------------------------------------
    args = apars.parse_args()

    return args


def writeProteinToGeneMap(fastaFile, outputFile):
    writer = csv.writer(open(outputFile, "w"), delimiter="\t")
    for proteinName, _ in readFastaTide(fastaFile, db="target"):
        proteinId = proteinName.split("|")[1]
        geneId = proteinName.split("|")[2].split(" ")[0]
        writer.writerow([proteinId, geneId])


def parseUntilFirstSpace(fastaId: str) -> str:
    return fastaId.split(" ")[0]


def readFastaTide(filePath, db="target", parseId=parseUntilFirstSpace):
    readFastaMaxQuant(filePath, db, parseId, specialAAs=[], decoyPrefix="decoy_")


def readFastaMaxQuant(
    filePath,
    db="target",
    parseId=parseUntilFirstSpace,
    specialAAs=["K", "R"],
    decoyPrefix="REV__",
):
    if db not in ["target", "decoy", "concat"]:
        sys.exit("unknown db mode: %s" % db)

    hasSpecialAAs = len(specialAAs) > 0
    name, seq = None, []
    with open(filePath, "r") as fp:
        for line in itertools.chain(fp, [">"]):
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    seq = "".join(seq)
                    if db in ["target", "concat"]:
                        yield (name, seq)

                    if db in ["decoy", "concat"]:
                        revSeq = seq[::-1]
                        if hasSpecialAAs:
                            revSeq = swapSpecialAAs(revSeq, specialAAs)
                        yield (decoyPrefix + name, revSeq)

                if len(line) > 1:
                    name, seq = parseId(line[1:]), []
            else:
                seq.append(line)


# from . import digestfast
# readFasta = digestfast.readFastaMaxQuant
readFasta = readFastaMaxQuant


# swaps the specialAAs with its preceding amino acid, as is done in MaxQuant
# e.g. specialAAs = ['R', 'K'] transforms ABCKDEFRK into ABKCDERKF
def swapSpecialAAs(seq, specialAAs):
    seq = list(seq)
    for i in range(1, len(seq)):
        if seq[i] in specialAAs:
            swapPositions(seq, i, i - 1)
    seq = "".join(seq)
    return seq


def swapPositions(seq, pos1, pos2):
    seq[pos1], seq[pos2] = seq[pos2], seq[pos1]


def getProteinIds(filePath):
    proteinIds = list()
    for proteinId, _ in readFasta(filePath):
        proteinIds.append(proteinId)
    return set(proteinIds)


def getProteinSequences(filePaths, parseId):
    proteinSequences = dict()
    for filePath in filePaths:
        for proteinId, proteinSeq in readFasta(filePath, db="concat", parseId=parseId):
            if (
                proteinId not in proteinSequences
            ):  # keep only first sequence per identifier
                proteinSequences[proteinId] = proteinSeq
    return proteinSequences


def filterFastaFile(fastaFile, filteredFastaFile, proteins, **kwargs):
    with open(filteredFastaFile, "w") as f:
        for prot, seq in readFasta(fastaFile, **kwargs):
            if prot in proteins:
                f.write(">" + prot + "\n" + seq + "\n")
                # f.write('>decoy_' + prot + '\n' + seq[::-1] + '\n')


def getPeptides(
    fastaFile,
    db="concat",
    min_len=6,
    max_len=50,
    pre=["K", "R"],
    not_post=["P"],
    post=[],
    digestion="full",
    miscleavages=0,
    methionineCleavage=True,
):
    for protein, seq in readFasta(fastaFile, db):
        if len(seq) == 0:
            raise ValueError(
                f"Found an empty sequence for protein id {protein}, please check your fasta file."
            )
        for peptide in getDigestedPeptides(
            seq,
            min_len,
            max_len,
            pre,
            not_post,
            digestion,
            miscleavages,
            methionineCleavage,
        ):
            yield peptide


# @profile
def getDigestedPeptides(
    seq,
    min_len=6,
    max_len=50,
    pre=["K", "R"],
    not_post=["P"],
    post=[],
    digestion="full",
    miscleavages=0,
    methionineCleavage=True,
):
    if digestion == "none":
        yield from nonSpecificDigest(seq, min_len, max_len)
    elif digestion == "semi":
        yield from semiSpecificDigest(
            seq, min_len, max_len, pre, not_post, post, miscleavages, methionineCleavage
        )
    else:
        yield from fullDigest(
            seq, min_len, max_len, pre, not_post, post, miscleavages, methionineCleavage
        )


def nonSpecificDigest(seq, min_len, max_len):
    lenS = len(seq)
    for i in range(lenS + 1):
        for j in range(i + min_len, min(lenS + 1, i + max_len + 1)):
            if j <= lenS:
                yield seq[i:j]


def semiSpecificDigest(
    seq, min_len, max_len, pre, not_post, post, miscleavages, methionineCleavage
):
    lenS, starts = len(seq), [0]
    methionineCleavage = methionineCleavage and seq[0] == "M"
    length_accepted = lambda x: x >= min_len and x <= max_len

    for i in range(lenS + 1):
        isCleavageSite = is_enzymatic(
            seq[min([lenS - 1, i])], seq[min([lenS - 1, i + 1])], pre, not_post, post
        )
        isMethionineCleavageSite = i == 0 and methionineCleavage
        if i == lenS or isCleavageSite or isMethionineCleavageSite:
            # peptides with enzymatic C-terminal (both enzymatic and non-enzymatic N-terminal)
            start = starts[0]
            for j in range(start, min([i + 1, lenS])):
                lenP = min([i, lenS - 1]) - j + 1
                if length_accepted(lenP):
                    yield (seq[j : i + 1])
            starts.append(i + 1)
            methionineCleaved = int(starts[0] == 0 and methionineCleavage)
            if len(starts) > miscleavages + 1 + methionineCleaved or i == lenS:
                starts = starts[1 + methionineCleaved :]
        else:  # peptides with non enzymatic C-terminal
            for start in starts:
                lenP = i - start + 1
                if length_accepted(lenP) and i + 1 not in starts:
                    yield (seq[start : i + 1])


def fullDigest(
    seq, min_len, max_len, pre, not_post, post, miscleavages, methionineCleavage
):
    lenS, starts = len(seq), [0]
    methionineCleavage = methionineCleavage and seq[0] == "M"
    length_accepted = lambda x: x >= min_len and x <= max_len

    cleavageSites = [0] if methionineCleavage else []
    cleavageSites.extend(
        [
            i
            for i in range(lenS)
            if is_enzymatic(seq[i], seq[min([lenS - 1, i + 1])], pre, not_post, post)
        ]
    )
    cleavageSites.append(lenS)
    for i in cleavageSites:
        for start in starts:
            lenP = i - start + 1
            if length_accepted(lenP):
                yield (seq[start : i + 1])
        starts.append(i + 1)
        methionineCleaved = int(starts[0] == 0 and methionineCleavage)
        if len(starts) > miscleavages + 1 + methionineCleaved:
            starts = starts[1 + methionineCleaved :]


def getPeptideToProteinMapWithEnzyme(
    fastaFile, min_len, max_len, enzyme, miscleavages, specialAAs, db
):
    if len(fastaFile) == 0:
        return dict()

    pre, not_post, post = getCleavageSites(enzyme)
    return getPeptideToProteinMap(
        fastaFile,
        db,
        digestion="full",
        min_len=min_len,
        max_len=max_len,
        pre=pre,
        not_post=not_post,
        post=post,
        miscleavages=miscleavages,
        methionineCleavage=True,
        specialAAs=specialAAs,
    )


def get_peptide_to_protein_map_from_params(
    fasta_files: List[str], digestion_params_list: List[DigestionParams]
):
    peptideToProteinMap = collections.defaultdict(list)
    for fasta_file in fasta_files:
        for params in digestion_params_list:
            pre, not_post, post = getCleavageSites(params.enzyme)
            for peptide, proteins in getPeptideToProteinMap(
                fasta_file,
                params.db,
                digestion=params.digestion,
                min_len=params.min_length,
                max_len=params.max_length,
                pre=pre,
                not_post=not_post,
                post=post,
                miscleavages=params.cleavages,
                methionineCleavage=params.methionine_cleavage,
                specialAAs=params.special_aas,
            ).items():
                peptideToProteinMap[peptide].extend(proteins)
    return peptideToProteinMap


def merge_peptide_to_protein_maps(peptide_protein_maps: Iterator[Dict[str, List[str]]]):
    peptideToProteinMap = collections.defaultdict(list)
    for peptide_protein_map in peptide_protein_maps:
        for peptide, proteins in peptide_protein_map.items():
            peptideToProteinMap[peptide].extend(proteins)
    return peptideToProteinMap


def getPeptideToProteinMap(
    fastaFile,
    db="concat",
    min_len=6,
    max_len=52,
    pre=["K", "R"],
    not_post=["P"],
    post=[],
    digestion="full",
    miscleavages=2,
    methionineCleavage=True,
    useHashKey=False,
    specialAAs=["K", "R"],
    parseId=parseUntilFirstSpace,
):
    peptideToProteinMap = collections.defaultdict(list)
    proteinToSeqMap = dict()

    logger.info(f"Parsing fasta file: {Path(fastaFile).name}")
    for proteinIdx, (protein, seq) in enumerate(
        readFasta(fastaFile, db, parseId, specialAAs=specialAAs)
    ):
        if proteinIdx % 10000 == 0:
            logger.info(f"Digesting protein {proteinIdx}")
        seenPeptides = set()
        proteinToSeqMap[protein] = seq
        # for peptide in digestfast.getDigestedPeptides(seq, min_len, max_len, pre, not_post, digestion, miscleavages, methionineCleavage):
        for peptide in getDigestedPeptides(
            seq,
            min_len,
            max_len,
            pre,
            not_post,
            post,
            digestion,
            miscleavages,
            methionineCleavage,
        ):
            peptide = peptide
            if useHashKey:
                hashKey = peptide[:6]
            else:
                hashKey = peptide
            if hashKey not in seenPeptides:
                seenPeptides.add(hashKey)
                peptideToProteinMap[hashKey].append(protein)

    if useHashKey:
        return (peptideToProteinMap, proteinToSeqMap)
    else:
        return peptideToProteinMap


def getPeptideToProteinMapFromFile(peptideToProteinMapFile, useHashKey=False):
    if useHashKey:
        logger.info("Hash key not supported yet, continuing without hash key...")
        useHashKey = False
    peptideToProteinMap = collections.defaultdict(list)
    reader = getTsvReader(peptideToProteinMapFile)
    for i, row in enumerate(reader):
        if (i + 1) % 1000000 == 0:
            logger.info(f"Processing peptide {i+1}")

        peptide, proteins = row[0], row[1].split(";")
        if useHashKey:
            sys.exit("Hash key not supported yet...")
            hashKey = peptide[:6]
        else:
            hashKey = peptide
        for protein in proteins:
            peptideToProteinMap[hashKey].append(protein)
    return peptideToProteinMap


def get_proteins(peptideToProteinMap, peptide):
    peptide = peptide  # .replace("I", "L")
    if len(peptideToProteinMap) == 2:
        hashKey = peptide[:6]
        proteins = list()
        if hashKey in peptideToProteinMap[0]:
            for protein in peptideToProteinMap[0][hashKey]:
                # TODO: This does not work correctly for full or partial digestion, since we might find the peptide with the wrong number of enzymatic terminals
                if peptide in peptideToProteinMap[1][protein]:
                    proteins.append(protein)
            proteins = sorted(proteins)
        # else:
        #    logger.warning("Could not find peptide " + peptide + " in fasta database")
        return proteins
    else:
        return peptideToProteinMap.get(peptide, [])


def getAllProteins(peptideToProteinMap):
    seenProteins = set()
    if len(peptideToProteinMap) == 2:
        for _, proteins in peptideToProteinMap[0].items():
            for protein in proteins:
                if protein not in seenProteins:
                    seenProteins.append(protein)
    else:
        for _, proteins in peptideToProteinMap.items():
            for protein in proteins:
                if protein not in seenProteins:
                    seenProteins.append(protein)
    return list(seenProteins)


def getIbaqPeptideToProteinMap(
    fasta_files: List[str], digestion_params_list: List[DigestionParams]
):
    digestion_params_list_ibaq = []
    for digestion_params in digestion_params_list:
        digestion_params.min_length = max([6, digestion_params.min_length])
        digestion_params.max_length = min([30, digestion_params.max_length])
        digestion_params.cleavages = 0
        digestion_params.methionine_cleavage = False
        digestion_params_list_ibaq.append(digestion_params)
    return get_peptide_to_protein_map_from_params(
        fasta_files, digestion_params_list_ibaq
    )


def getNumIbaqPeptidesPerProtein(
    fasta_files: List[str], digestion_params_list: List[DigestionParams]
):
    peptideToProteinMapIbaq = getIbaqPeptideToProteinMap(
        fasta_files, digestion_params_list
    )
    return getNumPeptidesPerProtein(peptideToProteinMapIbaq)


def getNumPeptidesPerProtein(peptideToProteinMap):
    numPeptidesPerProtein = collections.defaultdict(int)
    for peptide, proteins in peptideToProteinMap.items():
        for protein in proteins:
            numPeptidesPerProtein[protein] += 1

    return numPeptidesPerProtein


def getCleavageSites(enzyme):
    if enzyme not in ENZYME_CLEAVAGE_RULES:
        logger.error("Enzyme", enzyme, "not implemented yet")

    pre = ENZYME_CLEAVAGE_RULES[enzyme]["pre"]
    not_post = ENZYME_CLEAVAGE_RULES[enzyme]["not_post"]
    post = ENZYME_CLEAVAGE_RULES[enzyme]["post"]
    return pre, not_post, post


def is_enzymatic_advanced(
    aa1, aa2, pre=["K", "R"], not_post=["P"], post=[], methionineCleavage=True
):
    return (
        aa1 == "-"
        or aa2 == "-"
        or is_enzymatic(aa1, aa2, pre, not_post, post)
        or (methionineCleavage and aa1 == "M")
    )


def is_enzymatic(aa1, aa2, pre, not_post, post):
    return (aa1 in pre and aa2 not in not_post) or (aa2 in post)


def hasMiscleavage(seq, pre=["K", "R"], not_post=["P"], post=[]):
    for i in range(len(seq) - 1):
        if is_enzymatic_advanced(seq[i], seq[i + 1], pre, not_post, post):
            return True
    return False


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


if __name__ == "__main__":
    main(sys.argv[1:])
