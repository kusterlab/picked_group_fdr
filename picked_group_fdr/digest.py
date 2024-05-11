from __future__ import print_function, annotations
from pathlib import Path

import sys
import csv
import itertools
import collections
import logging
from typing import Dict, Iterator, List, Optional

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
    "clostripain": {"pre": ["R"], "not_post": [], "post": []},
    "cyanogen-bromide": {"pre": ["M"], "not_post": [], "post": []},
    "iodosobenzoate": {"pre": ["W"], "not_post": [], "post": []},
    "proline-endopeptidase": {"pre": ["P"], "not_post": [], "post": []},
    "staph-protease": {"pre": ["E"], "not_post": [], "post": []},
    "asp-n": {"pre": [], "not_post": [], "post": ["D"]},
    "lys-c": {"pre": ["K"], "not_post": ["P"], "post": []},
    "lys-cp": {"pre": ["K"], "not_post": [], "post": []},
    "lys-n": {"pre": [], "not_post": [], "post": ["K"]},
    "arg-c": {"pre": ["R"], "not_post": ["P"], "post": []},
    "glu-c": {"pre": ["E"], "not_post": ["P"], "post": []},
    "pepsin-a": {"pre": ["F", "L"], "not_post": ["P"], "post": []},
    "elastase-trypsin-chymotrypsin": {
        "pre": ["A", "L", "I", "V", "F", "K", "R", "W", "F", "Y"],
        "not_post": ["P"],
        "post": [],
    },
    "lysarginase": {"pre": [], "not_post": [], "post": ["K", "R"]},
    "v8-de": {"pre": ["N", "D", "E", "Q"], "not_post": ["P"], "post": []},
}

PeptideToProteinMap = Dict[str, List[str]]


def main(argv):  # pragma: no cover
    args = parse_args()

    digestion_params_list = get_digestion_params_list(args)

    if args.prosit_input:
        writer = get_tsv_writer(args.prosit_input, delimiter=",")
        writer.writerow(
            "modified_sequence,collision_energy,precursor_charge".split(",")
        )

        prosit_input_file_with_proteins = args.prosit_input.replace(
            ".csv", "_with_proteins.csv"
        )
        writer_with_proteins = get_tsv_writer(
            prosit_input_file_with_proteins, delimiter=","
        )
        writer_with_proteins.writerow(
            "modified_sequence,collision_energy,precursor_charge,protein".split(",")
        )

        for peptide, proteins in get_peptide_to_protein_map_from_params(
            args.fasta, digestion_params_list
        ).items():
            if not is_valid_prosit_peptide(peptide):
                continue

            for charge in [2, 3, 4]:
                writer.writerow([peptide, 30, charge])
                writer_with_proteins.writerow([peptide, 30, charge, proteins[0]])

    if args.peptide_protein_map:
        with open(args.peptide_protein_map + ".params.txt", "w") as f:
            f.write(" ".join(sys.argv))

        writer = get_tsv_writer(args.peptide_protein_map, delimiter="\t")
        for peptide, proteins in get_peptide_to_protein_map_from_params(
            args.fasta, digestion_params_list
        ).items():
            writer.writerow([peptide, ";".join(proteins)])

    if args.ibaq_map:
        writer = get_tsv_writer(args.ibaq_map, delimiter="\t")

        num_peptides_per_protein = get_num_ibaq_peptides_per_protein(
            args.fasta, digestion_params_list
        )
        for protein, num_peptides in num_peptides_per_protein.items():
            writer.writerow([protein, num_peptides])


def parse_args():  # pragma: no cover
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
        help="""Fasta file used as input.""",
    )

    apars.add_argument(
        "--prosit_input",
        default=None,
        metavar="M",
        required=False,
        help="""Path to file where to write the prosit input file.""",
    )

    apars.add_argument(
        "--peptide_protein_map",
        default=None,
        metavar="M",
        required=False,
        help="""Write mapping from peptides to all its proteins to 
                the specified file.""",
    )

    apars.add_argument(
        "--ibaq_map",
        default=None,
        metavar="M",
        required=False,
        help="""Write number of peptides per protein to the specified file that meet 
                the iBAQ criteria (6 <= pepLen <= 30, no miscleavages).""",
    )

    add_digestion_arguments(apars)

    # ------------------------------------------------
    args = apars.parse_args()

    return args


def is_valid_prosit_peptide(peptide):
    return len(peptide) <= 30 and "U" not in peptide and "X" not in peptide


def parse_until_first_space(fasta_id: str) -> str:
    return fasta_id.split(" ")[0]


def read_fasta_tide(
    file_path: str, db: str = "target", parse_id=parse_until_first_space
):
    read_fasta_maxquant(file_path, db, parse_id, special_aas=[], decoy_prefix="decoy_")


def read_fasta_maxquant(
    file_path: str,
    db: str = "target",
    parse_id=parse_until_first_space,
    special_aas: Optional[List[str]] = None,
    decoy_prefix: str = "REV__",
):
    if special_aas is None:
        special_aas = ["K", "R"]

    if db not in ["target", "decoy", "concat"]:
        raise ValueError("unknown db mode: %s" % db)

    has_special_aas = len(special_aas) > 0
    name, seq = None, []
    with open(file_path, "r") as fp:
        for line in itertools.chain(fp, [">"]):
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    seq = "".join(seq)
                    if db in ["target", "concat"]:
                        yield (name, seq)

                    if db in ["decoy", "concat"]:
                        rev_seq = seq[::-1]
                        if has_special_aas:
                            rev_seq = swap_special_aas(rev_seq, special_aas)
                        yield (decoy_prefix + name, rev_seq)

                if len(line) > 1:
                    name, seq = parse_id(line[1:]), []
            else:
                seq.append(line)


# from . import digestfast
# read_fasta = digestfast.readFastaMaxQuant
read_fasta = read_fasta_maxquant


def swap_special_aas(seq: str, special_aas: List[str]):
    """Swaps the special AAs with its preceding amino acid, as is done in MaxQuant.

    e.g. special_aas = ['R', 'K'] transforms ABCKDEFRK into ABKCDERKF
    """
    seq = list(seq)
    for i in range(1, len(seq)):
        if seq[i] in special_aas:
            swap_positions(seq, i, i - 1)
    seq = "".join(seq)
    return seq


def swap_positions(seq: str, pos1: int, pos2: int):
    seq[pos1], seq[pos2] = seq[pos2], seq[pos1]


def get_protein_sequences(file_paths: Optional[List[str]], **kwargs):
    if file_paths is None:
        return dict()

    protein_sequences = dict()
    for file_path in file_paths:
        for protein_id, protein_sequence in read_fasta(file_path, **kwargs):
            if (
                protein_id not in protein_sequences
            ):  # keep only first sequence per identifier
                protein_sequences[protein_id] = protein_sequence
    return protein_sequences


def filter_fasta_file(
    fasta_file: str, filtered_fasta_file: str, proteins: List[str], **kwargs
):
    with open(filtered_fasta_file, "w") as f:
        for prot, seq in read_fasta(fasta_file, **kwargs):
            if prot in proteins:
                f.write(">" + prot + "\n" + seq + "\n")
                # f.write('>decoy_' + prot + '\n' + seq[::-1] + '\n')


# @profile
def get_digested_peptides(
    seq: str,
    min_len: int = 6,
    max_len: int = 50,
    pre: List[str] = ["K", "R"],
    not_post: List[str] = ["P"],
    post: List[str] = [],
    digestion: str = "full",
    miscleavages: int = 0,
    methionine_cleavage: bool = True,
):
    if digestion == "none":
        yield from non_specific_digest(seq, min_len, max_len)
    elif digestion == "semi":
        yield from semi_specific_digest(
            seq,
            min_len,
            max_len,
            pre,
            not_post,
            post,
            miscleavages,
            methionine_cleavage,
        )
    else:
        yield from full_digest(
            seq,
            min_len,
            max_len,
            pre,
            not_post,
            post,
            miscleavages,
            methionine_cleavage,
        )


def non_specific_digest(seq, min_len, max_len):
    seq_len = len(seq)
    for i in range(seq_len + 1):
        for j in range(i + min_len, min(seq_len + 1, i + max_len + 1)):
            if j <= seq_len:
                yield seq[i:j]


def semi_specific_digest(
    seq: str,
    min_len: int,
    max_len: int,
    pre: List[str],
    not_post: List[str],
    post: List[str],
    miscleavages: int,
    methionine_cleavage: bool,
):
    seq_len, starts = len(seq), [0]
    methionine_cleavage = methionine_cleavage and seq[0] == "M"
    length_accepted = lambda x: x >= min_len and x <= max_len

    for i in range(seq_len + 1):
        is_cleavage_site = is_enzymatic(
            seq[min([seq_len - 1, i])],
            seq[min([seq_len - 1, i + 1])],
            pre,
            not_post,
            post,
        )
        is_methionine_cleavage_site = i == 0 and methionine_cleavage
        if i == seq_len or is_cleavage_site or is_methionine_cleavage_site:
            # peptides with enzymatic C-terminal (both enzymatic and non-enzymatic N-terminal)
            start = starts[0]
            for j in range(start, min([i + 1, seq_len])):
                pep_len = min([i, seq_len - 1]) - j + 1
                if length_accepted(pep_len):
                    yield (seq[j : i + 1])
            starts.append(i + 1)
            methionine_cleaved = int(starts[0] == 0 and methionine_cleavage)
            if len(starts) > miscleavages + 1 + methionine_cleaved or i == seq_len:
                starts = starts[1 + methionine_cleaved :]
        else:  # peptides with non enzymatic C-terminal
            for start in starts:
                pep_len = i - start + 1
                if length_accepted(pep_len) and i + 1 not in starts:
                    yield (seq[start : i + 1])


def full_digest(
    seq: str,
    min_len: int,
    max_len: int,
    pre: List[str],
    not_post: List[str],
    post: List[str],
    miscleavages: int,
    methionine_cleavage: bool,
):
    seq_len, starts = len(seq), [0]
    methionine_cleavage = methionine_cleavage and seq[0] == "M"

    check_pre = len(pre) > 0
    check_post = len(post) > 0

    cleavage_sites = [0] if methionine_cleavage else []
    # HACK: inline if statement instead of using is_enzymatic because it is ~20% faster
    cleavage_sites.extend(
        [
            i
            for i in range(seq_len)
            if (
                check_pre
                and seq[i] in pre
                and not seq[min([seq_len - 1, i + 1])] in not_post
            )
            or (check_post and seq[min([seq_len - 1, i + 1])] in post)
        ]
    )
    cleavage_sites.append(seq_len)
    for i in cleavage_sites:
        for start in starts:
            pep_len = i - start + 1
            if min_len <= pep_len <= max_len:
                yield (seq[start : i + 1])
        starts.append(i + 1)
        methionine_cleaved = int(starts[0] == 0 and methionine_cleavage)
        if len(starts) > miscleavages + 1 + methionine_cleaved:
            starts = starts[1 + methionine_cleaved :]


def get_peptide_to_protein_map_from_params(
    fasta_files: List[str], digestion_params_list: List[DigestionParams], **kwargs
):
    peptide_to_protein_map = collections.defaultdict(list)
    protein_to_seq_map = dict()
    for fasta_file in fasta_files:
        for params in digestion_params_list:
            # TODO: make sure we do not combine use_hash_key=True with use_hash_key=False
            peptide_to_protein_map_tmp = get_peptide_to_protein_map_from_params_single(
                fasta_file, params, **kwargs
            )

            # TODO: refactor peptide_to_protein_map as a class to get rid of this check
            if len(peptide_to_protein_map_tmp) == 2:
                (
                    peptide_to_protein_map_tmp,
                    protein_to_seq_map_tmp,
                ) = peptide_to_protein_map_tmp
                protein_to_seq_map |= protein_to_seq_map_tmp

            if len(peptide_to_protein_map) == 0:
                peptide_to_protein_map = peptide_to_protein_map_tmp
            else:
                for peptide, proteins in peptide_to_protein_map_tmp.items():
                    peptide_to_protein_map[peptide].extend(proteins)

    if len(protein_to_seq_map) > 0:
        return peptide_to_protein_map, protein_to_seq_map

    return peptide_to_protein_map


def get_peptide_to_protein_map_from_params_single(
    fasta_file: str, params: DigestionParams, **kwargs
):
    pre, not_post, post = get_cleavage_sites(params.enzyme)
    return get_peptide_to_protein_map(
        fasta_file,
        params.db,
        digestion=params.digestion,
        min_len=params.min_length,
        max_len=params.max_length,
        pre=pre,
        not_post=not_post,
        post=post,
        miscleavages=params.cleavages,
        methionine_cleavage=params.methionine_cleavage,
        use_hash_key=params.use_hash_key,
        special_aas=params.special_aas,
        **kwargs,
    )


def get_peptide_to_protein_map(
    fasta_file: str,
    db: str = "concat",
    min_len: int = 6,
    max_len: int = 52,
    pre: List[str] = ["K", "R"],
    not_post: List[str] = ["P"],
    post: List[str] = [],
    digestion: str = "full",
    miscleavages: int = 2,
    methionine_cleavage: bool = True,
    use_hash_key: bool = False,
    special_aas: List[str] = ["K", "R"],
    parse_id=parse_until_first_space,
):
    peptide_to_protein_map = collections.defaultdict(list)
    protein_to_seq_map = dict()

    logger.info(f"Parsing fasta file: {Path(fasta_file).name}")
    for protein_idx, (protein, seq) in enumerate(
        read_fasta(fasta_file, db, parse_id, special_aas=special_aas)
    ):
        if protein_idx % 10000 == 0:
            logger.info(f"Digesting protein {protein_idx}")
        seen_peptides = set()
        protein_to_seq_map[protein] = seq
        for peptide in get_digested_peptides(
            seq,
            min_len,
            max_len,
            pre,
            not_post,
            post,
            digestion,
            miscleavages,
            methionine_cleavage,
        ):
            peptide = peptide
            if use_hash_key:
                hash_key = peptide[:6]
            else:
                hash_key = peptide
            if hash_key not in seen_peptides:
                seen_peptides.add(hash_key)
                peptide_to_protein_map[hash_key].append(protein)

    if use_hash_key:
        return (peptide_to_protein_map, protein_to_seq_map)
    else:
        return peptide_to_protein_map


def merge_peptide_to_protein_maps(peptide_protein_maps: Iterator[PeptideToProteinMap]):
    peptide_to_protein_map = collections.defaultdict(list)
    for peptide_protein_map in peptide_protein_maps:
        for peptide, proteins in peptide_protein_map.items():
            peptide_to_protein_map[peptide].extend(proteins)
    return peptide_to_protein_map


def get_peptide_to_protein_map_from_file(
    peptide_to_protein_map_file, use_hash_key=False
):
    if use_hash_key:
        logger.info("Hash key not supported yet, continuing without hash key...")
        use_hash_key = False
    peptide_to_protein_map = collections.defaultdict(list)
    reader = get_tsv_reader(peptide_to_protein_map_file)
    for i, row in enumerate(reader):
        if (i + 1) % 1000000 == 0:
            logger.info(f"Processing peptide {i+1}")

        peptide, proteins = row[0], row[1].split(";")
        if use_hash_key:
            raise NotImplementedError("Hash key not supported yet...")
            hash_key = peptide[:6]
        else:
            hash_key = peptide
        for protein in proteins:
            peptide_to_protein_map[hash_key].append(protein)
    return peptide_to_protein_map


def get_proteins(peptide_to_protein_map, peptide: str):
    if len(peptide_to_protein_map) == 2:
        hash_key = peptide[:6]
        proteins = list()
        if hash_key in peptide_to_protein_map[0]:
            for protein in peptide_to_protein_map[0][hash_key]:
                # TODO: This does not work correctly for full or partial digestion,
                # since we might find the peptide with the wrong number of enzymatic terminals
                if peptide in peptide_to_protein_map[1][protein]:
                    proteins.append(protein)
            proteins = sorted(proteins)
        # else:
        #    logger.warning("Could not find peptide " + peptide + " in fasta database")
        return proteins
    else:
        return peptide_to_protein_map.get(peptide, [])


def get_all_proteins(peptide_to_protein_map):
    seen_proteins = set()
    if len(peptide_to_protein_map) == 2:
        peptide_to_protein_map = peptide_to_protein_map[0]
    for _, proteins in peptide_to_protein_map.items():
        for protein in proteins:
            if protein not in seen_proteins:
                seen_proteins.append(protein)
    return list(seen_proteins)


def get_ibaq_peptide_to_protein_map(
    fasta_files: List[str], digestion_params_list: List[DigestionParams], **kwargs
):
    digestion_params_list_ibaq = []
    for digestion_params in digestion_params_list:
        digestion_params.min_length = max([6, digestion_params.min_length])
        digestion_params.max_length = min([30, digestion_params.max_length])
        digestion_params.cleavages = 0
        digestion_params.methionine_cleavage = False
        digestion_params_list_ibaq.append(digestion_params)
    return get_peptide_to_protein_map_from_params(
        fasta_files, digestion_params_list_ibaq, **kwargs
    )


def get_num_ibaq_peptides_per_protein_from_args(args, peptide_to_protein_maps, **kwargs):
    digestion_params_list = get_digestion_params_list(args)
    if args.fasta:
        logger.info("In silico protein digest for iBAQ")
        num_ibaq_peptides_per_protein = get_num_ibaq_peptides_per_protein(
            args.fasta, digestion_params_list, **kwargs
        )
    elif args.peptide_protein_map:
        logger.warning("Found peptide_protein_map (instead of fasta input): ")
        logger.warning(
            "- calculating iBAQ values using all peptides in peptide_protein_map."
        )
        logger.warning("- cannot compute sequence coverage.")
        num_ibaq_peptides_per_protein = get_num_peptides_per_protein(
            merge_peptide_to_protein_maps(peptide_to_protein_maps)
        )
    else:
        raise ValueError(
            "No fasta or peptide to protein mapping file detected, please specify either the --fasta or --peptide_protein_map flags"
        )
    return num_ibaq_peptides_per_protein


def get_num_ibaq_peptides_per_protein(
    fasta_files: List[str], digestion_params_list: List[DigestionParams], **kwargs
) -> Dict[str, int]:
    peptide_to_protein_map_ibaq = get_ibaq_peptide_to_protein_map(
        fasta_files, digestion_params_list, **kwargs
    )
    return get_num_peptides_per_protein(peptide_to_protein_map_ibaq)


def get_num_peptides_per_protein(peptide_to_protein_map) -> Dict[str, int]:
    num_peptides_per_protein = collections.defaultdict(int)
    for _, proteins in peptide_to_protein_map.items():
        for protein in proteins:
            num_peptides_per_protein[protein] += 1

    return num_peptides_per_protein


def get_cleavage_sites(enzyme):
    if enzyme not in ENZYME_CLEAVAGE_RULES:
        logger.error("Enzyme", enzyme, "not implemented yet")

    pre = ENZYME_CLEAVAGE_RULES[enzyme]["pre"]
    not_post = ENZYME_CLEAVAGE_RULES[enzyme]["not_post"]
    post = ENZYME_CLEAVAGE_RULES[enzyme]["post"]
    return pre, not_post, post


def is_enzymatic_advanced(
    aa1: str,
    aa2: str,
    pre: List[str] = ["K", "R"],
    not_post: List[str] = ["P"],
    post: List[str] = [],
    methionine_cleavage: bool = True,
):
    return (
        aa1 == "-"
        or aa2 == "-"
        or is_enzymatic(aa1, aa2, pre, not_post, post)
        or (methionine_cleavage and aa1 == "M")
    )


def is_enzymatic(aa1, aa2, pre, not_post, post):
    return (aa1 in pre and aa2 not in not_post) or (aa2 in post)


def has_miscleavage(seq, pre=["K", "R"], not_post=["P"], post=[]):
    for i in range(len(seq) - 1):
        if is_enzymatic_advanced(seq[i], seq[i + 1], pre, not_post, post):
            return True
    return False


def get_tsv_reader(filename, delimiter="\t"):
    # Python 3
    if sys.version_info[0] >= 3:
        return csv.reader(open(filename, "r", newline=""), delimiter=delimiter)
    # Python 2
    else:
        return csv.reader(open(filename, "rb"), delimiter=delimiter)


def get_tsv_writer(filename, delimiter="\t"):
    # Python 3
    if sys.version_info[0] >= 3:
        return csv.writer(open(filename, "w", newline=""), delimiter=delimiter)
    # Python 2
    else:
        return csv.writer(open(filename, "wb"), delimiter=delimiter)


if __name__ == "__main__":
    main(sys.argv[1:])
