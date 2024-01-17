import json
import sys
import os
import logging
from typing import Dict, List


from .. import __version__, __copyright__
from .. import digest
from .. import protein_annotation
from ..picked_group_fdr import ArgumentParserWithLogger
from ..protein_annotation import ProteinAnnotation, parse_uniprot_id

# hacky way to get package logger when running as module
logger = logging.getLogger(__package__ + "." + __file__)


def parse_args(argv):
    import argparse

    apars = ArgumentParserWithLogger(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    apars.add_argument(
        "--log_file",
        default=None,
        metavar="L",
        required=True,
        help="""File with log output of Picked Group FDR.""",
    )

    apars.add_argument(
        "--fasta",
        default=None,
        metavar="F",
        nargs="+",
        required=True,
        help="""Fasta file(s) to create mapping from peptides to proteins.
                This should not contain the decoy sequences, unless you set the 
                --fasta_contains_decoys flag.""",
    )

    digest.add_digestion_arguments(apars)

    # ------------------------------------------------
    args = apars.parse_args(argv)

    return args


def main(argv):
    logger.info(f"AnalyzeMissingPeptides version {__version__}\n{__copyright__}")
    logger.info(
        f'Issued command: {os.path.basename(__file__)} {" ".join(map(str, argv))}'
    )

    args = parse_args(argv)

    db = "target" if args.fasta_contains_decoys else "concat"
    protein_annotations = protein_annotation.get_protein_annotations_multiple(
        args.fasta, db=db, parse_id=parse_uniprot_id
    )
    protein_sequences = digest.get_protein_sequences(
        args.fasta, db=db, parse_id=parse_uniprot_id
    )

    missing_peptides: Dict[str, List[ProteinAnnotation]] = {}
    with open(args.log_file, 'r') as f:
        for line in f:
            if "Missing peptide: " in line:
                peptide_protein_string = line.split("Missing peptide: ")[-1]
                peptide = peptide_protein_string.split()[0]
                proteins = " ".join(peptide_protein_string.split()[1:])
                proteins = json.loads(proteins.replace("'", "\""))

                mapped_protein_annotations = []
                for fasta_header in proteins:
                    mapped_protein_annotations.append(
                        ProteinAnnotation(
                            id=digest.parse_until_first_space(fasta_header),
                            fasta_header=fasta_header,
                            uniprot_id=parse_uniprot_id(fasta_header)
                        )
                    )
                missing_peptides[peptide] = mapped_protein_annotations
    
    missing_uniprot_id = 0
    for peptide, mapped_protein_annotations in missing_peptides.items():
        found_in_fasta = False
        for pa in mapped_protein_annotations:
            if pa.uniprot_id in protein_annotations:
                found_in_fasta = True
            
        if not found_in_fasta:
            missing_uniprot_id += 1

    num_missing_peptides = len(missing_peptides)
    print(f"Missing unique peptides: {num_missing_peptides}")
    print(f"- Protein identifier missing in fasta: {missing_uniprot_id} ({get_percentage_string(missing_uniprot_id, num_missing_peptides)})")


def get_percentage_string(subset_size, total_size):
    return f"{subset_size / total_size * 100:.2g}%"


if __name__ == "__main__":
    main(sys.argv[1:])
