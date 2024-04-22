"""Adds iBAQ intensity columns to combined_protein.tsv of FragPipe.

Example:

    python -m picked_group_fdr.pipeline.add_ibaq_columns \
        --protein_groups_in data/fragpipe_example/results_from_combined_ion/combined_protein.tsv \
        --protein_groups_out data/fragpipe_example/results_from_combined_ion/combined_protein_with_ibaq.tsv \
        --fasta data/fragpipe_example/iprg2016_with_labels.fasta

"""

import sys
import os
import logging
from typing import Dict

import pandas as pd

from .. import digest
from .. import digestion_params
from .. import protein_annotation

# hacky way to get package logger when running as module
logger = logging.getLogger(__package__ + "." + __file__)


def parse_args(argv):
    import argparse

    apars = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    apars.add_argument(
        "--protein_groups_in",
        default=None,
        metavar="PG",
        required=True,
        help="""Protein groups file in combined_protein.tsv MSFragger format.""",
    )

    apars.add_argument(
        "--protein_groups_out",
        default=None,
        metavar="PG",
        required=True,
        help="""Protein groups output file, mimicks a subset of the MQ protein groups columns.""",
    )

    apars.add_argument(
        "--peptide_protein_map",
        default=None,
        metavar="M",
        help="""TSV file with mapping from peptides to proteins.""",
    )

    apars.add_argument(
        "--fasta",
        default=None,
        metavar="F",
        nargs="+",
        help="""Fasta file(s) to create mapping from peptides to proteins.
                This should not contain the decoy sequences, unless you set the 
                --fasta_contains_decoys flag.""",
    )

    apars.add_argument(
        "--fasta_use_uniprot_id",
        help="""Parse protein identifiers in the fasta file as UniProt IDs, 
                i.e. Q9UM47 for the protein identifier sp|Q9UM47|NOTC3_HUMAN""",
        action="store_true",
    )

    apars.add_argument(
        "--gene_level",
        help="""Do quantification on gene-level instead of on protein group level""",
        action="store_true",
    )

    digestion_params.add_digestion_arguments(apars)

    # ------------------------------------------------
    args = apars.parse_args(argv)

    return args


def main(argv) -> None:
    logger.info(
        f'Issued command: {os.path.basename(__file__)} {" ".join(map(str, argv))}'
    )

    args = parse_args(argv)

    parse_id = digest.parse_until_first_space
    if args.gene_level:
        parse_id = protein_annotation.parse_gene_name_func
    elif args.fasta_use_uniprot_id:
        parse_id = protein_annotation.parse_uniprot_id

    num_ibaq_peptides_per_protein = digest.get_num_ibaq_peptides_per_protein_from_args(
        args, [], parse_id=parse_id
    )

    df = pd.read_csv(args.protein_groups_in, sep="\t")

    df = add_ibaq_columns(df, num_ibaq_peptides_per_protein)

    df.to_csv(args.protein_groups_out, sep="\t", index=False)

    logger.info(
        f"Protein group results have been written to: {args.protein_groups_out}"
    )


def add_ibaq_columns(df: pd.DataFrame, num_ibaq_peptides_per_protein: Dict[str, int]):
    df["Num iBAQ peptides"] = df["Protein"].apply(
        get_num_ibaq_peptides,
        num_ibaq_peptides_per_protein=num_ibaq_peptides_per_protein,
    )

    summed_intensity_cols = df.filter(regex="^(?!.*MaxLFQ )\S+ Intensity$").columns
    ibaq_intensity_cols = list(
        summed_intensity_cols.map(lambda x: x.replace(" Intensity", " iBAQ Intensity"))
    )
    df.loc[:, ibaq_intensity_cols] = (
        df.loc[:, summed_intensity_cols].values
        / df["Num iBAQ peptides"].values[:, None]
    )
    return df


def get_num_ibaq_peptides(
    protein_ids: str, num_ibaq_peptides_per_protein: Dict[str, int]
):
    numTheoreticalPeptides = [
        num_ibaq_peptides_per_protein.get(p, 1) for p in protein_ids.split(";")
    ]
    leadingProteinNumPeptides = numTheoreticalPeptides[0]
    return leadingProteinNumPeptides


if __name__ == "__main__":
    main(sys.argv[1:])
