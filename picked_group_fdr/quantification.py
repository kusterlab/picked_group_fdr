import sys
import os
import logging
from typing import Dict, List, Tuple

from . import writers
from . import digest
from . import digestion_params
from . import protein_annotation
from . import picked_group_fdr
from .parsers import maxquant
from .parsers import parsers
from .quant import maxquant as mq_quant
from .columns import protein_annotations as pa
from .columns.triqler import init_triqler_params

logger = logging.getLogger(__name__)


def parse_args(argv):
    import argparse

    apars = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    apars.add_argument(
        "--mq_evidence",
        default=None,
        metavar="EV",
        required=True,
        nargs="+",
        help="""MaxQuant evidence file.""",
    )

    apars.add_argument(
        "--mq_protein_groups",
        default=None,
        metavar="PG",
        required=True,
        help="""MaxQuant protein groups file.""",
    )

    apars.add_argument(
        "--protein_groups_out",
        default=None,
        metavar="PG",
        required=True,
        help="""Protein groups output file, mimicks a subset of the MQ protein groups columns.
                                """,
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

    add_quant_arguments(apars)

    # ------------------------------------------------
    args = apars.parse_args(argv)

    return args


def add_quant_arguments(apars) -> None:
    apars.add_argument(
        "--psm_fdr_cutoff",
        default=0.01,
        metavar="C",
        type=float,
        help="PSM-level FDR threshold used for filtering PSMs for quantification.",
    )

    apars.add_argument(
        "--lfq_min_peptide_ratios",
        default=2,
        type=int,
        metavar="M",
        help="""Minimum number of common peptides between two samples
                to qualify for calculating a peptide ratio in LFQ.""",
    )

    apars.add_argument(
        "--lfq_stabilize_large_ratios",
        help="""Apply stabilization of large ratios in LFQ as described
                in the MaxLFQ paper.""",
        action="store_false",
    )

    apars.add_argument(
        "--num_threads",
        default=1,
        type=int,
        metavar="T",
        help="""Maximum number of threads to use. Currently only speeds up the MaxLFQ part.""",
    )

    apars.add_argument(
        "--experimental_design_file",
        metavar="L",
        help="""Tab separated file in MaxQuant experimental design format. It should
                contain the column names: "Name", "Fraction" and "Experiment". 
                Optionally, a column named "Condition" can be added to enable 
                differential expression analysis with Triqler.""",
    )

    apars.add_argument(
        "--file_list_file",
        metavar="L",
        help="""[deprecated in favor of --experimental_design_file] Tab separated file 
                with lines of the format (third and fourth columns are optional): 
                raw_file <tab> condition <tab> experiment <tab> fraction.""",
    )


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
    db = "target" if args.fasta_contains_decoys else "concat"

    (
        peptide_to_protein_maps,
        num_ibaq_peptides_per_protein,
    ) = get_peptide_to_protein_maps(args, parse_id=parse_id)

    protein_annotations = protein_annotation.get_protein_annotations_multiple(
        args.fasta, db=db, parse_id=parse_id
    )
    protein_sequences = digest.get_protein_sequences(
        args.fasta, db=db, parse_id=parse_id
    )

    protein_group_results = maxquant.parse_mq_protein_groups_file(
        args.mq_protein_groups,
        additional_headers=pa.MQ_PROTEIN_ANNOTATION_HEADERS,
    )

    experimental_design = get_experimental_design(args)

    params = init_triqler_params(experimental_design)
    protein_groups_writer = writers.MaxQuantProteinGroupsWriter(
        num_ibaq_peptides_per_protein,
        protein_annotations,
        protein_sequences,
        args.lfq_min_peptide_ratios,
        args.lfq_stabilize_large_ratios,
        args.num_threads,
        params,
    )

    protein_group_results, post_err_probs = mq_quant.add_precursor_quants(
        protein_group_results,
        args.mq_evidence,
        peptide_to_protein_maps,
        experimental_design,
        discard_shared_peptides=True,
    )

    protein_group_results = protein_groups_writer.append_quant_columns(
        protein_group_results, post_err_probs, args.psm_fdr_cutoff
    )

    protein_group_results.write(args.protein_groups_out)

    logger.info(
        f"Protein group results have been written to: {args.protein_groups_out}"
    )


def get_experimental_design(args):
    experimental_design = None
    if args.experimental_design_file:
        experimental_design = parsers.parse_mq_experimental_design(
            args.experimental_design_file
        )
    elif args.file_list_file:
        experimental_design = parsers.parse_triqler_file_list(args.file_list_file)
    return experimental_design


def get_peptide_to_protein_maps(
    args, **kwargs
) -> Tuple[List[digest.PeptideToProteinMap], Dict[str, int]]:
    logger.info("Loading peptide to protein map...")

    digestion_params_list = digestion_params.get_digestion_params_list(args)
    peptide_to_protein_maps = picked_group_fdr.get_peptide_to_protein_maps(
        args.fasta,
        args.peptide_protein_map,
        digestion_params_list,
        args.mq_protein_groups,
        **kwargs,
    )

    num_ibaq_peptides_per_protein = digest.get_num_ibaq_peptides_per_protein_from_args(
        args, peptide_to_protein_maps
    )

    return peptide_to_protein_maps, num_ibaq_peptides_per_protein


if __name__ == "__main__":
    main(sys.argv[1:])
