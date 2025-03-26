"""Performs quantification of previously generated protein groups and adds columns to a new output file.

Currently only works with MaxQuant proteinGroups.txt.

Usage:

    python -m picked_group_fdr.quantification \
        --mq_evidence evidence.txt \
        --mq_protein_groups proteinGroups.txt \
        --protein_groups_out proteinGroups_with_quant.txt \
        --fasta db.fasta
"""

import sys
import os
import argparse
import logging
from typing import Dict, List, Tuple

from . import writers
from . import digest
from . import digestion_params
from . import protein_annotation
from . import peptide_protein_map
from .parsers import maxquant
from .parsers import fragpipe
from .parsers import parsers
from .protein_groups import ProteinGroups
from .scoring_strategy import ProteinScoringStrategy
from .results import ProteinGroupResults

# hacky way to get package logger when running as module
logger = logging.getLogger(__package__ + "." + __file__)


def parse_args(argv):
    apars = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    apars.add_argument(
        "--mq_evidence",
        default=None,
        metavar="EV",
        nargs="+",
        help="""MaxQuant evidence file.""",
    )

    apars.add_argument(
        "--fragpipe_psm",
        default=None,
        metavar="PSM",
        nargs="+",
        help="""Fragpipe psm.tsv output file(s); alternative for 
                --mq_evidence.""",
    )

    apars.add_argument(
        "--combined_ion",
        default=None,
        metavar="I",
        nargs="+",
        help="""Path to combined_ion.tsv produced by IonQuant/FragPipe. This enables
                quantification of protein groups by PickedGroupFDR.""",
    )

    apars.add_argument(
        "--mq_protein_groups",
        default=None,
        metavar="PG",
        help="""MaxQuant protein groups file.""",
    )

    apars.add_argument(
        "--combined_protein",
        default=None,
        metavar="PG",
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
        "--output_format",
        default="auto",
        metavar="PG",
        help="""Protein groups output format. Options are "auto", "maxquant" and 
                "fragpipe". "auto": decide based on input file format, uses
                "maxquant" if no suitable format is known; "maxquant": proteinGroups.txt
                format; "fragpipe": combined_protein.tsv format.""",
    )

    apars.add_argument(
        "--protein_group_fdr_threshold",
        default=0.01,
        metavar="C",
        type=float,
        help="""Protein group-level FDR threshold used for computing the number of
                protein groups at the given threshold. Note that this does not filter the
                protein_groups_out file for the specified FDR. Use the
                picked_group_fdr.pipeline.filter_fdr_maxquant module for this.""",
    )

    apars.add_argument(
        "--peptide_protein_map",
        default=None,
        metavar="M",
        nargs="+",
        help="""Tab-separated file with mapping from peptides to proteins; alternative for 
                --fasta flag if digestion is time consuming.""",
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
        "--gene_level",
        help="""Do quantification on gene-level instead of on protein group level""",
        action="store_true",
    )

    apars.add_argument(
        "--suppress_missing_peptide_warning",
        help="Suppress missing peptide warning when mapping peptides to proteins.",
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
        "--skip_lfq",
        help="Skip LFQ quantification, this significantly speeds up quantification.",
        action="store_true",
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


def do_quantification(
    score_type: ProteinScoringStrategy,
    args: argparse.Namespace,
    protein_group_results: ProteinGroupResults,
    protein_annotations: Dict[str, protein_annotation.ProteinAnnotation],
    use_pseudo_genes: bool,
    peptide_to_protein_maps: List[digest.PeptideToProteinMap],
    suppress_missing_peptide_warning: bool,
) -> Tuple[ProteinGroupResults, writers.ProteinGroupsWriter, List]:
    if not score_type.can_do_quantification():
        logger.warning(
            "Skipping quantification: need input file with precursor quantifications."
        )
        return protein_group_results

    logger.info("Preparing for quantification")

    score_origin = score_type.score_origin.long_description().lower()
    if args.output_format == "auto":
        if score_origin in ["maxquant", "fragpipe"]:
            args.output_format = score_origin
        else:
            args.output_format = "maxquant"

    parse_id = digest.parse_until_first_space
    if args.gene_level and not use_pseudo_genes:
        parse_id = protein_annotation.parse_gene_name_func
    elif args.fasta_use_uniprot_id:
        parse_id = protein_annotation.parse_uniprot_id

    protein_groups_writer = writers.get_protein_groups_output_writer(
        protein_group_results,
        args.output_format,
        args,
        protein_annotations,
        parse_id,
        peptide_to_protein_maps,
    )

    protein_groups = ProteinGroups.from_protein_group_results(protein_group_results)
    experimental_design = get_experimental_design(args)
    discard_shared_peptides = True
    protein_group_results, post_err_probs = score_type.get_quantification_parser()(
        score_type.get_evidence_file(args),
        score_type.get_quantification_file(args),
        protein_groups,
        protein_group_results,
        peptide_to_protein_maps,
        experimental_design,
        discard_shared_peptides,
        score_type=score_type,
        suppress_missing_peptide_warning=suppress_missing_peptide_warning,
    )

    return protein_group_results, protein_groups_writer, post_err_probs


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

    peptide_to_protein_maps = peptide_protein_map.get_peptide_to_protein_maps_from_args(
        args, use_pseudo_genes=False
    )

    protein_annotations = dict()
    if args.fasta:
        protein_annotations = protein_annotation.get_protein_annotations_multiple(
            args.fasta, db=db, parse_id=parse_id
        )

    score_type_string = "bestPEP"
    if peptide_to_protein_maps == [None]:
        score_type_string = "no_remap bestPEP"
    
    if args.mq_protein_groups:
        protein_group_results = maxquant.parse_mq_protein_groups_file(
            args.mq_protein_groups
        )
    elif args.combined_protein:
        protein_group_results = fragpipe.parse_fragpipe_combined_protein_file(
            args.combined_protein
        )
        score_type_string = "FragPipe " + score_type_string
    else:
        raise ValueError("Need either --mq_protein_groups or --combined_protein as input to read protein groups from.")

    score_type = ProteinScoringStrategy(score_type_string)

    use_pseudo_genes = False
    (
        protein_group_results,
        protein_groups_writer,
        post_err_probs,
    ) = do_quantification(
        score_type,
        args,
        protein_group_results,
        protein_annotations,
        use_pseudo_genes,
        peptide_to_protein_maps,
        args.suppress_missing_peptide_warning,
    )

    apply_filename_suffix = False
    method_config = None
    writers.finalize_output(
        protein_group_results,
        protein_groups_writer,
        post_err_probs,
        args.protein_groups_out,
        args.psm_fdr_cutoff,
        apply_filename_suffix,
        method_config,
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


if __name__ == "__main__":
    main(sys.argv[1:])
