import sys
import os
import logging



from .. import __version__, __copyright__
from .. import digest
from .. import writers
from .. import protein_annotation
from ..picked_group_fdr import ArgumentParserWithLogger
from ..quantification import add_quant_arguments
from ..parsers import maxquant
from ..quant.sage import add_precursor_quants_multiple
from ..protein_groups import ProteinGroups

# hacky way to get package logger when running as module
logger = logging.getLogger(__package__ + "." + __file__)


def parseArgs(argv):
    import argparse

    apars = ArgumentParserWithLogger(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    apars.add_argument(
        "--sage_results",
        default=None,
        metavar="PSM",
        nargs="+",
        help="""Sage results.sage.tsv output file(s); alternative for 
                --mq_evidence.""",
    )

    apars.add_argument(
        "--fasta",
        default=None,
        metavar="F",
        nargs="+",
        required=True,
        help="""Fasta file(s) for protein annotations (gene names, etc.).
                This should not contain the decoy sequences, unless you set the 
                --fasta_contains_decoys flag.""",
    )

    apars.add_argument(
        "--protein_groups",
        default=None,
        metavar="PG",
        required=True,
        help="Path to proteinGroups.txt produced by PickedGroupFDR.",
    )

    apars.add_argument(
        "--sage_lfq_tsv",
        default=None,
        metavar="I",
        required=True,
        help="""Path to lfq.tsv produced by Sage. This enables
                quantification of protein groups by PickedGroupFDR in 
                combined_protein.tsv.""",
    )

    apars.add_argument(
        "--protein_groups_out",
        default="./combined_protein.tsv",
        metavar="PG",
        help="""Protein groups output file.""",
    )

    apars.add_argument(
        "--output_format",
        default="maxquant",
        metavar="PG",
        help="""Protein groups output format. Options are "maxquant" (proteinGroups.txt 
                format) and "fragpipe" (combined_protein.tsv format).""",
    )

    add_quant_arguments(apars)
    digest.add_digestion_arguments(apars)

    # ------------------------------------------------
    args = apars.parse_args(argv)

    return args


def main(argv):
    logger.info(f"SageQuantification version {__version__}\n{__copyright__}")
    logger.info(
        f'Issued command: {os.path.basename(__file__)} {" ".join(map(str, argv))}'
    )

    args = parseArgs(argv)

    db = "target" if args.fasta_contains_decoys else "concat"
    parse_id = digest.parse_until_first_space
    protein_annotations = protein_annotation.get_protein_annotations_multiple(
        args.fasta, db=db, parse_id=parse_id
    )

    protein_group_results = maxquant.parse_mq_protein_groups_file(args.protein_groups)

    discard_shared_peptides = True
    psm_fdr_cutoff = 0.01
    peptide_to_protein_maps = dict()

    protein_groups = ProteinGroups.from_protein_group_results(protein_group_results)
    protein_groups.create_index()

    protein_group_results, post_err_probs = add_precursor_quants_multiple(
        args.sage_results,
        args.sage_lfq_tsv,
        protein_groups,
        protein_group_results,
        peptide_to_protein_maps=None,
        file_list_file=args.file_list_file,
        discard_shared_peptides=discard_shared_peptides,
    )

    protein_groups_writer = writers.get_protein_groups_output_writer(
        protein_group_results,
        args.output_format,
        args,
        protein_annotations,
        parse_id,
        peptide_to_protein_maps,
    )

    protein_group_results = protein_groups_writer.append_quant_columns(
        protein_group_results, post_err_probs, psm_fdr_cutoff
    )

    protein_groups_writer.write(protein_group_results, args.protein_groups_out)

    logger.info(f"Protein group results have been written to: {args.protein_groups_out}")


if __name__ == "__main__":
    main(sys.argv[1:])
