import sys
import os
import logging
from typing import Dict, List, Optional

from .update_fragpipe_results import (
    write_fragpipe_combined_protein_file,
)
from .. import __version__, __copyright__
from .. import digest
from .. import protein_annotation
from ..picked_group_fdr import ArgumentParserWithLogger
from ..parsers import maxquant
from ..quant.sage import add_precursor_quants_multiple
from ..protein_annotation import ProteinAnnotation
from ..protein_groups import ProteinGroups
from ..results import ProteinGroupResults

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
        help="""Fasta file(s) to create mapping from peptides to proteins.
                This should not contain the decoy sequences, unless you set the 
                --fasta_contains_decoys flag.""",
    )

    apars.add_argument(
        "--fasta_contains_decoys",
        help="Set this flag if your fasta file already contains decoy protein sequences.",
        action="store_true",
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
        "--output_folder",
        default="./",
        metavar="DIR",
        help="""Path to output folder. Folder is created if it does not exist yet.""",
    )

    # ------------------------------------------------
    args = apars.parse_args(argv)

    return args


def main(argv):
    logger.info(f"SageQuantification version {__version__}\n{__copyright__}")
    logger.info(
        f'Issued command: {os.path.basename(__file__)} {" ".join(map(str, argv))}'
    )

    args = parseArgs(argv)

    protein_groups = ProteinGroups.from_mq_protein_groups_file(args.protein_groups)
    protein_groups.create_index()

    db = "target" if args.fasta_contains_decoys else "concat"
    protein_annotations = protein_annotation.get_protein_annotations_multiple(
        args.fasta, db=db, parse_id=digest.parse_until_first_space
    )

    # create a fresh ProteinGroupResults object for combined_protein.tsv
    protein_group_results = maxquant.parse_mq_protein_groups_file(args.protein_groups)
    generate_fragpipe_combined_protein_file(
        args.sage_results,
        args.sage_lfq_tsv,
        protein_groups,
        protein_group_results,
        protein_annotations,
        args.output_folder,
    )


def generate_fragpipe_combined_protein_file(
    sage_results_files: List[str],
    combined_ion_file: str,
    protein_groups: ProteinGroups,
    protein_group_results: ProteinGroupResults,
    protein_annotations: Dict[str, ProteinAnnotation],
    output_folder: Optional[str] = None,
    psm_fdr_cutoff: float = 0.01,
    discard_shared_peptides: bool = True,
):
    """Generate protein group results in FragPipe's combined_protein.tsv format."""

    protein_group_results, post_err_probs = add_precursor_quants_multiple(
        sage_results_files,
        combined_ion_file,
        protein_groups,
        protein_group_results,
        discard_shared_peptides,
    )

    fragpipe_combined_protein_file_out = f"{output_folder}/combined_protein.tsv"

    write_fragpipe_combined_protein_file(
        fragpipe_combined_protein_file_out,
        protein_groups,
        protein_group_results,
        protein_annotations,
        post_err_probs,
        psm_fdr_cutoff,
    )


if __name__ == "__main__":
    main(sys.argv[1:])
