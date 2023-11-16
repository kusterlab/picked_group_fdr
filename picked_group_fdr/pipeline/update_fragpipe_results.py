from pathlib import Path
import sys
import os
import logging
from typing import Union

from picked_group_fdr.protein_groups import ProteinGroups

from ..parsers import tsv, fragpipe

from .. import __version__, __copyright__
from .. import helpers
from ..picked_group_fdr import ArgumentParserWithLogger

# hacky way to get the package logger instead of just __main__ when running as python -m picked_group_fdr.pipeline.update_evidence_from_pout ...
logger = logging.getLogger(__package__ + "." + __file__)


def parseArgs(argv):
    import argparse

    apars = ArgumentParserWithLogger(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    apars.add_argument(
        "--fragpipe_psm",
        default=None,
        metavar="PSM",
        nargs="+",
        required=True,
        help="""Fragpipe psm.tsv output file(s).""",
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
        "--protein_groups",
        default=None,
        metavar="PG",
        required=True,
        help="proteinGroups.txt produced by PickedGroupFDR",
    )

    apars.add_argument(
        "--output_folder",
        default=None,
        metavar="DIR",
        help="""Output folder (optional). If this argument is not specified, the
                original psm.tsv and protein.tsv are overwritten, keeping copies
                of the original files as psm.original.tsv and protein.original.tsv.""",
    )

    # ------------------------------------------------
    args = apars.parse_args(argv)

    return args


def main(argv):
    logger.info(f"UpdateFragPipeResults version {__version__}\n{__copyright__}")
    logger.info(
        f'Issued command: {os.path.basename(__file__)} {" ".join(map(str, argv))}'
    )

    args = parseArgs(argv)

    protein_groups = ProteinGroups.from_mq_protein_groups_file(args.protein_groups)
    protein_groups.create_index()

    for fragpipe_psm_file in args.fragpipe_psm:
        update_fragpipe_psm_file_single(
            fragpipe_psm_file, protein_groups, args.output_folder
        )


def update_fragpipe_psm_file_single(
    fragpipe_psm_file: str,
    protein_groups: ProteinGroups,
    output_folder: Union[str, None] = None,
    discard_shared_peptides: bool = True,
) -> None:
    missing_peptides_in_protein_groups = 0
    peptides_not_mapping_to_leading_protein = 0
    shared_peptide_precursors, unique_peptide_precursors = 0, 0

    delimiter = tsv.get_delimiter(fragpipe_psm_file)
    reader = tsv.get_tsv_reader(fragpipe_psm_file, delimiter)
    headers = next(reader)

    # these columns will be overwritten
    protein_col = tsv.get_column_index(headers, "Protein")
    other_proteins_col = tsv.get_column_index(headers, "Mapped Proteins")

    fragpipe_psm_file_out = fragpipe_psm_file + ".tmp"
    if output_folder is not None:
        output_folder = f"{output_folder}/{Path(fragpipe_psm_file).parts[-2]}"
        Path(output_folder).mkdir(parents=True, exist_ok=True)
        fragpipe_psm_file_out = f"{output_folder}/{Path(fragpipe_psm_file).name}.tmp"

    writer = tsv.get_tsv_writer(fragpipe_psm_file_out)
    writer.writerow(headers)
    for row, proteins in fragpipe.parse_fragpipe_psm_file_for_peptide_remapping(
        reader, headers
    ):
        row_protein_groups = protein_groups.get_protein_groups(proteins)

        if len(row_protein_groups) == 0:
            logger.debug(
                f"Could not find any of the proteins {proteins} in proteinGroups.txt"
            )
            missing_peptides_in_protein_groups += 1
            continue

        if discard_shared_peptides and helpers.is_shared_peptide(row_protein_groups):
            shared_peptide_precursors += 1
            continue

        unique_peptide_precursors += 1

        leading_protein = row_protein_groups[0][0]
        if leading_protein not in proteins:
            peptides_not_mapping_to_leading_protein += 1
            continue

        row[protein_col] = leading_protein
        row[other_proteins_col] = ", ".join(
            [p for p in proteins if p != leading_protein]
        )
        writer.writerow(row)

    if missing_peptides_in_protein_groups > 0:
        logger.debug(
            f"Skipped {missing_peptides_in_protein_groups} precursors from "
            f"proteins not present in proteinGroups.txt file"
        )

    logger.info(
        f"Skipped {peptides_not_mapping_to_leading_protein} precursors not "
        f"mapping to leading protein in protein group."
    )

    logger.info(
        f"Found {unique_peptide_precursors} precursors from unique and "
        f"{shared_peptide_precursors} precursors from shared peptides"
    )

    if output_folder is None:
        os.rename(fragpipe_psm_file, fragpipe_psm_file.replace(".tsv", ".original.tsv"))
        os.rename(fragpipe_psm_file_out, fragpipe_psm_file)
    else:
        os.rename(fragpipe_psm_file_out, fragpipe_psm_file_out.replace(".tmp", ""))


if __name__ == "__main__":
    main(sys.argv[1:])
