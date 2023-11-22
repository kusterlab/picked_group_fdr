import collections
from pathlib import Path
import sys
import os
import logging
from typing import Dict, List, Optional

import numpy as np

from .. import __version__, __copyright__
from .. import helpers
from ..picked_group_fdr import ArgumentParserWithLogger
from .. import digest
from .. import protein_annotation
from ..parsers import maxquant, tsv, fragpipe
from ..protein_annotation import ProteinAnnotation
from ..protein_groups import ProteinGroups
from ..quant.precursor_quant import PrecursorQuant
from ..quant.base import ProteinGroupColumns
from ..quant.fragpipe_protein_annotations import (
    FragpipeProteinAnnotationsColumns,
)
from ..quant.peptide_count import UniquePeptideCountColumns
from ..quant.sequence_coverage import SequenceCoverageColumns
from ..quant.protein_probability import ProteinProbabilityColumns
from ..quant.top_peptide import TopPeptideProbabilityColumns
from ..quant.spectral_count import SpectralCountColumns
from ..quant.sum_and_ibaq import SummedIntensityAndIbaqColumns
from ..quant.modifications import ModificationsColumns
from ..quant.indistinguishable_proteins import IndistinguishableProteinsColumns
from ..results import ProteinGroupResults


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

    protein_annotations = protein_annotation.get_protein_annotations_multiple(
        args.fasta, parseId=digest.parseUntilFirstSpace
    )
    protein_sequences = digest.getProteinSequences(
        args.fasta, parseId=digest.parseUntilFirstSpace
    )

    for fragpipe_psm_file in args.fragpipe_psm:
        fragpipe_psm_file_out = update_fragpipe_psm_file_single(
            fragpipe_psm_file, protein_groups, protein_annotations, args.output_folder
        )

        # create a fresh ProteinGroupResults object for each psm.tsv
        protein_group_results = maxquant.parse_mq_protein_groups_file(
            args.protein_groups
        )
        generate_fragpipe_protein_file(
            fragpipe_psm_file_out,
            protein_groups,
            protein_group_results,
            protein_annotations,
            protein_sequences,
            args.output_folder,
        )


def update_fragpipe_psm_file_single(
    fragpipe_psm_file: str,
    protein_groups: ProteinGroups,
    protein_annotations: Dict[str, ProteinAnnotation],
    output_folder: Optional[str] = None,
    discard_shared_peptides: bool = True,
) -> str:
    """Update protein mappings for each peptide using the PickedGroupFDR protein groups.

    These columns are updated:
    - Protein
    - Protein ID (P00167) - from fasta (TODO)
    - Entry Name (CYB5_HUMAN) - from fasta (TODO)
    - Gene (CYB5A) - from fasta (TODO)
    - Protein Description (Cytochrome b5) - from fasta (TODO)
    - Mapped Genes (TODO)
    - Mapped Proteins

    Args:
        fragpipe_psm_file (str): _description_
        protein_groups (ProteinGroups): _description_
        output_folder (Union[str, None], optional): _description_. Defaults to None.
        discard_shared_peptides (bool, optional): _description_. Defaults to True.

    Returns:
        str: path to updated psm.tsv file
    """
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

    missing_peptides_in_protein_groups = 0
    peptides_not_mapping_to_leading_protein = 0
    shared_peptide_precursors, unique_peptide_precursors = 0, 0

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

    return fragpipe_psm_file_out.replace(".tmp", "")


def generate_fragpipe_protein_file(
    fragpipe_psm_file: str,
    protein_groups: ProteinGroups,
    protein_group_results: ProteinGroupResults,
    protein_annotations: Dict[str, ProteinAnnotation],
    protein_sequences: Dict[str, str],
    output_folder: Optional[str] = None,
    discard_shared_peptides: bool = True,
):
    """Generate experiment specific protein.tsv file from psm.tsv and fasta file.

    https://fragpipe.nesvilab.org/docs/tutorial_fragpipe_outputs.html

    Output columns (LFQ):
    - Protein (sp|P00167|CYB5_HUMAN) - from fasta
    - Protein ID (P00167) - from fasta
    - Entry Name (CYB5_HUMAN) - from fasta
    - Gene (CYB5A) - from fasta
    - Length - from fasta
    - Organism (Homo sapiens OX=9606) - from fasta
    - Protein Description (Cytochrome b5) - from fasta
    - Protein Existence (1:Experimental evidence at protein level) - from fasta, but PE field only contains integer
    - Coverage
    - Protein Probability (1.000) - update this with 1 - protein-level PEP
    - Top Peptide Probability (0.990)
    - Total Peptides
    - Unique Peptides
    - Razor Peptides
    - Total Spectral Count
    - Unique Spectral Count
    - Razor Spectral Count
    - Total Intensity (0)
    - Unique Intensity (0)
    - Razor Intensity (0)
    - Razor Assigned Modifications (19M(15.9949),19M(15.9949),22C(57.0215)) - recurring mod_peptide+charge pairs are not included
    - Razor Observed Modifications - seems to always be empty
    - Indistinguishable Proteins (sp|P0CE48|EFTU2_ECOLI)

    Args:
        fragpipe_psm_file (str): file in Fragpipe's psm.tsv format
        fasta_file (str): fasta file with all protein sequences
    """
    delimiter = tsv.get_delimiter(fragpipe_psm_file)
    reader = tsv.get_tsv_reader(fragpipe_psm_file, delimiter)
    headers = next(reader)

    for (
        peptide,
        charge,
        post_err_prob,
        assigned_mods,
        observed_mods,
        proteins,
    ) in fragpipe.parse_fragpipe_psm_file_for_protein_tsv(reader, headers):
        protein_group_idxs = protein_groups.get_protein_group_idxs(proteins)

        if len(protein_group_idxs) == 0:
            logger.debug(
                f"Could not find any of the proteins {proteins} in proteinGroups.txt"
            )
            continue

        if discard_shared_peptides and helpers.is_shared_peptide(protein_group_idxs):
            continue

        for proteinGroupIdx in protein_group_idxs:
            precursorQuant = PrecursorQuant(
                peptide=peptide,
                charge=charge,
                experiment="1",
                fraction=-1,
                intensity=np.nan,
                post_err_prob=post_err_prob,
                tmt_intensities=None,  # TODO: add TMT support
                silac_intensities=None,  # TODO: add SILAC support
                evidence_id=-1,
                assigned_mods=assigned_mods,
                observed_mods=observed_mods,
            )
            protein_group_results[proteinGroupIdx].precursorQuants.append(
                precursorQuant
            )
    
    protein_group_results.remove_protein_groups_without_precursors()

    for header in maxquant.MQ_PROTEIN_ANNOTATION_HEADERS:
        protein_group_results.remove_column(header)

    silac_channels = []
    num_ibaq_peptides_per_protein = collections.defaultdict(lambda: 1)

    columns: List[ProteinGroupColumns] = [
        FragpipeProteinAnnotationsColumns(protein_groups, protein_annotations),
        SequenceCoverageColumns(protein_sequences),
        ProteinProbabilityColumns(),
        TopPeptideProbabilityColumns(),
        UniquePeptideCountColumns(),
        SpectralCountColumns(),
        SummedIntensityAndIbaqColumns(silac_channels, num_ibaq_peptides_per_protein),
        ModificationsColumns(),
        IndistinguishableProteinsColumns(),
    ]

    experiments = ["1"]
    experiment_to_idx_map = {"1": 0}
    post_err_prob_cutoff = 1.01
    for c in columns:
        c.append_headers(protein_group_results, experiments)
        c.append_columns(
            protein_group_results, experiment_to_idx_map, post_err_prob_cutoff
        )

    fragpipe_protein_file_out = fragpipe_psm_file.replace("psm.tsv", "protein.tsv")
    if output_folder is not None:
        output_folder = f"{output_folder}/{Path(fragpipe_psm_file).parts[-2]}"
        Path(output_folder).mkdir(parents=True, exist_ok=True)
        fragpipe_protein_file_out = Path(fragpipe_psm_file).name.replace(
            "psm.tsv", "protein.tsv"
        )
        fragpipe_protein_file_out = f"{output_folder}/{fragpipe_protein_file_out}"
    else:
        os.rename(
            fragpipe_protein_file_out,
            fragpipe_protein_file_out.replace(".tsv", ".original.tsv"),
        )

    header_dict = {
        "Protein": "Protein",
        "Protein ID": "Protein ID",
        "Entry Name": "Entry Name",
        "Gene": "Gene",
        "Length": "Length",
        "Organism": "Organism",
        "Protein Description": "Protein Description",
        "Protein Existence": "Protein Existence",
        "Coverage": "Sequence coverage [%]",
        "Protein Probability": "Protein Probability",
        "Top Peptide Probability": "Top Peptide Probability",
        "Total Peptides": "Unique peptides 1",
        "Unique Peptides": "Unique peptides 1",
        "Razor Peptides": "Unique peptides 1",
        "Total Spectral Count": "Spectral count 1",
        "Unique Spectral Count": "Spectral count 1",
        "Razor Spectral Count": "Spectral count 1",
        "Total Intensity": "Intensity 1",
        "Unique Intensity": "Intensity 1",
        "Razor Intensity": "Intensity 1",
        "Razor Assigned Modifications": "Razor Assigned Modifications",
        "Razor Observed Modifications": "Razor Observed Modifications",
        "Indistinguishable Proteins": "Indistinguishable Proteins"
    }
    protein_group_results.write(fragpipe_protein_file_out, header_dict=header_dict)


if __name__ == "__main__":
    main(sys.argv[1:])
