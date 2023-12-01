from pathlib import Path
import sys
import os
import logging
from typing import Dict, List, Optional, Tuple

import numpy as np

from .. import __version__, __copyright__
from .. import helpers
from .. import digest
from .. import protein_annotation
from .. import writers
from ..picked_group_fdr import ArgumentParserWithLogger
from ..parsers import maxquant
from ..parsers import tsv
from ..parsers import fragpipe
from ..protein_annotation import ProteinAnnotation
from ..protein_groups import ProteinGroups
from ..results import ProteinGroupResults
from ..precursor_quant import PrecursorQuant

# hacky way to get package logger when running as module
logger = logging.getLogger(__package__ + "." + __file__)


def parse_args(argv):
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
        "--combined_ion",
        default=None,
        metavar="I",
        help="""Path to combined_ion.tsv produced by IonQuant/FragPipe. This enables
                quantification of protein groups by PickedGroupFDR in 
                combined_protein.tsv.""",
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

    args = parse_args(argv)

    protein_groups = ProteinGroups.from_mq_protein_groups_file(args.protein_groups)
    protein_groups.create_index()

    db = "target" if args.fasta_contains_decoys else "concat"
    protein_annotations = protein_annotation.get_protein_annotations_multiple(
        args.fasta, db=db, parse_id=digest.parse_until_first_space
    )
    protein_sequences = digest.get_protein_sequences(
        args.fasta, db=db, parse_id=digest.parse_until_first_space
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

    if args.combined_ion is not None:
        # create a fresh ProteinGroupResults object for combined_protein.tsv
        protein_group_results = maxquant.parse_mq_protein_groups_file(
            args.protein_groups
        )
        generate_fragpipe_combined_protein_file(
            args.fragpipe_psm,
            args.combined_ion,
            protein_groups,
            protein_group_results,
            protein_annotations,
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
    - Protein ID (P00167) - from fasta
    - Entry Name (CYB5_HUMAN) - from fasta
    - Gene (CYB5A) - from fasta
    - Protein Description (Cytochrome b5) - from fasta
    - Mapped Genes
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
    protein_id_col = tsv.get_column_index(headers, "Protein ID")
    entry_name_col = tsv.get_column_index(headers, "Entry Name")
    gene_col = tsv.get_column_index(headers, "Gene")
    description_col = tsv.get_column_index(headers, "Protein Description")
    mapped_genes_col = tsv.get_column_index(headers, "Mapped Genes")
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
        protein_annotation = protein_annotations.get(
            leading_protein,
            ProteinAnnotation(id=leading_protein, fasta_header=leading_protein),
        )
        row[protein_id_col] = protein_annotation.uniprot_id
        row[entry_name_col] = protein_annotation.entry_name
        row[gene_col] = protein_annotation.gene_name
        row[description_col] = protein_annotation.description

        other_genes = [
            protein_annotations.get(
                p,
                ProteinAnnotation(id=leading_protein, fasta_header=leading_protein),
            ).gene_name
            for p in proteins
            if p != leading_protein
        ]
        other_genes = [g for g in other_genes if g is not None and len(g) > 0]
        row[mapped_genes_col] = ", ".join(other_genes)
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
        if os.path.isfile(fragpipe_psm_file):
            os.rename(
                fragpipe_psm_file, fragpipe_psm_file.replace(".tsv", ".original.tsv")
            )
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
    psm_fdr_cutoff: float = 0.01,
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
    experiment = "1"
    protein_group_results.experiments.append(experiment)
    protein_group_results, post_err_probs = add_precursor_quants(
        fragpipe_psm_file,
        protein_group_results,
        protein_groups,
        experiment,
        discard_shared_peptides,
    )

    fragpipe_protein_file_out = get_fragpipe_protein_out_path(
        output_folder, fragpipe_psm_file
    )

    fragpipe_columns = writers.FragPipeSingleProteinWriter(
        protein_groups, protein_annotations, protein_sequences
    )

    protein_group_results = writers.append_quant_columns(
        protein_group_results, fragpipe_columns, post_err_probs, psm_fdr_cutoff
    )

    protein_group_results.write(
        fragpipe_protein_file_out,
        header_dict=fragpipe_columns.get_header_dict(),
        format_extra_columns=fragpipe_columns.get_extra_columns_formatter(),
    )


def get_fragpipe_protein_out_path(output_folder, fragpipe_psm_file):
    fragpipe_protein_file_out = fragpipe_psm_file.replace("psm.tsv", "protein.tsv")
    if output_folder is not None:
        output_folder = f"{output_folder}/{Path(fragpipe_psm_file).parts[-2]}"
        Path(output_folder).mkdir(parents=True, exist_ok=True)
        fragpipe_protein_file_out = Path(fragpipe_psm_file).name.replace(
            "psm.tsv", "protein.tsv"
        )
        fragpipe_protein_file_out = f"{output_folder}/{fragpipe_protein_file_out}"
    else:
        fragpipe_protein_file_backup = fragpipe_protein_file_out.replace(
            ".tsv", ".original.tsv"
        )
        if os.path.isfile(fragpipe_protein_file_out) and not os.path.isfile(
            fragpipe_protein_file_backup
        ):
            os.rename(
                fragpipe_protein_file_out,
                fragpipe_protein_file_backup,
            )
    return fragpipe_protein_file_out


def add_precursor_quants(
    fragpipe_psm_file: str,
    protein_group_results: ProteinGroupResults,
    protein_groups: ProteinGroups,
    experiment: str,
    discard_shared_peptides: bool,
):
    delimiter = tsv.get_delimiter(fragpipe_psm_file)
    reader = tsv.get_tsv_reader(fragpipe_psm_file, delimiter)
    headers = next(reader)

    post_err_probs = []
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

        if not helpers.is_decoy(proteins):
            post_err_probs.append((post_err_prob, "", experiment, peptide))

        for protein_group_idx in protein_group_idxs:
            precursor_quant = PrecursorQuant(
                peptide=peptide,
                charge=charge,
                experiment=experiment,
                fraction=-1,
                intensity=np.nan,
                post_err_prob=post_err_prob,
                tmt_intensities=None,  # TODO: add TMT support
                silac_intensities=None,  # TODO: add SILAC support
                evidence_id=-1,
                assigned_mods=assigned_mods,
                observed_mods=observed_mods,
            )
            protein_group_results[protein_group_idx].precursorQuants.append(
                precursor_quant
            )
    return protein_group_results, post_err_probs


def generate_fragpipe_combined_protein_file(
    fragpipe_psm_files: List[str],
    combined_ion_file: str,
    protein_groups: ProteinGroups,
    protein_group_results: ProteinGroupResults,
    protein_annotations: Dict[str, ProteinAnnotation],
    output_folder: Optional[str] = None,
    psm_fdr_cutoff: float = 0.01,
    discard_shared_peptides: bool = True,
):
    """Generate combined_protein.tsv file from psm.tsv and protein annotations.

    Args:
        fragpipe_psm_file (str): file in Fragpipe's psm.tsv format
        fasta_file (str): fasta file with all protein sequences
    """
    experiments = []
    post_err_probs_combined = []
    for fragpipe_psm_file in fragpipe_psm_files:
        experiment = Path(fragpipe_psm_file).parent.name
        experiments.append(experiment)
        protein_group_results, post_err_probs = add_precursor_quants(
            fragpipe_psm_file,
            protein_group_results,
            protein_groups,
            experiment,
            discard_shared_peptides,
        )
        post_err_probs_combined.extend(post_err_probs)

    if combined_ion_file is not None:
        protein_group_results = update_precursor_quants(
            protein_group_results,
            protein_groups,
            combined_ion_file,
            discard_shared_peptides,
        )

    if output_folder is None:
        output_folder = Path(fragpipe_psm_files[0]).parents[1]
    fragpipe_combined_protein_file_out = f"{output_folder}/combined_protein.tsv"

    write_fragpipe_combined_protein_file(
        fragpipe_combined_protein_file_out,
        protein_groups,
        protein_group_results,
        protein_annotations,
        experiments,
        post_err_probs,
        psm_fdr_cutoff,
    )


def write_fragpipe_combined_protein_file(
    protein_groups_out_file: str,
    protein_groups: ProteinGroups,
    protein_group_results: ProteinGroupResults,
    protein_annotations: Dict[str, ProteinAnnotation],
    experiments: List[str],
    post_err_probs: List,
    psm_fdr_cutoff: float = 0.01,
):
    """Write combined_protein.tsv file using FragPipe's format.

    An additional column "Protein group" is added that is not present in
    FragPipe's output format.

    https://fragpipe.nesvilab.org/docs/tutorial_fragpipe_outputs.html

    Output columns (LFQ):
    - Protein (sp|P00167|CYB5_HUMAN) - from fasta
    - Protein ID (P00167) - from fasta
    - Entry Name (CYB5_HUMAN) - from fasta
    - Gene (CYB5A) - from fasta
    - Protein Length - from fasta
    - Organism (Homo sapiens OX=9606) - from fasta
    - Protein Existence (1:Experimental evidence at protein level) - from fasta, but PE field only contains integer
    - Description (Cytochrome b5) - from fasta
    - Protein Probability (1.000) - update this with 1 - protein-level PEP
    - Top Peptide Probability (0.990)
    - Combined Total Peptides
    - Combined Spectral Count
    - Combined Unique Spectral Count
    - Combined Total Spectral Count
    - <Experiment> Spectral Count
    - <Experiment> Unique Spectral Count
    - <Experiment> Total Spectral Count
    - <Experiment> Intensity
    - <Experiment> MaxLFQ Intensity
    - Indistinguishable Proteins (sp|P0CE48|EFTU2_ECOLI)

    Args:
        fragpipe_psm_file (str): file in Fragpipe's psm.tsv format
        fasta_file (str): fasta file with all protein sequences
    """
    fragpipe_columns = writers.FragPipeCombinedProteinWriter(
        protein_groups, protein_annotations
    )

    protein_group_results.experiments = experiments
    protein_group_results = writers.append_quant_columns(
        protein_group_results, fragpipe_columns, post_err_probs, psm_fdr_cutoff
    )

    if os.path.isfile(protein_groups_out_file) and not os.path.isfile(
        protein_groups_out_file.replace(".tsv", ".original.tsv")
    ):
        os.rename(
            protein_groups_out_file,
            protein_groups_out_file.replace(".tsv", ".original.tsv"),
        )

    protein_group_results.write(
        protein_groups_out_file,
        header_dict=fragpipe_columns.get_header_dict(experiments),
        format_extra_columns=fragpipe_columns.get_extra_columns_formatter(),
    )


def update_precursor_quants(
    protein_group_results: ProteinGroupResults,
    protein_groups: ProteinGroups,
    combined_ion_file: str,
    discard_shared_peptides: bool,
):
    delimiter = tsv.get_delimiter(combined_ion_file)
    reader = tsv.get_tsv_reader(combined_ion_file, delimiter)
    headers = next(reader)

    for (
        peptide,
        charge,
        assigned_mods,
        proteins,
        intensities,
    ) in fragpipe.parse_fragpipe_combined_ion_file(reader, headers):
        protein_group_idxs = protein_groups.get_protein_group_idxs(proteins)

        if len(protein_group_idxs) == 0:
            logger.debug(
                f"Could not find any of the proteins {proteins} in proteinGroups.txt"
            )
            continue

        if discard_shared_peptides and helpers.is_shared_peptide(protein_group_idxs):
            continue

        for protein_group_idx in protein_group_idxs:
            precursors_to_update: Dict[str, Tuple[float, int]] = {}
            for pq_idx, pq in enumerate(
                protein_group_results[protein_group_idx].precursorQuants
            ):
                if (
                    pq.peptide == peptide
                    and pq.charge == charge
                    and pq.assigned_mods == assigned_mods
                ):
                    if (
                        pq.post_err_prob
                        < precursors_to_update.get(pq.experiment, (1.01, np.nan))[0]
                    ):
                        precursors_to_update[pq.experiment] = (pq.post_err_prob, pq_idx)

            for experiment, intensity in intensities:
                if intensity == 0.0:
                    continue

                if experiment in precursors_to_update:
                    pq_idx = precursors_to_update[experiment][1]
                    protein_group_results[protein_group_idx].precursorQuants[
                        pq_idx
                    ].intensity = intensity
                else:
                    # match-between-runs hit
                    precursor_quant = PrecursorQuant(
                        peptide=peptide,
                        charge=charge,
                        experiment=experiment,
                        fraction=-1,
                        intensity=intensity,
                        post_err_prob=np.nan,
                        tmt_intensities=None,  # TODO: add TMT support
                        silac_intensities=None,  # TODO: add SILAC support
                        evidence_id=-1,
                        assigned_mods=assigned_mods,
                        observed_mods="",
                    )
                    protein_group_results[protein_group_idx].precursorQuants.append(
                        precursor_quant
                    )
    return protein_group_results


if __name__ == "__main__":
    main(sys.argv[1:])
