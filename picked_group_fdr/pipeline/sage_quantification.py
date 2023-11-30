import collections
from pathlib import Path
import sys
import os
import logging
from typing import Dict, List, Optional, Tuple

import numpy as np

from .. import __version__, __copyright__
from .. import helpers
from .. import fdr
from ..parsers import sage
from .update_fragpipe_results import (
    get_fragpipe_combined_protein_headers,
    fragpipe_format_extra_columns,
)
from ..scoring import ProteinScoringStrategy
from ..quantification import retain_only_identified_precursors
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
from ..quant.protein_probability import ProteinProbabilityColumns
from ..quant.top_peptide import TopPeptideProbabilityColumns
from ..quant.spectral_count import SpectralCountColumns
from ..quant.sum_and_ibaq import SummedIntensityAndIbaqColumns
from ..quant.lfq import LFQIntensityColumns
from ..quant.indistinguishable_proteins import IndistinguishableProteinsColumns
from ..results import ProteinGroupResults


# hacky way to get the package logger instead of just __main__ when running as python -m picked_group_fdr.pipeline.update_evidence_from_pout ...
logger = logging.getLogger(__package__ + "." + __file__)

FRAGPIPE_PROTEIN_OUTPUT_DICT = {
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
    "Indistinguishable Proteins": "Indistinguishable Proteins",
}

FRAGPIPE_COMBINED_PROTEIN_OUTPUT_DICT = {
    "Protein": "Protein",
    "Protein ID": "Protein ID",
    "Entry Name": "Entry Name",
    "Gene": "Gene",
    "Protein Length": "Length",
    "Organism": "Organism",
    "Protein Existence": "Protein Existence",
    "Description": "Protein Description",
    "Combined Total Peptides": "Combined Total Peptides",
    "Combined Spectral Count": "Combined Spectral Count",
    "Combined Unique Spectral Count": "Combined Spectral Count",
    "Combined Total Spectral Count": "Combined Spectral Count",
}


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

    protein_annotations = protein_annotation.get_protein_annotations_multiple(
        args.fasta, parse_id=digest.parse_until_first_space
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


def add_precursor_quants(
    fragpipe_psm_file: str,
    protein_group_results: ProteinGroupResults,
    protein_groups: ProteinGroups,
    discard_shared_peptides: bool,
):
    delimiter = tsv.get_delimiter(fragpipe_psm_file)
    reader = tsv.get_tsv_reader(fragpipe_psm_file, delimiter)
    headers = next(reader)

    get_proteins = lambda peptide, proteins: proteins
    score_type = ProteinScoringStrategy("Sage bestPEP")

    post_err_probs = []
    for (
        peptide,
        proteins,
        experiment,
        post_err_prob,
        charge,
    ) in sage.parse_sage_results_file(
        reader, headers, get_proteins, score_type, for_quantification=True
    ):
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
            precursorQuant = PrecursorQuant(
                peptide=peptide,
                charge=charge,
                experiment=experiment,
                fraction=-1,
                intensity=np.nan,
                post_err_prob=post_err_prob,
                tmt_intensities=None,  # TODO: add TMT support
                silac_intensities=None,  # TODO: add SILAC support
                evidence_id=-1,
            )
            protein_group_results[protein_group_idx].precursorQuants.append(
                precursorQuant
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
    """Generate experiment specific protein.tsv file from psm.tsv and fasta file.

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
    post_err_probs_combined = []
    for fragpipe_psm_file in fragpipe_psm_files:
        protein_group_results, post_err_probs = add_precursor_quants(
            fragpipe_psm_file,
            protein_group_results,
            protein_groups,
            discard_shared_peptides,
        )
        post_err_probs_combined.extend(post_err_probs)

    protein_group_results, experiments = update_precursor_quants(
        protein_group_results,
        protein_groups,
        combined_ion_file,
        discard_shared_peptides,
    )

    protein_group_results.remove_protein_groups_without_precursors()

    silac_channels = []
    num_ibaq_peptides_per_protein = collections.defaultdict(lambda: 1)

    columns: List[ProteinGroupColumns] = [
        FragpipeProteinAnnotationsColumns(protein_groups, protein_annotations),
        ProteinProbabilityColumns(),
        TopPeptideProbabilityColumns(),
        UniquePeptideCountColumns(),
        SpectralCountColumns(),
        SummedIntensityAndIbaqColumns(silac_channels, num_ibaq_peptides_per_protein),
        IndistinguishableProteinsColumns(),
    ]

    min_peptide_ratios_lfq = 1
    stabilize_large_ratios_lfq = True
    num_threads = 1
    columns.append(
        LFQIntensityColumns(
            silac_channels,
            min_peptide_ratios_lfq,
            stabilize_large_ratios_lfq,
            num_threads,
        )
    )

    experiment_to_idx_map = {
        experiment: idx for idx, experiment in enumerate(experiments)
    }
    post_err_prob_cutoff = fdr.calc_post_err_prob_cutoff(
        [x[0] for x in post_err_probs if not helpers.is_mbr(x[0])], psm_fdr_cutoff
    )
    logger.info(
        f"PEP-cutoff corresponding to {psm_fdr_cutoff*100:g}% PSM-level FDR: {post_err_prob_cutoff}"
    )

    logger.info("Filtering for identified precursors")
    # precursor = (peptide, charge) tuple
    # this filter also ensures that MBR precursors which were matched to
    # unidentified precursors are removed
    for pgr in protein_group_results:
        pgr.precursorQuants = retain_only_identified_precursors(
            pgr.precursorQuants, post_err_prob_cutoff
        )

    for c in columns:
        c.append_headers(protein_group_results, experiments)
        c.append_columns(
            protein_group_results, experiment_to_idx_map, post_err_prob_cutoff
        )

    if output_folder is not None:
        fragpipe_combined_protein_file_out = f"{output_folder}/combined_protein.tsv"
    else:
        fragpipe_combined_protein_file_out = (
            f"{Path(fragpipe_psm_files[0]).parents[1]}/combined_protein.tsv"
        )
        if os.path.isfile(fragpipe_combined_protein_file_out):
            os.rename(
                fragpipe_combined_protein_file_out,
                fragpipe_combined_protein_file_out.replace(".tsv", ".original.tsv"),
            )

    protein_group_results.write(
        fragpipe_combined_protein_file_out,
        header_dict=get_fragpipe_combined_protein_headers(experiments),
        format_extra_columns=fragpipe_format_extra_columns,
    )


def update_precursor_quants(
    protein_group_results: ProteinGroupResults,
    protein_groups: ProteinGroups,
    sage_lfq_tsv: str,
    discard_shared_peptides: bool,
):
    delimiter = tsv.get_delimiter(sage_lfq_tsv)
    reader = tsv.get_tsv_reader(sage_lfq_tsv, delimiter)
    headers = next(reader)

    experiments = sage.get_experiments_from_sage_lfq_headers(headers)

    for (
        peptide,
        charge,
        proteins,
        intensities,
    ) in sage.parse_sage_lfq_file(reader, headers):
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
                # if quant.lfq_settings.combine_charge_state is set to true, charge is
                # always -1
                if pq.peptide == peptide and (pq.charge == charge or charge == -1):
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
                    protein_group_results[protein_group_idx].precursorQuants[
                        pq_idx
                    ].charge = charge
                else:
                    # match-between-runs hit
                    precursorQuant = PrecursorQuant(
                        peptide=peptide,
                        charge=charge,
                        experiment=experiment,
                        fraction=-1,
                        intensity=intensity,
                        post_err_prob=np.nan,
                        tmt_intensities=None,  # TODO: add TMT support
                        silac_intensities=None,  # TODO: add SILAC support
                        evidence_id=-1,
                    )
                    protein_group_results[protein_group_idx].precursorQuants.append(
                        precursorQuant
                    )
    return protein_group_results, experiments


if __name__ == "__main__":
    main(sys.argv[1:])
