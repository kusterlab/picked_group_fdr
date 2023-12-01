import sys
import os
import logging
from typing import Dict, List, Optional, Tuple

import numpy as np

from .update_fragpipe_results import (
    write_fragpipe_combined_protein_file,
)
from .. import __version__, __copyright__
from .. import helpers
from .. import quant
from .. import digest
from .. import protein_annotation
from ..parsers import sage
from ..scoring import ProteinScoringStrategy
from ..picked_group_fdr import ArgumentParserWithLogger
from ..parsers import maxquant, tsv
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
            precursorQuant = quant.PrecursorQuant(
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
    """Generate protein group results in FragPipe's combined_protein.tsv format.
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
                    if charge == -1:
                        protein_group_results[protein_group_idx].precursorQuants[
                            pq_idx
                        ].charge = charge
                else:
                    # match-between-runs hit
                    precursorQuant = quant.PrecursorQuant(
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
