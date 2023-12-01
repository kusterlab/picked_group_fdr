import sys
import os
import logging
from typing import Dict, List, Tuple

import numpy as np

from . import serializers
from . import digest
from . import digestion_params
from . import protein_annotation
from . import helpers
from . import quant
from . import picked_group_fdr
from .parsers import maxquant
from .parsers import psm
from .parsers import parsers
from .quant import protein_annotations
from .protein_groups import ProteinGroups
from .scoring import ProteinScoringStrategy

# for type hints only
from .results import ProteinGroupResults

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
        help="""Fasta file to create mapping from peptides to proteins. This should not
                contain the decoy sequences, unless you set the --fasta_contains_decoys 
                flag.""",
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

    apars.add_argument(
        "--file_list_file",
        metavar="L",
        help="""Tab separated file with lines of the format (third and 
                fourth columns are optional): raw_file <tab> condition 
                <tab> experiment <tab> fraction.
                """,
        required=False,
    )

    apars.add_argument(
        "--lfq_min_peptide_ratios",
        default=2,
        type=int,
        metavar="M",
        help="""Minimum number of common peptides between two samples
                to qualify for calculating a peptide ratio in LFQ
                """,
    )

    apars.add_argument(
        "--num_threads",
        default=1,
        type=int,
        metavar="T",
        help="""Maximum number of threads to use.""",
    )

    apars.add_argument(
        "--lfq_stabilize_large_ratios",
        help="""Apply stabilization of large ratios in LFQ as described
                in the MaxLFQ paper.""",
        action="store_false",
    )

    digestion_params.add_digestion_arguments(apars)

    # ------------------------------------------------
    args = apars.parse_args(argv)

    return args


def main(argv):
    logger.info(
        f'Issued command: {os.path.basename(__file__)} {" ".join(map(str, argv))}'
    )

    args = parse_args(argv)

    parse_id = digest.parse_until_first_space
    if args.gene_level:
        parse_id = protein_annotation.parse_gene_name_func
    elif args.fasta_use_uniprot_id:
        parse_id = protein_annotation.parse_uniprot_id

    (
        peptide_to_protein_maps,
        num_ibaq_peptides_per_protein,
    ) = get_peptide_to_protein_maps(args)
    protein_sequences = digest.get_protein_sequences(args.fasta, parse_id)

    protein_group_results = maxquant.parse_mq_protein_groups_file(
        args.mq_protein_groups,
        additional_headers=protein_annotations.MQ_PROTEIN_ANNOTATION_HEADERS,
    )

    protein_group_results = do_quantification(
        args.mq_evidence,
        protein_group_results,
        protein_sequences,
        peptide_to_protein_maps,
        num_ibaq_peptides_per_protein,
        args.file_list_file,
        min_peptide_ratios_lfq=args.lfq_min_peptide_ratios,
        stabilize_large_ratios_lfq=args.lfq_stabilize_large_ratios,
        num_threads=args.num_threads,
    )

    protein_group_results.write(args.protein_groups_out)

    logger.info(
        f"Protein group results have been written to: {args.protein_groups_out}"
    )


def get_peptide_to_protein_maps(args):
    logger.info("Loading peptide to protein map...")

    digestion_params_list = digestion_params.get_digestion_params_list(args)

    peptide_to_protein_maps = picked_group_fdr.get_peptide_to_protein_maps(
        args.fasta,
        args.peptide_protein_map,
        digestion_params_list,
        args.mq_protein_groups,
    )
    if args.fasta:
        logger.info("In silico protein digest for iBAQ")
        num_ibaq_peptides_per_protein = digest.get_num_ibaq_peptides_per_protein(
            args.fasta, digestion_params_list
        )
    elif args.peptide_protein_map:
        logger.warning("Found peptide_protein_map (instead of fasta input): ")
        logger.warning(
            "- calculating iBAQ values using all peptides in peptide_protein_map."
        )
        logger.warning("- cannot compute sequence coverage.")
        num_ibaq_peptides_per_protein = digest.get_num_peptides_per_protein(
            digest.merge_peptide_to_protein_maps(peptide_to_protein_maps)
        )
    else:
        raise ValueError(
            "No peptide to protein map found, use either the --fasta or the --peptide_protein_map arguments"
        )

    return peptide_to_protein_maps, num_ibaq_peptides_per_protein


def do_quantification(
    mq_evidence_files: List[str],
    protein_group_results: ProteinGroupResults,
    protein_sequences: Dict[str, str],
    peptide_to_protein_maps: List[Dict[str, List[str]]],
    num_ibaq_peptides_per_protein: Dict[str, int],
    file_list_file: str,
    psm_fdr_cutoff: float = 0.01,
    discard_shared_peptides: bool = True,
    min_peptide_ratios_lfq: int = 2,
    stabilize_large_ratios_lfq: bool = True,
    num_threads: int = 1,
):
    params = init_triqler_params()

    logger.info("Preparing for quantification")
    file_mapping = None
    if file_list_file:
        (
            protein_group_results.experiments,
            file_mapping,
            params,
        ) = parsers.parse_file_list(file_list_file, params)

    protein_group_results, post_err_probs = parse_evidence_files(
        protein_group_results,
        mq_evidence_files,
        peptide_to_protein_maps,
        file_mapping,
        discard_shared_peptides,
    )

    columns = serializers.get_mq_protein_groups_columns(
        num_ibaq_peptides_per_protein,
        protein_sequences,
        min_peptide_ratios_lfq,
        stabilize_large_ratios_lfq,
        num_threads,
        params,
    )

    return serializers.append_quant_columns(
        protein_group_results, columns, post_err_probs, psm_fdr_cutoff
    )


def parse_evidence_files(
    protein_group_results: ProteinGroupResults,
    mq_evidence_files: List[str],
    peptide_to_protein_maps: List[Dict[str, List[str]]],
    file_mapping: Dict[str, Tuple[str, str]],
    discard_shared_peptides: bool,
):
    protein_groups = ProteinGroups.from_protein_group_results(protein_group_results)
    protein_groups.create_index()

    post_err_probs = list()
    shared_peptide_precursors, unique_peptide_precursors = 0, 0
    parsed_experiments = set()
    missing_peptides_in_protein_groups = 0

    for (
        peptide,
        proteins,
        charge,
        raw_file,
        experiment,
        fraction,
        intensity,
        post_err_prob,
        tmt_cols,
        silac_cols,
        evidence_id,
    ) in psm.parse_evidence_file_multiple(
        mq_evidence_files,
        peptide_to_protein_maps=peptide_to_protein_maps,
        score_type=ProteinScoringStrategy("bestPEP"),
        for_quantification=True,
    ):
        if protein_group_results.num_tmt_channels is None:
            # There are 3 columns per TMT channel:
            #     Reporter intensity corrected,
            #     Reporter intensity
            #     Reporter intensity count
            protein_group_results.num_tmt_channels = int(len(tmt_cols) / 3)
        if protein_group_results.num_silac_channels is None:
            protein_group_results.num_silac_channels = len(silac_cols)

        # override the parsed experiment and fraction if --file_list_file option is used
        if file_mapping:
            experiment, fraction = file_mapping[raw_file]
        elif experiment not in parsed_experiments:
            parsed_experiments.add(experiment)

        protein_group_idxs = protein_groups.get_protein_group_idxs(proteins)

        # removes peptides not present in the proteinGroups.txt file
        if len(protein_group_idxs) == 0:
            logger.debug(
                f"Could not find any of the proteins {proteins} in proteinGroups.txt"
            )
            missing_peptides_in_protein_groups += 1
            continue

        if discard_shared_peptides and helpers.is_shared_peptide(protein_group_idxs):
            shared_peptide_precursors += 1
            continue

        unique_peptide_precursors += 1

        if not helpers.is_decoy(proteins):
            post_err_probs.append((post_err_prob, raw_file, experiment, peptide))

        if len(tmt_cols) > 0:
            tmt_cols = np.array(tmt_cols, dtype="float64")
        else:
            tmt_cols = None

        if len(silac_cols) > 0:
            silac_cols = np.array(silac_cols, dtype="float64")
        else:
            silac_cols = None

        for protein_group_idx in protein_group_idxs:
            precursor_quant = quant.PrecursorQuant(
                peptide,
                charge,
                experiment,
                fraction,
                intensity,
                post_err_prob,
                tmt_cols,
                silac_cols,
                evidence_id,
            )
            protein_group_results[protein_group_idx].precursorQuants.append(
                precursor_quant
            )

    if missing_peptides_in_protein_groups > 0:
        logger.debug(
            f"Skipped {missing_peptides_in_protein_groups} precursors from proteins not present in proteinGroups.txt file"
        )

    logger.info(
        f"Found {unique_peptide_precursors} precursors from unique and {shared_peptide_precursors} precursors from shared peptides"
    )

    if len(parsed_experiments) > 0:
        protein_group_results.experiments = sorted(list(parsed_experiments))

    return protein_group_results, post_err_probs


def init_triqler_params():
    params = dict()
    # TODO: make these parameters configurable from the command line
    params["decoyPattern"] = "REV__"
    params["groups"] = []
    params["groupLabels"] = []
    params["numThreads"] = 4
    params["warningFilter"] = "ignore"
    params["foldChangeEval"] = 0.8
    params["returnPosteriors"] = False
    params["minSamples"] = 5
    return params


if __name__ == "__main__":
    main(sys.argv[1:])
