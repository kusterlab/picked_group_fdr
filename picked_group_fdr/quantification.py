import sys
import os
import logging
import collections
from typing import Dict, List, Tuple

import numpy as np
import triqler.parsers

from . import digest
from . import digestion_params
from . import protein_annotation
from . import helpers
from . import quant
from . import fdr
from .parsers import maxquant
from .parsers import psm
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

    peptide_to_protein_map, num_ibaq_peptides_per_protein = get_peptide_to_protein_maps(
        args, parse_id
    )
    protein_sequences = digest.get_protein_sequences(args.fasta, parse_id)
    protein_group_results = maxquant.parse_mq_protein_groups_file(
        args.mq_protein_groups,
        additional_headers=maxquant.MQ_PROTEIN_ANNOTATION_HEADERS,
    )

    do_quantification(
        args.mq_evidence,
        protein_group_results,
        protein_sequences,
        peptide_to_protein_map,
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


def get_peptide_to_protein_maps(args, parse_id):
    min_len_ibaq = max([6, args.min_length])
    max_len_ibaq = min([30, args.max_length])

    logger.info("Loading peptide to protein map...")
    if args.fasta:
        pre, not_post, post = digest.get_cleavage_sites(args.enzyme)

        db = "concat"
        if args.fasta_contains_decoys:
            db = "target"

        peptide_to_protein_map = digest.get_peptide_to_protein_map(
            args.fasta,
            db=db,
            digestion=args.digestion,
            min_len=args.min_length,
            max_len=args.max_length,
            pre=pre,
            not_post=not_post,
            post=post,
            miscleavages=args.cleavages,
            methionine_cleavage=True,
            special_aas=list(args.special_aas),
            parse_od=parse_id,
            use_hash_key=(args.digestion == "none"),
        )

        peptide_to_protein_map_ibaq = digest.get_peptide_to_protein_map(
            args.fasta,
            db=db,
            digestion=args.digestion,
            min_len=min_len_ibaq,
            max_len=max_len_ibaq,
            pre=pre,
            not_post=not_post,
            post=post,
            miscleavages=0,
            methionine_cleavage=False,
            special_aas=list(args.special_aas),
            parse_od=parse_id,
            use_hash_key=(args.digestion == "none"),
        )
    elif args.peptide_protein_map:
        pre, not_post, post = digest.get_cleavage_sites(args.enzyme)

        peptide_to_protein_map = digest.get_peptide_to_protein_map_from_file(
            args.peptide_protein_map, useHashKey=True
        )

        peptide_to_protein_map_ibaq = dict()
        for peptide, proteins in peptide_to_protein_map.items():
            peptide_len = len(peptide)
            if (
                peptide_len >= min_len_ibaq
                and peptide_len <= max_len_ibaq
                and not digest.has_miscleavage(peptide, pre, not_post, post)
            ):
                peptide_to_protein_map_ibaq[peptide] = proteins
    else:
        raise ValueError(
            "No peptide to protein map found, use either the --fasta or the --peptide_protein_map arguments"
        )

    num_ibaq_peptides_per_protein = digest.get_num_peptides_per_protein(
        peptide_to_protein_map_ibaq
    )

    return peptide_to_protein_map, num_ibaq_peptides_per_protein


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
        experiments, file_mapping, params = parse_file_list(file_list_file, params)

    (
        protein_group_results,
        post_err_probs,
        num_tmt_channels,
        num_silac_channels,
        parsed_experiments,
    ) = parse_evidence_files(
        protein_group_results,
        mq_evidence_files,
        peptide_to_protein_maps,
        file_mapping,
        discard_shared_peptides,
    )

    silac_channels = get_silac_channels(num_silac_channels)

    if len(parsed_experiments) > 0:
        experiments = sorted(list(parsed_experiments))
    experiment_to_idx_map = dict([(v, k) for k, v in enumerate(experiments)])

    # (1) technically this is a precursor-level FDR and not a PSM-level FDR
    # (2) in contrast to MaxQuant, we set a global precursor-level FDR
    #         instead of a per raw file PSM-level FDR
    post_err_prob_cutoff = fdr.calc_post_err_prob_cutoff(
        [x[0] for x in post_err_probs if not helpers.is_mbr(x[0])], psm_fdr_cutoff
    )
    logger.info(
        f"PEP-cutoff corresponding to {psm_fdr_cutoff*100:g}% PSM-level FDR: {post_err_prob_cutoff}"
    )

    print_num_peptides_at_fdr(post_err_probs, post_err_prob_cutoff)

    logger.info("Filtering for identified precursors")
    # precursor = (peptide, charge) tuple
    # this filter also ensures that MBR precursors which were matched to
    # unidentified precursors are removed
    for pgr in protein_group_results:
        pgr.precursorQuants = retain_only_identified_precursors(
            pgr.precursorQuants, post_err_prob_cutoff
        )

    columns: List[quant.ProteinGroupColumns] = [
        quant.UniquePeptideCountColumns(),
        quant.IdentificationTypeColumns(),
        quant.SummedIntensityAndIbaqColumns(silac_channels, num_ibaq_peptides_per_protein),
        quant.SequenceCoverageColumns(protein_sequences),
        quant.EvidenceIdsColumns(),
    ]

    if num_tmt_channels > 0:
        columns.append(quant.TMTIntensityColumns(num_tmt_channels))
    else:
        columns.append(
            quant.LFQIntensityColumns(
                silac_channels,
                min_peptide_ratios_lfq,
                stabilize_large_ratios_lfq,
                num_threads,
            )
        )
        # TODO: add SILAC functionality of Triqler
        if num_silac_channels == 0:
            columns.append(quant.TriqlerIntensityColumns(params))

    for c in columns:
        c.append_headers(protein_group_results, experiments)
        c.append_columns(protein_group_results, experiment_to_idx_map, post_err_prob_cutoff)


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
    num_tmt_channels, num_silac_channels = -1, -1
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
        if num_tmt_channels == -1:
            # There are 3 columns per TMT channel:
            #     Reporter intensity corrected,
            #     Reporter intensity
            #     Reporter intensity count
            num_tmt_channels = int(len(tmt_cols) / 3)
        if num_silac_channels == -1:
            num_silac_channels = len(silac_cols)

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

    return (
        protein_group_results,
        post_err_probs,
        num_tmt_channels,
        num_silac_channels,
        parsed_experiments,
    )


def print_num_peptides_at_fdr(post_err_probs: List, post_err_prob_cutoff: float):
    surviving_mod_peptides = set(
        [x[3] for x in post_err_probs if x[0] <= post_err_prob_cutoff]
    )

    peptides_per_rawfile = collections.defaultdict(list)
    peptides_per_experiment = collections.defaultdict(list)
    peptides_per_rawfile_mbr = collections.defaultdict(list)
    peptides_per_experiment_mbr = collections.defaultdict(list)
    for post_err_prob, rawfile, experiment, peptide in post_err_probs:
        if post_err_prob <= post_err_prob_cutoff:
            peptides_per_rawfile[rawfile].append(peptide)
            peptides_per_experiment[experiment].append(peptide)
        elif helpers.is_mbr(post_err_prob) and peptide in surviving_mod_peptides:
            peptides_per_rawfile_mbr[rawfile].append(peptide)
            peptides_per_experiment_mbr[experiment].append(peptide)

    logger.info("Precursor counts per rawfile (1% PSM-level FDR):")
    for rawfile, peptides in sorted(peptides_per_rawfile.items()):
        num_peptides = len(set(peptides))
        num_peptides_with_mbr = len(set(peptides + peptides_per_rawfile_mbr[rawfile]))
        logger.info(
            f"    {rawfile}: {num_peptides} {'(' + str(num_peptides_with_mbr) + ' with MBR)' if num_peptides_with_mbr > num_peptides else ''}"
        )

    logger.info("Precursor counts per experiment (1% PSM-level FDR):")
    for experiment, peptides in sorted(peptides_per_experiment.items()):
        num_peptides = len(set(peptides))
        num_peptides_with_mbr = len(set(peptides + peptides_per_experiment_mbr[experiment]))
        logger.info(
            f"    {experiment}: {num_peptides} {'(' + str(num_peptides_with_mbr) + ' with MBR)' if num_peptides_with_mbr > num_peptides else ''}"
        )


def get_silac_channels(num_silac_channels: int):
    silac_channels = list()
    if num_silac_channels == 3:
        silac_channels = ["L", "M", "H"]
    elif num_silac_channels == 2:
        silac_channels = ["L", "H"]
    elif num_silac_channels != 0:
        sys.exit("ERROR: Found a number of SILAC channels not equal to 2 or 3")
    return silac_channels


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


def parse_file_list(file_list_file: str, params: Dict):
    file_info_list = triqler.parsers.parseFileList(file_list_file)
    file_mapping = dict()
    experiments = list()
    for rawfile, condition, experiment, fraction in file_info_list:
        if experiment not in experiments:
            experiments.append(experiment)
        file_mapping[rawfile] = (experiment, fraction)
        # Note that params["groupLabels"] and params["groups"] are only used by Triqler
        if condition not in params["groupLabels"]:
            params["groupLabels"].append(condition)
            params["groups"].append([])
        params["groups"][params["groupLabels"].index(condition)].append(
            experiments.index(experiment)
        )
    return experiments, file_mapping, params


def retain_only_identified_precursors(
    precursor_list: List[quant.PrecursorQuant], post_err_prob_cutoff
):
    identified_precursors = set()
    for precursor in precursor_list:
        if precursor.post_err_prob <= post_err_prob_cutoff:
            identified_precursors.add((precursor.peptide, precursor.charge))
    return [
        precursor_row
        for precursor_row in precursor_list
        if (precursor_row.peptide, precursor_row.charge) in identified_precursors
    ]


if __name__ == "__main__":
    main(sys.argv[1:])
