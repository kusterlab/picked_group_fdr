from pathlib import Path
import sys
import os
import logging
import argparse
from timeit import default_timer as timer
from typing import Any, Callable, List, Dict, Union

import numpy as np

from picked_group_fdr.protein_groups import ProteinGroups

from . import __version__, __copyright__
from . import digest
from . import protein_annotation
from . import helpers
from . import proteotypicity
from . import entrapment
from . import methods
from . import fdr
from . import results
from . import writers
from . import quantification
from .quant import maxquant as mq_quant
from .digestion_params import add_digestion_arguments, get_digestion_params_list
from .parsers import psm
from .grouping import PseudoGeneGrouping
from .results import ProteinGroupResults
from .plotter import PlotterFactory
from .peptide_info import PeptideInfoList

# for type hints only
from .scoring_strategy import ProteinScoringStrategy
from .grouping import ProteinGroupingStrategy
from .competition import ProteinCompetitionStrategy
from .plotter import Plotter, NoPlotter

logger = logging.getLogger(__name__)

GREETER = f"PickedGroupFDR version {__version__}\n{__copyright__}"


class ArgumentParserWithLogger(argparse.ArgumentParser):
    def error(self, message):
        logger.error(f"Error parsing input arguments: {message}")
        super().error(message)


def parse_args(argv):
    apars = ArgumentParserWithLogger(
        description=GREETER, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    apars.add_argument(
        "--mq_evidence",
        default=None,
        metavar="EV",
        nargs="+",
        help="""MaxQuant evidence file(s).""",
    )

    apars.add_argument(
        "--perc_evidence",
        default=None,
        metavar="POUT",
        nargs="+",
        help="""Percolator output file(s) with PSMs or peptides; alternative for 
                --mq_evidence.""",
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
        help="""Path to combined_ion.tsv produced by IonQuant/FragPipe. This enables
                quantification of protein groups by PickedGroupFDR.""",
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
        "--sage_lfq_tsv",
        default=None,
        metavar="I",
        help="""Path to lfq.tsv produced by Sage. This enables
                quantification of protein groups by PickedGroupFDR.""",
    )

    apars.add_argument(
        "--protein_groups_out",
        default=None,
        metavar="PG",
        help="""Protein groups output file, mimicks a subset of the MQ protein groups columns.""",
    )

    apars.add_argument(
        "--output_format",
        default="auto",
        metavar="PG",
        help="""Protein groups output format. Options are "auto", "maxquant" and 
                "fragpipe". "auto": decide based on input file format, will default to 
                "maxquant" if no suitable format is known; "maxquant" (proteinGroups.txt
                format) and "fragpipe" (combined_protein.tsv format).""",
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
        "--mq_protein_groups",
        default=None,
        metavar="PG",
        help="""MaxQuant protein groups file; only specify if we want to keep 
                MaxQuant's original grouping instead of Picked Grouping.""",
    )

    apars.add_argument(
        "--methods",
        default="picked_protein_group",
        metavar="M1,M2",
        help="""Use one or more protein group FDR estimation methods, separated by 
                commas. Examples of builtin methods: picked_protein_group, 
                picked_protein_group_mq_input, savitski, maxquant. Alternatively, 
                specify one or more paths to toml files with a .toml extension 
                following the same format as the builtin method toml files.""",
    )

    apars.add_argument(
        "--peptide_protein_map",
        default=None,
        metavar="M",
        nargs="+",
        help="""File with mapping from peptides to proteins; alternative for 
                --fasta flag if digestion is time consuming.""",
    )

    apars.add_argument(
        "--peptide_proteotypicity_map",
        default=None,
        metavar="M",
        help="""File with mapping from peptides to proteotypicity.""",
    )

    apars.add_argument(
        "--keep_all_proteins",
        default=None,
        help="""Keep proteins that do not have peptides below the PSM FDR filter.""",
        action="store_true",
    )

    apars.add_argument(
        "--gene_level",
        help="""Report gene-level statistics instead of protein group-level. This 
                requires the GN= field to be present in the fasta file.""",
        action="store_true",
    )

    apars.add_argument(
        "--do_quant",
        help="""Do protein quantification, will calculate summed intensity, iBAQ 
                and LFQ intensities.""",
        action="store_true",
    )

    apars.add_argument(
        "--suppress_missing_peptide_warning",
        help="Suppress missing peptide warning when mapping peptides to proteins.",
        action="store_true",
    )

    apars.add_argument(
        "--figure_base_fn",
        default=None,
        metavar="F",
        help="""Base file name for calibration and performance figures.""",
    )

    apars.add_argument(
        "--plot_figures", help="Plot figures with matplotlib.", action="store_true"
    )

    add_digestion_arguments(apars)

    quantification.add_quant_arguments(apars)

    # ------------------------------------------------
    args = apars.parse_args(argv)

    return args


def write_protein_groups(
    protein_groups_writer: writers.ProteinGroupsWriter,
    protein_group_results: results.ProteinGroupResults,
    protein_groups_out: str,
    config: Dict[str, Any],
    apply_suffix: bool,
):
    if apply_suffix:
        base, ext = os.path.splitext(protein_groups_out)
        label = config.get("label", None)
        if label is None:
            label = methods.short_description(
                config["scoreType"], config["grouping"], config["pickedStrategy"], True
            )
        else:
            label = label.lower().replace(" ", "_")
        protein_groups_out = f"{base}_{label}{ext}"
    
    Path(protein_groups_out).parent.mkdir(parents=True, exist_ok=True)

    protein_groups_writer.write(protein_group_results, protein_groups_out)
    logger.info(f"Protein group results have been written to: {protein_groups_out}")


def main(argv):
    logger.info(GREETER)
    logger.info(
        f'Issued command: {os.path.basename(__file__)} {" ".join(map(str, argv))}'
    )

    args = parse_args(argv)

    # set seed for random shuffling of protein groups, see competition.do_competition()
    np.random.seed(1)

    configs = methods.get_methods(args)
    digestion_params_list = get_digestion_params_list(args)

    start = timer()

    plotter = PlotterFactory.get_plotter(args.figure_base_fn, args.plot_figures)
    plotter.init_plots()

    parse_id = digest.parse_until_first_space
    db = "target" if args.fasta_contains_decoys else "concat"
    protein_annotations = dict()
    if args.fasta:
        protein_annotations = protein_annotation.get_protein_annotations_multiple(
            args.fasta, db=db, parse_id=parse_id
        )
        if args.gene_level:
            if protein_annotation.has_gene_names(
                protein_annotations, min_ratio_with_genes=0.5
            ):
                parse_id = protein_annotation.parse_gene_name_func
                protein_annotations = (
                    protein_annotation.get_protein_annotations_multiple(
                        args.fasta, db=db, parse_id=parse_id
                    )
                )
            else:
                logger.warning(
                    (
                        "Found fewer than 50% of proteins without gene names in the "
                        "fasta file, will attempt to infer pseudo-genes based on "
                        "shared peptides instead."
                    )
                )
                for config in configs:
                    config["grouping"] = PseudoGeneGrouping()

    peptide_to_protein_maps = list()
    peptide_to_proteotypicity_map = dict()
    for config in configs:
        method_description_long = methods.long_description(
            config["scoreType"], config["grouping"], config["pickedStrategy"], True
        )
        label = config.get("label", "")
        logger.info(
            f"Protein group level estimation method: {label} ({method_description_long})"
        )

        evidence_files = config["scoreType"].get_evidence_file(args)
        if not evidence_files:
            logger.warning(
                (
                    f'No evidence input file found, skipping method "{label}". Check '
                    "if an appropriate method was specified by the --methods flag."
                )
            )
            continue

        if len(peptide_to_protein_maps) == 0 and (
            config["grouping"].needs_peptide_to_protein_map()
            or config["scoreType"].remaps_peptides_to_proteins()
        ):
            peptide_to_protein_maps = get_peptide_to_protein_maps(
                args.fasta,
                args.peptide_protein_map,
                digestion_params_list,
                args.mq_protein_groups,
            )

        if (
            len(peptide_to_proteotypicity_map) == 0
            and config["scoreType"].use_proteotypicity
        ):
            peptide_to_proteotypicity_map = (
                proteotypicity.get_peptide_to_proteotypicity_from_file(
                    args.peptide_proteotypicity_map
                )
            )

        peptide_info_list = parse_evidence_files(
            evidence_files,
            peptide_to_protein_maps,
            config["scoreType"],
            args.suppress_missing_peptide_warning,
        )

        plotter.set_series_label_base(config.get("label", None))
        protein_group_results = get_protein_group_results(
            peptide_info_list,
            args.mq_protein_groups,
            peptide_to_proteotypicity_map,
            config["pickedStrategy"],
            config["scoreType"],
            config["grouping"],
            plotter,
            args.keep_all_proteins,
        )

        protein_groups_writer = writers.MinimalProteinGroupsWriter(protein_annotations)

        post_err_probs = None
        psm_fdr_cutoff = 0.01
        if args.do_quant:
            (
                protein_group_results,
                protein_groups_writer,
                post_err_probs,
            ) = do_quantification(
                config["scoreType"],
                args,
                protein_group_results,
                protein_annotations,
                parse_id,
                peptide_to_protein_maps,
            )

        protein_groups_writer.append_quant_columns(
            protein_group_results, post_err_probs, psm_fdr_cutoff
        )

        if args.protein_groups_out:
            write_protein_groups(
                protein_groups_writer,
                protein_group_results,
                args.protein_groups_out,
                config,
                apply_suffix=len(configs) > 1,
            )

    end = timer()
    logger.info(
        f"PickedGroupFDR execution took {'%.1f' % (end - start)} seconds wall clock time"
    )

    plotter.decorate_plots()
    plotter.save_plots()
    plotter.show()


def get_peptide_to_protein_maps(
    fasta, peptide_protein_map, digestion_params_list, mq_protein_groups
):
    peptide_to_protein_maps = list()
    if fasta:
        for digestion_params in digestion_params_list:
            peptide_to_protein_maps.append(
                digest.get_peptide_to_protein_map_from_params(fasta, [digestion_params])
            )
            entrapment.mark_entrapment_proteins(
                peptide_to_protein_maps[-1], mq_protein_groups
            )
    elif peptide_protein_map:
        logger.info("Loading peptide to protein map")
        for peptide_protein_map in peptide_protein_map:
            peptide_to_protein_maps.append(
                digest.get_peptide_to_protein_map_from_file(
                    peptide_protein_map, use_hash_key=False
                )
            )
    else:
        sys.exit(
            (
                "No fasta or peptide to protein mapping file detected, please"
                "specify either the --fasta or --peptide_protein_map flags."
            )
        )
    return peptide_to_protein_maps


def get_protein_group_results(
    peptide_info_list: PeptideInfoList,
    mq_protein_groups_file: str,
    peptide_to_proteotypicity_map,
    picked_strategy: ProteinCompetitionStrategy,
    score_type: ProteinScoringStrategy,
    grouping_strategy: ProteinGroupingStrategy,
    plotter: Union[Plotter, NoPlotter],
    keep_all_proteins: bool,
):
    protein_groups = grouping_strategy.group_proteins(
        peptide_info_list, mq_protein_groups_file
    )

    # for razor peptide strategy
    score_type.set_peptide_counts_per_protein(peptide_info_list)

    for rescue_step in grouping_strategy.get_rescue_steps():
        if rescue_step:
            if not score_type.can_do_protein_group_rescue():
                raise Exception(
                    "Cannot do rescue step for other score types than bestPEP"
                )
            protein_groups = grouping_strategy.rescue_protein_groups(
                peptide_info_list,
                protein_group_results,
                protein_groups,
                protein_group_peptide_infos,
            )

        protein_group_peptide_infos = score_type.collect_peptide_scores_per_protein(
            protein_groups,
            peptide_info_list,
            suppress_missing_protein_warning=rescue_step,
        )

        # find optimal division factor for multPEP score
        score_type.optimize_hyperparameters(protein_groups, protein_group_peptide_infos)

        if rescue_step and picked_strategy.short_description() == "pgT":
            (
                protein_groups,
                protein_group_peptide_infos,
            ) = grouping_strategy.update_protein_groups(
                protein_groups, protein_group_peptide_infos
            )

        (
            picked_protein_groups,
            picked_protein_group_peptide_infos,
            protein_scores,
        ) = picked_strategy.do_competition(
            protein_groups, protein_group_peptide_infos, score_type
        )

        reported_qvals, observed_qvals = fdr.calculate_protein_fdrs(
            picked_protein_groups, protein_scores
        )

        peptide_score_cutoff = (
            score_type.peptide_score_cutoff if rescue_step else float("inf")
        )
        protein_group_results = ProteinGroupResults.from_protein_groups(
            picked_protein_groups,
            picked_protein_group_peptide_infos,
            protein_scores,
            reported_qvals,
            peptide_score_cutoff,
            keep_all_proteins,
        )

    if score_type.use_proteotypicity:
        proteotypicity.calculate_proteotypicity_scores(
            picked_protein_groups,
            picked_protein_group_peptide_infos,
            peptide_to_proteotypicity_map,
            score_type,
            peptide_score_cutoff,
        )

    plotter.set_series_label(
        score_type, grouping_strategy, picked_strategy, rescue_step=rescue_step
    )
    plotter.plot_qval_calibration_and_performance(
        reported_qvals, observed_qvals, absent_ratio=1.0
    )

    return protein_group_results


def parse_evidence_files(
    evidence_files: List[str],
    peptide_to_protein_maps: List[Dict[str, List[str]]],
    score_type: ProteinScoringStrategy,
    suppress_missing_peptide_warning: bool,
) -> PeptideInfoList:
    """Returns best score per peptide"""
    if not score_type.remaps_peptides_to_proteins():
        peptide_to_protein_maps = [None]

    if len(peptide_to_protein_maps) == 1:
        peptide_to_protein_maps = peptide_to_protein_maps * len(evidence_files)

    peptide_info_list = dict()
    for peptide, proteins, _, score in psm.parse_evidence_file_multiple(
        evidence_files,
        peptide_to_protein_maps=peptide_to_protein_maps,
        score_type=score_type,
        suppress_missing_peptide_warning=suppress_missing_peptide_warning,
    ):
        peptide = helpers.clean_peptide(peptide)
        if np.isnan(score) or score >= peptide_info_list.get(peptide, [np.inf])[0]:
            continue

        peptide_info_list[peptide] = [score, proteins]

    return peptide_info_list


def do_quantification(
    score_type: ProteinScoringStrategy,
    args: argparse.Namespace,
    protein_group_results: ProteinGroupResults,
    protein_annotations: Dict[str, protein_annotation.ProteinAnnotation],
    parse_id: Callable,
    peptide_to_protein_maps: List[Dict[str, List[str]]],
):
    if not score_type.can_do_quantification():
        logger.warning(
            "Skipping quantification: need input file with precursor quantifications."
        )
        return protein_group_results

    logger.info("Preparing for quantification")

    score_origin = score_type.score_origin.long_description().lower()
    if args.output_format == "auto" and score_origin in ["maxquant", "fragpipe"]:
        args.output_format = score_origin

    protein_groups_writer = writers.get_protein_groups_output_writer(
        protein_group_results,
        args.output_format,
        args,
        protein_annotations,
        parse_id,
        peptide_to_protein_maps,
    )

    discard_shared_peptides = True

    protein_groups = ProteinGroups.from_protein_group_results(protein_group_results)
    protein_groups.create_index()
    
    protein_group_results, post_err_probs = score_type.get_quantification_parser()(
        score_type.get_evidence_file(args),
        score_type.get_quantification_file(args),
        protein_groups,
        protein_group_results,
        peptide_to_protein_maps,
        args.file_list_file,
        discard_shared_peptides,
    )

    return protein_group_results, protein_groups_writer, post_err_probs


if __name__ == "__main__":
    main(sys.argv[1:])
