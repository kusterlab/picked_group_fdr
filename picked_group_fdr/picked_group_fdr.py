from pathlib import Path
import sys
import os
import logging
import argparse
from timeit import default_timer as timer
from typing import List, Dict, Tuple, Union

import numpy as np

from . import __version__, __copyright__
from . import digest
from . import protein_annotation
from . import helpers
from . import entrapment
from . import methods
from . import fdr
from . import writers
from . import quantification
from .digestion_params import (
    add_digestion_arguments,
    get_digestion_params_list,
    DigestionParams,
)
from .parsers import psm
from .protein_groups import ProteinGroups
from .results import ProteinGroupResults
from .plotter import PlotterFactory
from .peptide_info import PeptideInfoList

# for type hints only
from .scoring_strategy import ProteinScoringStrategy
from .plotter import Plotter, NoPlotter

logger = logging.getLogger(__name__)

GREETER = f"PickedGroupFDR version {__version__}\n{__copyright__}"


def main(argv: List[str]) -> None:
    """Main function."""
    logger.info(GREETER)
    logger.info(
        f'Issued command: {os.path.basename(__file__)} {" ".join(map(str, argv))}'
    )

    args = parse_args(argv)
    run_picked_group_fdr(args)


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
        nargs="+",
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
        nargs="+",
        help="""Path to lfq.tsv file(s) produced by Sage. This enables
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
                "fragpipe". "auto": decide based on input file format, uses
                "maxquant" if no suitable format is known; "maxquant": proteinGroups.txt
                format; "fragpipe": combined_protein.tsv format.""",
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
        help="""Tab-separated file with mapping from peptides to proteins; alternative for 
                --fasta flag if digestion is time consuming.""",
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


def run_picked_group_fdr(args: argparse.Namespace) -> None:
    """Run PickedGroupFDR algorithm."""
    start = timer()

    # set seed for random shuffling of protein groups, see competition.do_competition()
    np.random.seed(1)

    protein_annotations, use_pseudo_genes = protein_annotation.get_protein_annotations(
        args.fasta, args.fasta_contains_decoys, args.gene_level
    )
    method_configs = methods.get_methods(args.methods, use_pseudo_genes)

    peptide_to_protein_maps = [None]
    if requires_peptide_to_protein_map(method_configs):
        peptide_to_protein_maps = get_peptide_to_protein_maps_from_args(
            args, use_pseudo_genes
        )

    plotter = PlotterFactory.get_plotter(args.figure_base_fn, args.plot_figures)
    for method_config in method_configs:
        run_method(
            args,
            method_config,
            peptide_to_protein_maps,
            protein_annotations,
            plotter,
            use_pseudo_genes,
            apply_filename_suffix=len(method_configs) > 1,
        )

    end = timer()
    logger.info(
        f"PickedGroupFDR execution took {'%.1f' % (end - start)} seconds wall clock time"
    )

    plotter.decorate_plots()
    plotter.save_plots()
    plotter.show()


def run_method(
    args: argparse.Namespace,
    method_config: methods.MethodConfig,
    peptide_to_protein_maps: List[digest.PeptideToProteinMap],
    protein_annotations: Dict[str, protein_annotation.ProteinAnnotation],
    plotter: Plotter,
    use_pseudo_genes: bool,
    apply_filename_suffix: bool,
) -> None:
    logger.info(
        f"Protein group level estimation method: {method_config.label} ({method_config.long_description()})"
    )

    evidence_files = method_config.score_type.get_evidence_file(args)
    if not evidence_files:
        logger.warning(
            (
                f'No evidence input file found, skipping method "{method_config.label}". Check '
                "if an appropriate method was specified by the --methods flag."
            )
        )
        return

    peptide_info_list = parse_evidence_files(
        evidence_files,
        peptide_to_protein_maps,
        method_config.score_type,
        args.suppress_missing_peptide_warning,
    )

    plotter.set_series_label_base(method_config.label)
    protein_group_results = get_protein_group_results(
        peptide_info_list,
        args.mq_protein_groups,
        method_config,
        plotter,
        args.keep_all_proteins,
    )

    protein_groups_writer = writers.MinimalProteinGroupsWriter(protein_annotations)
    post_err_probs = None
    if args.do_quant:
        (
            protein_group_results,
            protein_groups_writer,
            post_err_probs,
        ) = do_quantification(
            method_config.score_type,
            args,
            protein_group_results,
            protein_annotations,
            use_pseudo_genes,
            peptide_to_protein_maps,
            args.suppress_missing_peptide_warning,
        )

    protein_groups_writer.append_quant_columns(
        protein_group_results, post_err_probs, args.psm_fdr_cutoff
    )

    if args.protein_groups_out:
        protein_groups_out = get_output_filename(
            args.protein_groups_out, apply_filename_suffix, method_config
        )
        write_protein_groups(
            protein_groups_writer,
            protein_group_results,
            protein_groups_out,
        )


def requires_peptide_to_protein_map(method_configs: List[methods.MethodConfig]) -> bool:
    """Check if any configuration requires a peptide-to-protein map."""
    for method_config in method_configs:
        if (
            method_config.grouping_strategy.needs_peptide_to_protein_map()
            or method_config.score_type.remaps_peptides_to_proteins()
        ):
            return True
    return False


def get_peptide_to_protein_maps_from_args(
    args: argparse.Namespace,
    use_pseudo_genes: bool,
) -> List[digest.PeptideToProteinMap]:
    parse_id = digest.parse_until_first_space
    if args.gene_level and not use_pseudo_genes:
        parse_id = protein_annotation.parse_gene_name_func

    digestion_params_list = get_digestion_params_list(args)
    return get_peptide_to_protein_maps(
        args.fasta,
        args.peptide_protein_map,
        digestion_params_list,
        args.mq_protein_groups,
        parse_id=parse_id,
    )


def get_peptide_to_protein_maps(
    fasta_file: str,
    peptide_protein_map_files: List[str],
    digestion_params_list: List[DigestionParams],
    mq_protein_groups_file: str,
    **kwargs,
):
    peptide_to_protein_maps = list()
    if fasta_file:
        for digestion_params in digestion_params_list:
            peptide_to_protein_maps.append(
                digest.get_peptide_to_protein_map_from_params(
                    fasta_file, [digestion_params], **kwargs
                )
            )
            entrapment.mark_entrapment_proteins(
                peptide_to_protein_maps[-1], mq_protein_groups_file
            )
    elif peptide_protein_map_files:
        logger.info("Loading peptide to protein map")
        for peptide_protein_map_file in peptide_protein_map_files:
            peptide_to_protein_maps.append(
                digest.get_peptide_to_protein_map_from_file(
                    peptide_protein_map_file, use_hash_key=False
                )
            )
    else:
        raise ValueError(
            (
                "No fasta or peptide to protein mapping file detected, please"
                "specify either the --fasta or --peptide_protein_map flags."
            )
        )
    return peptide_to_protein_maps


def get_protein_group_results(
    peptide_info_list: PeptideInfoList,
    mq_protein_groups_file: str,
    method_config: methods.MethodConfig,
    plotter: Union[Plotter, NoPlotter],
    keep_all_proteins: bool,
) -> ProteinGroupResults:
    picked_strategy, score_type, grouping_strategy = (
        method_config.picked_strategy,
        method_config.score_type,
        method_config.grouping_strategy,
    )

    protein_groups = grouping_strategy.group_proteins(
        peptide_info_list, mq_protein_groups_file
    )

    # for razor peptide strategy
    score_type.set_peptide_counts_per_protein(peptide_info_list)

    for rescue_step in grouping_strategy.get_rescue_steps():
        if rescue_step:
            if not score_type.can_do_protein_group_rescue():
                raise NotImplementedError(
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

    plotter.set_series_label(
        method_config, rescue_step=rescue_step
    )
    plotter.plot_qval_calibration_and_performance(
        reported_qvals, observed_qvals, absent_ratio=1.0
    )

    return protein_group_results


def parse_evidence_files(
    evidence_files: List[str],
    peptide_to_protein_maps: List[digest.PeptideToProteinMap],
    score_type: ProteinScoringStrategy,
    suppress_missing_peptide_warning: bool,
) -> PeptideInfoList:
    """Returns best score per peptide"""
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
    use_pseudo_genes: bool,
    peptide_to_protein_maps: List[digest.PeptideToProteinMap],
    suppress_missing_peptide_warning: bool,
) -> Tuple[ProteinGroupResults, writers.ProteinGroupsWriter, List]:
    if not score_type.can_do_quantification():
        logger.warning(
            "Skipping quantification: need input file with precursor quantifications."
        )
        return protein_group_results

    logger.info("Preparing for quantification")

    score_origin = score_type.score_origin.long_description().lower()
    if args.output_format == "auto" and score_origin in ["maxquant", "fragpipe"]:
        args.output_format = score_origin

    parse_id = digest.parse_until_first_space
    if args.gene_level and not use_pseudo_genes:
        parse_id = protein_annotation.parse_gene_name_func

    protein_groups_writer = writers.get_protein_groups_output_writer(
        protein_group_results,
        args.output_format,
        args,
        protein_annotations,
        parse_id,
        peptide_to_protein_maps,
    )

    protein_groups = ProteinGroups.from_protein_group_results(protein_group_results)
    experimental_design = quantification.get_experimental_design(args)
    discard_shared_peptides = True
    protein_group_results, post_err_probs = score_type.get_quantification_parser()(
        score_type.get_evidence_file(args),
        score_type.get_quantification_file(args),
        protein_groups,
        protein_group_results,
        peptide_to_protein_maps,
        experimental_design,
        discard_shared_peptides,
        score_type=score_type,
        suppress_missing_peptide_warning=suppress_missing_peptide_warning,
    )

    return protein_group_results, protein_groups_writer, post_err_probs


def write_protein_groups(
    protein_groups_writer: writers.ProteinGroupsWriter,
    protein_group_results: ProteinGroupResults,
    protein_groups_out: str,
) -> None:
    Path(protein_groups_out).parent.mkdir(parents=True, exist_ok=True)

    protein_groups_writer.write(protein_group_results, protein_groups_out)
    logger.info(f"Protein group results have been written to: {protein_groups_out}")


def get_output_filename(
    protein_groups_out: str,
    apply_filename_suffix: bool,
    method_config: methods.MethodConfig,
):
    if not apply_filename_suffix:
        return protein_groups_out

    base, ext = os.path.splitext(protein_groups_out)
    label = method_config.label
    if label is None:
        label = method_config.short_description(rescue_step=True)
    else:
        label = label.lower().replace(" ", "_")
    return f"{base}_{label}{ext}"


if __name__ == "__main__":
    main(sys.argv[1:])
