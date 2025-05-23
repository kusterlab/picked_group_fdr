import sys
import os
import logging
import argparse
from timeit import default_timer as timer
from typing import List, Dict, Optional, Union

import numpy as np

from . import __version__, __copyright__
from . import digest
from . import protein_annotation
from . import peptide_protein_map
from . import methods
from . import fdr
from .parsers import evidence
from . import writers
from . import quantification
from .digestion_params import (
    add_digestion_arguments,
)
from .results import ProteinGroupResults
from .plotter import PlotterFactory

# for type hints only
from .plotter import Plotter, NoPlotter
from .protein_groups import ProteinGroups
from .peptide_info import PeptideInfoList, ProteinGroupPeptideInfos

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
        "--diann_reports",
        default=None,
        metavar="R",
        nargs="+",
        help="""DIA-NN report file(s) in .parquet or .tsv format.""",
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
        "--protein_group_fdr_threshold",
        default=0.01,
        metavar="C",
        type=float,
        help="""Protein group-level FDR threshold used for setting the peptide level 
                cutoff in the rescuing procedure and computing the number of
                protein groups at the given threshold. Note that this does not filter the
                protein_groups_out file for the specified FDR. Use the
                picked_group_fdr.pipeline.filter_fdr_maxquant module for this.""",
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
        args.fasta,
        args.fasta_contains_decoys,
        args.gene_level,
        args.fasta_use_uniprot_id,
    )
    method_configs = methods.get_methods(args.methods, use_pseudo_genes)

    peptide_to_protein_maps = [None]
    if methods.requires_peptide_to_protein_map(method_configs):
        peptide_to_protein_maps = (
            peptide_protein_map.get_peptide_to_protein_maps_from_args(
                args, use_pseudo_genes
            )
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

    peptide_info_list = evidence.parse_evidence_files(
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
        args.protein_group_fdr_threshold,
        args.psm_fdr_cutoff,
    )

    protein_groups_writer = writers.MinimalProteinGroupsWriter(protein_annotations)
    post_err_probs = None
    if args.do_quant and method_config.score_type.can_do_quantification():
        (
            protein_group_results,
            protein_groups_writer,
            post_err_probs,
        ) = quantification.do_quantification(
            method_config.score_type,
            args,
            protein_group_results,
            protein_annotations,
            use_pseudo_genes,
            peptide_to_protein_maps,
            args.suppress_missing_peptide_warning,
        )

    writers.finalize_output(
        protein_group_results,
        protein_groups_writer,
        post_err_probs,
        args.protein_groups_out,
        args.psm_fdr_cutoff,
        apply_filename_suffix,
        method_config,
    )


def get_protein_group_results(
    peptide_info_list: PeptideInfoList,
    mq_protein_groups_file: str = "",
    method_config: Union[methods.MethodConfig, None] = None,
    plotter: Union[Plotter, NoPlotter, None] = None,
    keep_all_proteins: bool = False,
    protein_group_fdr_threshold: float = 0.01,
    psm_fdr_cutoff: float = 0.01,
) -> ProteinGroupResults:
    """Performs protein inference and protein group-level FDR estimation.

    Args:
        peptide_info_list (PeptideInfoList = Dict[str, Tuple[float, List[str]]]):
            Dictionary of (peptide_sequence) => (peptide_search_engine_score,
            [protein1, ...]).
        mq_protein_groups_file (str, optional): Path to protein groups file in
            MaxQuant proteinGroups.txt format, only used if
            method_config.grouping_strategy is MQNativeGrouping. Defaults to "".
        method_config (Union[methods.MethodConfig, None], optional):
            methods.MethodConfig object that specifies the protein grouping,
            scoring and picked strategies to use. If set to `None`, uses
            `RescuedSubsetGrouping`, `BestPEPScore` and `PickedGroupStrategy`.
            Defaults to None.
        plotter (Union[Plotter, NoPlotter, None], optional): Plotter object for
            plotting FDR calibration curves if given entrapment search results.
            If set to `None`, `NoPlotter` (i.e. no plots are generated) is used.
            Defaults to None.
        keep_all_proteins (bool, optional): Keep proteins that do not have
            peptides below `psm_fdr_cutoff`. Defaults to False.
        protein_group_fdr_threshold (float, optional): Protein group-level FDR
            threshold to aim for in rescuing procedure. Note that this does not
            filter out protein groups above this threshold. Defaults to 0.01.
        psm_fdr_cutoff (float, optional): PSM-level FDR used solely for counting
            the number of peptides for each protein group. Defaults to 0.01.

    Raises:
        NotImplementedError: Raised if an invalid combination of
            method.MethodConfig is passed.

    Returns:
        ProteinGroupResults: `ProteinGroupResults` object containing a list
            of `ProteinGroupResult` objects with the proteins in the group with
            associated information, such as protein group-level FDR and number
            of peptides.
    """
    if method_config is None:
        method_config = methods.parse_method_toml(
            method_name="picked_protein_group", use_pseudo_genes=False
        )

    if plotter is None:
        plotter = NoPlotter()

    _, score_type, grouping_strategy = (
        method_config.picked_strategy,
        method_config.score_type,
        method_config.grouping_strategy,
    )

    protein_groups = grouping_strategy.group_proteins(
        peptide_info_list, mq_protein_groups_file
    )

    # for razor peptide strategy
    score_type.set_peptide_counts_per_protein(peptide_info_list)

    if grouping_strategy.has_rescue_step():
        method_config.score_type.save_protein_score_strategy()
        method_config.score_type.set_protein_score_strategy("bestPEP")

    (
        protein_group_results,
        protein_groups,
        protein_group_peptide_infos,
        reported_qvals,
        observed_qvals,
    ) = get_protein_group_results_single_step(
        peptide_info_list,
        method_config,
        protein_groups,
        keep_all_proteins=keep_all_proteins,
        protein_group_fdr_threshold=protein_group_fdr_threshold,
        psm_fdr_cutoff=psm_fdr_cutoff,
    )

    if grouping_strategy.has_rescue_step():
        method_config.score_type.reset_protein_score_strategy()
        (
            protein_group_results,
            _,
            _,
            reported_qvals,
            observed_qvals,
        ) = get_protein_group_results_single_step(
            peptide_info_list,
            method_config,
            protein_groups,
            protein_group_peptide_infos=protein_group_peptide_infos,
            protein_group_results=protein_group_results,
            rescue_step=True,
            keep_all_proteins=keep_all_proteins,
            protein_group_fdr_threshold=protein_group_fdr_threshold,
            psm_fdr_cutoff=psm_fdr_cutoff,
        )

    plotter.set_series_label(
        method_config, rescue_step=grouping_strategy.has_rescue_step()
    )
    plotter.plot_qval_calibration_and_performance(
        reported_qvals, observed_qvals, absent_ratio=1.0
    )

    return protein_group_results


def get_protein_group_results_single_step(
    peptide_info_list: PeptideInfoList,
    method_config: methods.MethodConfig,
    protein_groups: ProteinGroups,
    protein_group_peptide_infos: Optional[ProteinGroupPeptideInfos] = None,
    protein_group_results: Optional[ProteinGroupResults] = None,
    rescue_step: bool = False,
    keep_all_proteins: bool = False,
    protein_group_fdr_threshold: float = 0.01,
    psm_fdr_cutoff: float = 0.01,
) -> ProteinGroupResults:
    if protein_group_peptide_infos is None:
        protein_group_peptide_infos = []

    if protein_group_results is None:
        protein_group_results = ProteinGroupResults()

    picked_strategy, score_type, grouping_strategy = (
        method_config.picked_strategy,
        method_config.score_type,
        method_config.grouping_strategy,
    )

    obsolete_protein_groups, obsolete_protein_group_peptide_infos = [], []
    if rescue_step:
        if not score_type.can_do_protein_group_rescue():
            raise NotImplementedError(
                "Cannot do rescue step for other score types than bestPEP"
            )
        (
            protein_groups,
            obsolete_protein_groups,
            obsolete_protein_group_peptide_infos,
        ) = grouping_strategy.rescue_protein_groups(
            peptide_info_list,
            protein_group_results,
            protein_group_fdr_threshold,
            protein_groups,
            protein_group_peptide_infos,
        )

    protein_group_peptide_infos = score_type.collect_peptide_scores_per_protein(
        protein_groups,
        peptide_info_list,
        psm_fdr_cutoff,
        suppress_missing_protein_warning=rescue_step,
    )

    # find optimal division factor for multPEP score
    score_type.optimize_hyperparameters(
        protein_groups, protein_group_peptide_infos, protein_group_fdr_threshold
    )

    (
        picked_protein_groups,
        picked_protein_group_peptide_infos,
        protein_scores,
    ) = picked_strategy.do_competition(
        protein_groups,
        protein_group_peptide_infos,
        obsolete_protein_groups,
        obsolete_protein_group_peptide_infos,
        score_type,
    )

    reported_qvals, observed_qvals = fdr.calculate_protein_fdrs(
        picked_protein_groups, protein_scores, protein_group_fdr_threshold
    )

    # peptide-level score cutoff for counting number of peptides per protein
    # this cutoff is computed above in score_type.collect_peptide_scores_per_protein()
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

    return (
        protein_group_results,
        protein_groups,
        protein_group_peptide_infos,
        reported_qvals,
        observed_qvals,
    )


if __name__ == "__main__":
    main(sys.argv[1:])
