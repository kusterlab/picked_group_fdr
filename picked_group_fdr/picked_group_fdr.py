import sys
import os
import logging
import argparse
from timeit import default_timer as timer
from typing import List, Dict, Union

import numpy as np

from . import __version__, __copyright__
from . import digest
from . import helpers
from . import parsers
from . import proteotypicity
from . import quantification
from . import entrapment
from . import methods
from . import fdr

from .digestion_params import add_digestion_arguments, get_digestion_params_list
from .grouping import PseudoGeneGrouping
from .results import ProteinGroupResults
from .plotter import PlotterFactory
from .peptide_info import PeptideInfoList

# for type hints only
from .scoring import ProteinScoringStrategy
from .grouping import ProteinGroupingStrategy
from .competition import ProteinCompetitionStrategy
from .plotter import Plotter, NoPlotter

logger = logging.getLogger(__name__)

GREETER = f"PickedGroupFDR version {__version__}\n{__copyright__}"


class ArgumentParserWithLogger(argparse.ArgumentParser):
    def error(self, message):
        logger.error(f"Error parsing input arguments: {message}")
        super().error(message)


def parseArgs(argv):
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
        "--protein_groups_out",
        default=None,
        metavar="PG",
        help="""Protein groups output file, mimicks a subset of the MQ protein groups columns.""",
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
        "--perc_evidence",
        default=None,
        metavar="POUT",
        nargs="+",
        help="""Percolator output file(s) with PSMs or peptides; alternative for 
                --mq_evidence if we want to use Percolator PEPs instead of 
                MaxQuant's PEPs.""",
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
        "--lfq_min_peptide_ratios",
        default=2,
        type=int,
        metavar="M",
        help="""Minimum number of common peptides between two samples
                to qualify for calculating a peptide ratio in LFQ.""",
    )

    apars.add_argument(
        "--num_threads",
        default=1,
        type=int,
        metavar="T",
        help="""Maximum number of threads to use. Currently only speeds up the MaxLFQ part.""",
    )

    apars.add_argument(
        "--file_list_file",
        metavar="L",
        help="""Tab separated file with lines of the format (third and fourth 
                columns are optional): raw_file <tab> condition <tab> experiment 
                <tab> fraction.""",
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

    # ------------------------------------------------
    args = apars.parse_args(argv)

    return args


def main(argv):
    logger.info(GREETER)
    logger.info(
        f'Issued command: {os.path.basename(__file__)} {" ".join(map(str, argv))}'
    )

    args = parseArgs(argv)

    # set seed for random shuffling of protein groups, see competition.do_competition()
    np.random.seed(1)

    configs = methods.get_methods(args)
    digestion_params_list = get_digestion_params_list(args)

    start = timer()

    plotter = PlotterFactory.get_plotter(args.figure_base_fn, args.plot_figures)
    plotter.initPlots()

    parseId = digest.parseUntilFirstSpace
    proteinAnnotations = dict()
    if args.fasta:
        proteinAnnotations = digest.get_protein_annotations_multiple(
            args.fasta, parseId
        )
        if args.gene_level:
            if digest.hasGeneNames(proteinAnnotations, minRatioWithGenes=0.5):
                parseId = digest.parseGeneNameFunc
                proteinAnnotations = digest.get_protein_annotations_multiple(
                    args.fasta, parseId
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

    peptideToProteinMaps = list()
    peptideToProteotypicityMap = dict()
    for config in configs:
        methodDescriptionLong = methods.long_description(
            config["scoreType"], config["grouping"], config["pickedStrategy"], True
        )
        label = config.get("label", "")
        logger.info(
            f"Protein group level estimation method: {label} ({methodDescriptionLong})"
        )

        evidenceFiles = config["scoreType"].get_evidence_file(args)
        if not evidenceFiles:
            logger.warning(
                (
                    f'No evidence input file found, skipping method "{label}". Check '
                    "if an appropriate method was specified by the --methods flag."
                )
            )
            continue

        if len(peptideToProteinMaps) == 0 and (
            config["grouping"].needs_peptide_to_protein_map()
            or config["scoreType"].remaps_peptides_to_proteins()
        ):
            if args.fasta:
                for digestion_params in digestion_params_list:
                    peptideToProteinMaps.append(
                        digest.get_peptide_to_protein_map_from_params(
                            args.fasta, [digestion_params]
                        )
                    )
                    entrapment.markEntrapmentProteins(
                        peptideToProteinMaps[-1], args.mq_protein_groups
                    )
            elif args.peptide_protein_map:
                logger.info("Loading peptide to protein map")
                for peptide_protein_map in args.peptide_protein_map:
                    peptideToProteinMaps.append(
                        digest.getPeptideToProteinMapFromFile(
                            peptide_protein_map, useHashKey=False
                        )
                    )
            else:
                sys.exit(
                    (
                        "No fasta or peptide to protein mapping file detected, please"
                        "specify either the --fasta or --peptide_protein_map flags."
                    )
                )

        if (
            len(peptideToProteotypicityMap) == 0
            and config["scoreType"].use_proteotypicity
        ):
            peptideToProteotypicityMap = (
                proteotypicity.getPeptideToProteotypicityFromFile(
                    args.peptide_proteotypicity_map
                )
            )

        peptideInfoList = parseEvidenceFiles(
            evidenceFiles,
            peptideToProteinMaps,
            config["scoreType"],
            args.suppress_missing_peptide_warning,
        )

        plotter.set_series_label_base(config.get("label", None))
        proteinGroupResults = getProteinGroupResults(
            peptideInfoList,
            args.mq_protein_groups,
            proteinAnnotations,
            peptideToProteotypicityMap,
            config["pickedStrategy"],
            config["scoreType"],
            config["grouping"],
            plotter,
            args.keep_all_proteins,
        )

        if args.do_quant:
            doQuantification(
                config, args, proteinGroupResults, parseId, peptideToProteinMaps
            )

        if args.protein_groups_out:
            writeProteinGroups(
                proteinGroupResults,
                args.protein_groups_out,
                config,
                apply_suffix=len(configs) > 1,
            )

    end = timer()
    logger.info(
        f"PickedGroupFDR execution took {'%.1f' % (end - start)} seconds wall clock time"
    )

    plotter.decoratePlots()
    plotter.savePlots()
    plotter.show()


def getProteinGroupResults(
    peptideInfoList: PeptideInfoList,
    mqProteinGroupsFile: str,
    proteinAnnotations,
    peptideToProteotypicityMap,
    pickedStrategy: ProteinCompetitionStrategy,
    scoreType: ProteinScoringStrategy,
    groupingStrategy: ProteinGroupingStrategy,
    plotter: Union[Plotter, NoPlotter],
    keep_all_proteins: bool,
):
    proteinGroups = groupingStrategy.group_proteins(
        peptideInfoList, mqProteinGroupsFile
    )

    # for razor peptide strategy
    scoreType.set_peptide_counts_per_protein(peptideInfoList)

    for rescue_step in groupingStrategy.get_rescue_steps():
        if rescue_step:
            if not scoreType.can_do_protein_group_rescue():
                raise Exception(
                    "Cannot do rescue step for other score types than bestPEP"
                )
            proteinGroups = groupingStrategy.rescue_protein_groups(
                peptideInfoList,
                proteinGroupResults,
                proteinGroups,
                proteinGroupPeptideInfos,
            )

        proteinGroupPeptideInfos = scoreType.collect_peptide_scores_per_protein(
            proteinGroups, peptideInfoList, suppressMissingProteinWarning=rescue_step
        )

        # find optimal division factor for multPEP score
        scoreType.optimize_hyperparameters(proteinGroups, proteinGroupPeptideInfos)

        if rescue_step and pickedStrategy.short_description() == "pgT":
            (
                proteinGroups,
                proteinGroupPeptideInfos,
            ) = groupingStrategy.update_protein_groups(
                proteinGroups, proteinGroupPeptideInfos
            )

        (
            pickedProteinGroups,
            pickedProteinGroupPeptideInfos,
            proteinScores,
        ) = pickedStrategy.do_competition(
            proteinGroups, proteinGroupPeptideInfos, scoreType
        )

        reportedQvals, observedQvals = fdr.calculateProteinFDRs(
            pickedProteinGroups, proteinScores
        )

        scoreCutoff = groupingStrategy.score_cutoff if rescue_step else float("inf")
        proteinGroupResults = ProteinGroupResults.from_protein_groups(
            pickedProteinGroups,
            pickedProteinGroupPeptideInfos,
            proteinScores,
            reportedQvals,
            scoreCutoff,
            proteinAnnotations,
            keep_all_proteins,
        )

    if scoreType.use_proteotypicity:
        proteotypicity.calculateProteotypicityScores(
            pickedProteinGroups,
            pickedProteinGroupPeptideInfos,
            peptideToProteotypicityMap,
            scoreType,
            scoreCutoff,
        )

    plotter.set_series_label(
        scoreType, groupingStrategy, pickedStrategy, rescue_step=rescue_step
    )
    plotter.plotQvalCalibrationAndPerformance(
        reportedQvals, observedQvals, absentRatio=1.0
    )

    return proteinGroupResults


def parseEvidenceFiles(
    evidenceFiles: List[str],
    peptideToProteinMaps: List[Dict[str, List[str]]],
    scoreType,
    suppressMissingPeptideWarning: bool,
) -> PeptideInfoList:
    """Returns best score per peptide"""
    if not scoreType.remaps_peptides_to_proteins():
        peptideToProteinMaps = [None]

    if len(peptideToProteinMaps) == 1:
        peptideToProteinMaps = peptideToProteinMaps * len(evidenceFiles)

    peptideInfoList = dict()
    for peptide, proteins, _, score in parsers.parseEvidenceFiles(
        evidenceFiles,
        peptideToProteinMaps=peptideToProteinMaps,
        scoreType=scoreType,
        suppressMissingPeptideWarning=suppressMissingPeptideWarning,
    ):
        peptide = helpers.cleanPeptide(peptide)
        if np.isnan(score) or score >= peptideInfoList.get(peptide, [np.inf])[0]:
            continue

        peptideInfoList[peptide] = [score, proteins]

    return peptideInfoList


def doQuantification(config, args, proteinGroupResults, parseId, peptideToProteinMaps):
    if not config["scoreType"].can_do_quantification():
        logger.warning(
            "Skipping quantification... Cannot do quantification from percolator output file; MQ evidence file input needed"
        )
        return

    proteinSequences = {}
    if args.fasta:
        digestion_params_list = get_digestion_params_list(args)

        logger.info("In silico protein digest for iBAQ")
        numIbaqPeptidesPerProtein = digest.getNumIbaqPeptidesPerProtein(
            args.fasta, digestion_params_list
        )
        proteinSequences = digest.getProteinSequences(args.fasta, parseId)
    elif args.peptide_protein_map:
        logger.warning("Found peptide_protein_map (instead of fasta input): ")
        logger.warning(
            "- calculating iBAQ values using all peptides in peptide_protein_map."
        )
        logger.warning("- cannot compute sequence coverage.")
        numIbaqPeptidesPerProtein = digest.getNumPeptidesPerProtein(
            digest.merge_peptide_to_protein_maps(peptideToProteinMaps)
        )
    else:
        sys.exit(
            "No fasta or peptide to protein mapping file detected, please specify either the --fasta or --peptide_protein_map flags"
        )

    quantification.doQuantification(
        args.mq_evidence,
        proteinGroupResults,
        proteinSequences,
        peptideToProteinMaps,
        numIbaqPeptidesPerProtein,
        args.file_list_file,
        config["scoreType"],
        minPeptideRatiosLFQ=args.lfq_min_peptide_ratios,
        numThreads=args.num_threads,
    )


def writeProteinGroups(proteinGroupResults: ProteinGroupResults, protein_groups_out, config, apply_suffix):
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
    proteinGroupResults.write(protein_groups_out)
    logger.info(f"Protein group results have been written to: {protein_groups_out}")


if __name__ == "__main__":
    main(sys.argv[1:])
