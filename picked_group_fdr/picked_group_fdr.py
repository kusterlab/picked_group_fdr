import sys
import os
import collections
import logging
import argparse
from timeit import default_timer as timer
from typing import List, Dict, Tuple, Union

import numpy as np

from . import version
from . import digest
from . import helpers
from . import parsers
from . import results
from . import proteotypicity
from . import quantification
from . import entrapment
from . import methods
from . import fdr

from .grouping import PseudoGeneGrouping
from .results import ProteinGroupResult, ProteinGroupResults
from .plotter import PlotterFactory
from .peptide_info import PeptideInfoList, ProteinGroupPeptideInfos

# for type hints only
from .scoring import ProteinScoringStrategy
from .grouping import ProteinGroupingStrategy
from .competition import ProteinCompetitionStrategy
from .plotter import Plotter, NoPlotter

logger = logging.getLogger(__name__)

__version__ = version.get_version_from_pyproject()
__copyright__ = '''Copyright (c) 2020-2022 Matthew The. All rights reserved.
Written by Matthew The (matthew.the@tum.de) at the
Chair of Proteomics and Bioanalytics at the Technical University of Munich.'''


class ArgumentParserWithLogger(argparse.ArgumentParser):
    def error(self, message):
        logger.error(f"Error parsing input arguments: {message}")
        super().error(message)


def parseArgs(argv):
    apars = ArgumentParserWithLogger(
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    apars.add_argument('--mq_evidence', default=None, metavar = "EV",
                         help='''MaxQuant evidence file.''')

    apars.add_argument('--protein_groups_out', default=None, metavar = "PG",
                         help='''Protein groups output file, mimicks a subset of the MQ protein groups columns.
                                ''')

    apars.add_argument('--fasta', default=None, metavar = "F",
                         help='''Fasta file to create mapping from peptides to proteins.
                                                    ''')

    apars.add_argument('--mq_protein_groups', default=None, metavar = "PG",
                         help='''MaxQuant protein groups file; only specify if we want to keep MaxQuant's original grouping instead of Picked Grouping
                                ''')

    apars.add_argument('--perc_evidence', default=None, metavar = "POUT",
                         help='''Percolator output file with PSMs or peptides; alternative for --mq_evidence if we want to use Percolator PEPs instead of MaxQuant's PEPs
                                ''')

    apars.add_argument('--methods', default=None, metavar = "M1,M2",
                         help='''Use one or more predefined protein group FDR estimation methods, separated by commas.''')

    apars.add_argument('--peptide_protein_map', default=None, metavar = "M",
                         help='''File with mapping from peptides to proteins; alternative for --fasta flag if digestion is time consuming.
                                ''')

    apars.add_argument('--peptide_proteotypicity_map', default=None, metavar = "M",
                         help='''File with mapping from peptides to proteotypicity.
                                ''')

    apars.add_argument('--gene_level',
                         help='Report gene-level statistics instead of protein group-level. This requires the GN= field to be present in the fasta file.',
                         action='store_true')

    apars.add_argument('--do_quant',
                         help='Do protein quantification, will calculate summed intensity, iBAQ and LFQ intensities.',
                         action='store_true')

    apars.add_argument('--suppress_missing_peptide_warning',
                         help='Suppress missing peptide warning.',
                         action='store_true')

    apars.add_argument('--lfq_min_peptide_ratios', default=2, type=int, metavar='M',
                         help='''Minimum number of common peptides between two samples
                                 to qualify for calculating a peptide ratio in LFQ
                                 ''')

    apars.add_argument('--num_threads', default=1, type=int, metavar='T',
                         help='''Maximum number of threads to use.''')

    apars.add_argument('--file_list_file', metavar='L',
                         help='Tab separated file with lines of the format (third and fourth columns are optional): raw_file <tab> condition <tab> experiment <tab> fraction.')

    apars.add_argument('--figure_base_fn', default=None, metavar = "F",
                         help='''Base file name for calibration and performance figures.
                                ''')

    apars.add_argument('--plot_figures',
                         help='Plot figures with matplotlib.',
                         action='store_true')

    digest.addArguments(apars)

    # ------------------------------------------------
    args = apars.parse_args(argv)

    return args


def main(argv):
    logger.info(f'PickedGroupFDR version {__version__}\n{__copyright__}')
    logger.info(f'Issued command: {os.path.basename(__file__)} {" ".join(map(str, argv))}')
    
    args = parseArgs(argv)
    
    np.random.seed(1) # set seed for random shuffling of protein groups, see competition.do_competition()
    
    configs = methods.get_methods(args)

    start = timer()
    
    plotter = PlotterFactory.get_plotter(args.figure_base_fn, args.plot_figures)
    plotter.initPlots()
    
    parseId = digest.parseUntilFirstSpace
    proteinAnnotations = digest.getProteinAnnotations(args.fasta, parseId)
    if args.gene_level:
        if digest.hasGeneNames(proteinAnnotations, minRatioWithGenes=0.5):
            parseId = digest.parseGeneNameFunc
            proteinAnnotations = digest.getProteinAnnotations(args.fasta, parseId)
        else:
            logger.warning("Found fewer than 50% of proteins without gene names in the fasta file, will attempt to infer pseudo-genes based on shared peptides instead")
            for config in configs:
                config["grouping"] = PseudoGeneGrouping()
    
    peptideToProteinMap = dict()
    peptideToProteotypicityMap = dict()
    for config in configs:
        methodDescriptionLong = methods.long_description(config['scoreType'], config['grouping'], config['pickedStrategy'], True)
        label = config.get('label', "")
        logger.info(f"Protein group level estimation method: {label} ({methodDescriptionLong})")
        
        peptideFile = config['scoreType'].get_evidence_file(args)
        if not peptideFile:
            logger.warning("No evidence file provided, skipping...")
            continue
        
        if len(peptideToProteinMap) == 0 and (config['grouping'].needs_peptide_to_protein_map() or config['scoreType'].remaps_peptides_to_proteins()):
            peptideToProteinMap = digest.get_peptide_to_protein_map(args, parseId)
            entrapment.markEntrapmentProteins(peptideToProteinMap, args.mq_protein_groups)
            
        if len(peptideToProteotypicityMap) == 0 and config['scoreType'].use_proteotypicity:
            peptideToProteotypicityMap = proteotypicity.getPeptideToProteotypicityFromFile(args.peptide_proteotypicity_map)
        
        peptideInfoList = parseMqEvidenceFile(peptideFile, 
                                              peptideToProteinMap, 
                                              config['scoreType'], 
                                              args.suppress_missing_peptide_warning)
        
        plotter.set_series_label_base(config.get('label', None))
        proteinGroupResults = getProteinGroupResults(
                peptideInfoList, args.mq_protein_groups, 
                proteinAnnotations, peptideToProteotypicityMap,
                config['pickedStrategy'], config['scoreType'], config['grouping'], 
                plotter)
        
        if args.do_quant:
            doQuantification(config, args, proteinGroupResults, parseId, peptideToProteinMap)

        if args.protein_groups_out:
            writeProteinGroups(proteinGroupResults, args.protein_groups_out, config, apply_suffix=len(configs) > 1)
    
    end = timer()
    logger.info(f"PickedGroupFDR execution took {'%.1f' % (end - start)} seconds wall clock time")

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
        plotter: Union[Plotter, NoPlotter]):
    proteinGroups = groupingStrategy.group_proteins(peptideInfoList, mqProteinGroupsFile)
    
    # for razor peptide strategy
    scoreType.set_peptide_counts_per_protein(peptideInfoList)
    
    for rescue_step in groupingStrategy.get_rescue_steps():
        if rescue_step:
            if not scoreType.can_do_protein_group_rescue():
                raise Exception("Cannot do rescue step for other score types than bestPEP")
            proteinGroups = groupingStrategy.rescue_protein_groups(peptideInfoList, proteinGroupResults, proteinGroups, proteinGroupPeptideInfos)
        
        proteinGroupPeptideInfos = scoreType.collect_peptide_scores_per_protein(
            proteinGroups, peptideInfoList, suppressMissingProteinWarning=rescue_step)
        
        # find optimal division factor for multPEP score
        scoreType.optimize_hyperparameters(proteinGroups, proteinGroupPeptideInfos)
        
        if rescue_step and pickedStrategy.short_description() == "pgT":
            proteinGroups, proteinGroupPeptideInfos = groupingStrategy.update_protein_groups(proteinGroups, proteinGroupPeptideInfos)
        
        pickedProteinGroups, pickedProteinGroupPeptideInfos, proteinScores = pickedStrategy.do_competition(proteinGroups, proteinGroupPeptideInfos, scoreType)
        
        reportedQvals, observedQvals = fdr.calculateProteinFDRs(
                pickedProteinGroups, proteinScores)
        
        scoreCutoff = groupingStrategy.score_cutoff if rescue_step else float("inf")
        proteinGroupResults = ProteinGroupResults.from_protein_groups(
                pickedProteinGroups, pickedProteinGroupPeptideInfos, 
                proteinScores, reportedQvals, 
                scoreCutoff, proteinAnnotations)

    if scoreType.use_proteotypicity:
        proteotypicity.calculateProteotypicityScores(pickedProteinGroups, pickedProteinGroupPeptideInfos, peptideToProteotypicityMap, scoreType, scoreCutoff)
    
    plotter.set_series_label(scoreType, groupingStrategy, pickedStrategy, rescue_step=rescue_step)
    plotter.plotQvalCalibrationAndPerformance(reportedQvals, observedQvals, absentRatio=1.0)
        
    return proteinGroupResults


def parseMqEvidenceFile(mqEvidenceFile: str, peptideToProteinMap, scoreType, suppressMissingPeptideWarning: bool) -> PeptideInfoList:
    """Returns best score per peptide"""
    peptideInfoList = dict()
    for peptide, tmp_proteins, _, score in parsers.parseMqEvidenceFile(mqEvidenceFile, scoreType = scoreType):
        peptide = helpers.cleanPeptide(peptide)
        if np.isnan(score) or score >= peptideInfoList.get(peptide, [np.inf])[0]:
            continue
        
        if scoreType.remaps_peptides_to_proteins():
            proteins = digest.getProteins(peptideToProteinMap, peptide)
            if len(proteins) == 0:
                if not helpers.isContaminant(tmp_proteins) and not suppressMissingPeptideWarning:
                    logger.warning(f"Missing peptide: {peptide} {tmp_proteins}")
                continue
        else:
            proteins = tmp_proteins

        proteins = helpers.removeDecoyProteinsFromTargetPeptides(proteins)
        peptideInfoList[peptide] = [score, proteins]

    return peptideInfoList


def doQuantification(config, args, proteinGroupResults, parseId, peptideToProteinMap):
    if not config['scoreType'].can_do_quantification():
        logger.warning("Skipping quantification... Cannot do quantification from percolator output file; MQ evidence file input needed")
        return
    
    logger.info("In silico protein digest for iBAQ")
    numIbaqPeptidesPerProtein = digest.getNumIbaqPeptidesPerProtein(args)
    proteinSequences = digest.getProteinSequences(args.fasta, parseId)
    
    quantification.doQuantification(
            args.mq_evidence, proteinGroupResults, proteinSequences,
            peptideToProteinMap, numIbaqPeptidesPerProtein, 
            args.file_list_file, config['scoreType'], minPeptideRatiosLFQ=args.lfq_min_peptide_ratios, numThreads=args.num_threads)


def writeProteinGroups(proteinGroupResults, protein_groups_out, config, apply_suffix):
    if apply_suffix:
        base, ext = os.path.splitext(protein_groups_out)
        label = config.get('label', None)
        if label is None:
            label = methods.short_description(config['scoreType'], config['grouping'], config['pickedStrategy'], True)
        else:
            label = label.lower().replace(" ", "_")
        protein_groups_out = f"{base}_{label}{ext}" 
    proteinGroupResults.write(protein_groups_out)
    logger.info(f"Protein group results have been written to: {protein_groups_out}")


if __name__ == "__main__":
    main(sys.argv[1:])
