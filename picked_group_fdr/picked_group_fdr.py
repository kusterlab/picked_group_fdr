import sys
import os
import collections
import logging

import numpy as np

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


logger = logging.getLogger(__name__)

__version__ = "0.1.0"
__copyright__ = '''Copyright (c) 2020-2022 Matthew The. All rights reserved.
Written by Matthew The (matthew.the@tum.de) at the
Chair of Proteomics and Bioanalytics at the Technical University of Munich.'''


def parseArgs(argv):
  import argparse
  apars = argparse.ArgumentParser(
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  apars.add_argument('--mq_evidence', default=None, metavar = "EV",
                     help='''MaxQuant evidence file.
                          ''')
  
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
  
  '''
  pickedStrategy:     PickedStrategy() = picked FDR
                      ClassicStrategy() = classic FDR
  
  grouping:           MQNativeGrouping() = MaxQuant grouping from proteinGroups.txt
                      SubsetGrouping() = Emulate protein grouping of MaxQuant based on evidence.txt
                                   (currently does not work with simulated datasets since peptideToProteinMap does not contain entrapment labels)
                      NoGrouping() = No protein grouping, each protein is in its own group
                      +Rescued = Rescue protein groups by only considering peptides below 1% protein FDR threshold
  '''
  
  configs = methods.get_methods(args)
  
  from timeit import default_timer as timer

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
    config['scoreType'].set_peptide_counts_per_protein(peptideInfoList) # for razor peptides
    
    plotter.set_series_label_base(config.get('label', None))
    proteinGroupResults = getProteinGroupResults(
        peptideInfoList, args.mq_protein_groups, 
        proteinAnnotations, peptideToProteotypicityMap,
        config['pickedStrategy'], config['scoreType'], config['grouping'], 
        plotter)
    
    if args.do_quant:
      if not config['scoreType'].can_do_quantification():
        logger.warning("Skipping quantification... Cannot do quantification from percolator output file; MQ evidence file input needed")
        continue
      
      logger.info("In silico protein digest for iBAQ")
      numIbaqPeptidesPerProtein = digest.getNumIbaqPeptidesPerProtein(args)
      proteinSequences = digest.getProteinSequences(args.fasta, parseId)
      
      quantification.doQuantification(
          args.mq_evidence, proteinGroupResults, proteinSequences,
          peptideToProteinMap, numIbaqPeptidesPerProtein, 
          args.file_list_file, config['scoreType'], minPeptideRatiosLFQ = args.lfq_min_peptide_ratios)
    
    if args.protein_groups_out:
      protein_groups_out = args.protein_groups_out
      if len(configs) > 1:
        base, ext = os.path.splitext(protein_groups_out)
        label = config.get('label', None)
        if label is None:
          label = methods.short_description(config['scoreType'], config['grouping'], config['pickedStrategy'], True)
        else:
          label = label.lower().replace(" ", "_")
        protein_groups_out = f"{base}_{label}{ext}" 
      proteinGroupResults.write(protein_groups_out)
      logger.info(f"Protein group results have been written to: {protein_groups_out}")
  
  end = timer()
  logger.info(f"PickedGroupFDR execution took {'%.1f' % (end - start)} seconds wall clock time")

  plotter.decoratePlots()
  plotter.savePlots()
  plotter.show()


def getProteinGroupResults(
    peptideInfoList, mqProteinGroupsFile, proteinAnnotations, peptideToProteotypicityMap, 
    pickedStrategy, scoreType, groupingStrategy, plotter):    
  proteinGroups = groupingStrategy.group_proteins(peptideInfoList, mqProteinGroupsFile)
  
  for rescue_step in groupingStrategy.get_rescue_steps():
    if rescue_step:
      if not scoreType.can_do_protein_group_rescue():
        raise Exception("Cannot do rescue step for other score types than bestPEP")
      proteinGroups = groupingStrategy.rescue_protein_groups(peptideInfoList, proteinGroupResults, proteinGroups)
    
    proteinGroupScores = scoreType.collect_peptide_scores_per_protein(
      proteinGroups, peptideInfoList, suppressMissingProteinWarning=rescue_step)
    
    scoreType.optimize_hyperparameters(proteinGroups, proteinGroupScores)
    
    pickedProteinGroups, proteinGroupScores = pickedStrategy.do_competition(proteinGroups, proteinGroupScores, scoreType)
    
    reportedQvals, observedQvals, proteinScores = fdr.calculateProteinFDRs(
        pickedProteinGroups, proteinGroupScores, scoreType)
    
    if len(groupingStrategy.get_rescue_steps()) == 1 or rescue_step:
      plotter.set_series_label(scoreType, groupingStrategy, pickedStrategy, rescue_step=rescue_step)
      
      #absentRatio = 0.9 # entrapment / (pool + entrapment)
      #absentRatio = 96788.0 / (16503 + 96788.0)
      absentRatio = 1.0
      plotter.plotQvalCalibrationAndPerformance(reportedQvals, observedQvals, absentRatio)
    
    scoreCutoff = groupingStrategy.score_cutoff if rescue_step else float("inf")
    proteinGroupResults = ProteinGroupResults.from_protein_groups(
        pickedProteinGroups, proteinGroupScores, 
        proteinScores, reportedQvals, 
        scoreCutoff, proteinAnnotations)
    
    if scoreType.use_proteotypicity:
      proteotypicity.calculateProteotypicityScores(pickedProteinGroups, proteinGroupScores, peptideToProteotypicityMap, scoreType, scoreCutoff)
    
  return proteinGroupResults


def parseMqEvidenceFile(mqEvidenceFile, peptideToProteinMap, scoreType, suppressMissingPeptideWarning):
  peptideInfoList = dict()
  for peptide, tmp_proteins, _, score in parsers.parseMqEvidenceFile(mqEvidenceFile, scoreType = scoreType):
    peptide = helpers.cleanPeptide(peptide)
    if score < peptideInfoList.get(peptide, [np.inf])[0]:
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


if __name__ == "__main__":
  main(sys.argv[1:])
