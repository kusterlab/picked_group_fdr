import sys
import os
import logging
import collections
from typing import List

import numpy as np
import triqler.parsers

from . import digest
from . import digestion_params
from . import helpers
from . import fdr
from .parsers import maxquant
from .parsers import psm

from .protein_groups import ProteinGroups
from .scoring import ProteinScoringStrategy
from .quant.base import ProteinGroupColumns
from .quant.precursor_quant import PrecursorQuant
from .quant.lfq import LFQIntensityColumns
from .quant.sum_and_ibaq import SummedIntensityAndIbaqColumns
from .quant.peptide_count import UniquePeptideCountColumns
from .quant.id_type import IdentificationTypeColumns
from .quant.tmt import TMTIntensityColumns
from .quant.triqler import TriqlerIntensityColumns
from .quant.sequence_coverage import SequenceCoverageColumns
from .quant.evidence_ids import EvidenceIdsColumns


logger = logging.getLogger(__name__)


def parseArgs(argv):
    import argparse
    apars = argparse.ArgumentParser(
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    apars.add_argument('--mq_evidence', default=None, metavar = "EV", required = True, nargs="+",
                         help='''MaxQuant evidence file.''')
    
    apars.add_argument('--mq_protein_groups', default=None, metavar = "PG", required = True,
                         help='''MaxQuant protein groups file.''')
    
    apars.add_argument('--protein_groups_out', default=None, metavar = "PG", required = True,
                         help='''Protein groups output file, mimicks a subset of the MQ protein groups columns.
                                ''')
    
    apars.add_argument('--peptide_protein_map', default=None, metavar = "M",
                         help='''TSV file with mapping from peptides to proteins.''')
    
    apars.add_argument('--fasta', default=None, metavar = "F",
                         help='''Fasta file to create mapping from peptides to proteins. This should not contain the decoy sequences, unless you set the --fasta_contains_decoys flag.''')
    
    apars.add_argument('--fasta_use_uniprot_id',
                         help='''Parse protein identifiers in the fasta file as UniProt IDs, 
                                 i.e. Q9UM47 for the protein identifier sp|Q9UM47|NOTC3_HUMAN''',
                         action='store_true')

    apars.add_argument('--gene_level',
                         help='''Do quantification on gene-level instead of on protein group level''',
                         action='store_true')
    
    apars.add_argument('--file_list_file', metavar='L', 
                         help='''Tab separated file with lines of the format (third and 
                                 fourth columns are optional): raw_file <tab> condition 
                                 <tab> experiment <tab> fraction.
                                 ''',
                         required = False)
    
    apars.add_argument('--lfq_min_peptide_ratios', default=2, type=int, metavar='M',
                         help='''Minimum number of common peptides between two samples
                                 to qualify for calculating a peptide ratio in LFQ
                                 ''')
    
    apars.add_argument('--num_threads', default=1, type=int, metavar='T',
                         help='''Maximum number of threads to use.''')
    
    apars.add_argument('--lfq_stabilize_large_ratios',
                         help='''Apply stabilization of large ratios in LFQ as described
                                 in the MaxLFQ paper.''',
                         action='store_false')
    
    
    digestion_params.add_digestion_arguments(apars)
                                
    # ------------------------------------------------
    args = apars.parse_args(argv)
    
    return args


def main(argv):
    logger.info(f'Issued command: {os.path.basename(__file__)} {" ".join(map(str, argv))}')
    
    args = parseArgs(argv)
    
    parseId = digest.parseUntilFirstSpace
    if args.gene_level:
        parseId = digest.parseGeneNameFunc
    elif args.fasta_use_uniprot_id:
        parseId = digest.parseUniProtId
    
    peptideToProteinMap, numIbaqPeptidesPerProtein = getPeptideToProteinMaps(args, parseId)
    proteinSequences = digest.getProteinSequences(args.fasta, parseId)
    proteinGroupResults = maxquant.parse_mq_protein_groups_file(args.mq_protein_groups)
    
    scoreType = ProteinScoringStrategy("bestPEP")
    doQuantification(args.mq_evidence, proteinGroupResults, proteinSequences,
                     peptideToProteinMap, numIbaqPeptidesPerProtein, args.file_list_file, 
                     scoreType, minPeptideRatiosLFQ = args.lfq_min_peptide_ratios,
                     stabilizeLargeRatiosLFQ = args.lfq_stabilize_large_ratios,
                     numThreads = args.num_threads)
    
    proteinGroupResults.write(args.protein_groups_out)
    
    logger.info(f"Protein group results have been written to: {args.protein_groups_out}")


def getPeptideToProteinMaps(args, parseId):
    minLenIbaq = max([6, args.min_length])
    maxLenIbaq = min([30, args.max_length])
    
    logger.info("Loading peptide to protein map...")
    if args.fasta:
        pre, not_post, post = digest.getCleavageSites(args.enzyme)
        
        db = 'concat'
        if args.fasta_contains_decoys:
            db = 'target'
            
        peptideToProteinMap = digest.getPeptideToProteinMap(
                args.fasta, db = db, digestion = args.digestion, 
                min_len = args.min_length, max_len = args.max_length, 
                pre = pre, not_post = not_post, post = post, miscleavages = args.cleavages, 
                methionineCleavage = True, specialAAs = list(args.special_aas),
                parseId = parseId, useHashKey = (args.digestion == "none"))
        
        peptideToProteinMapIbaq = digest.getPeptideToProteinMap(
                args.fasta, db = db, digestion = args.digestion, 
                min_len = minLenIbaq, max_len = maxLenIbaq, 
                pre = pre, not_post = not_post, post = post, miscleavages = 0,
                methionineCleavage = False, specialAAs = list(args.special_aas),
                parseId = parseId, useHashKey = (args.digestion == "none"))
    elif args.peptide_protein_map:
        pre, not_post, post = digest.getCleavageSites(args.enzyme)
        
        peptideToProteinMap = digest.getPeptideToProteinMapFromFile(
                args.peptide_protein_map, useHashKey = True)
        
        peptideToProteinMapIbaq = dict()
        for peptide, proteins in peptideToProteinMap.items():
            peptideLen = len(peptide)
            if (peptideLen >= minLenIbaq and 
                    peptideLen <= maxLenIbaq and 
                    not digest.hasMiscleavage(peptide, pre, not_post, post)):
                peptideToProteinMapIbaq[peptide] = proteins
    else:
        sys.exit('No peptide to protein map found, use either the --fasta or the --peptide_protein_map arguments')
    
    numIbaqPeptidesPerProtein = digest.getNumPeptidesPerProtein(peptideToProteinMapIbaq)
    
    return peptideToProteinMap, numIbaqPeptidesPerProtein


def doQuantification(mqEvidenceFiles, proteinGroupResults, proteinSequences,
        peptideToProteinMaps, numIbaqPeptidesPerProtein, fileListFile, 
        scoreType, psmQvalCutoff = 0.01, 
        discardSharedPeptides = True, 
        minPeptideRatiosLFQ = 2, 
        stabilizeLargeRatiosLFQ = True,
        numThreads = 1):
    params = initTriqlerParams()
    
    logger.info("Preparing for quantification")
    fileMapping = None
    if fileListFile:
        experiments, fileMapping, params = parseFileList(fileListFile, params) 
    
    proteinGroupResults, postErrProbs, numTmtChannels, numSilacChannels, parsedExperiments = parseEvidenceFiles(
            proteinGroupResults, mqEvidenceFiles, peptideToProteinMaps, 
            fileMapping, scoreType, discardSharedPeptides)
    
    silacChannels = getSilacChannels(numSilacChannels)
    
    if len(parsedExperiments) > 0:
        experiments = sorted(list(parsedExperiments))
    experimentToIdxMap = dict([(v,k) for k, v in enumerate(experiments)])
    
    # (1) technically this is a precursor-level FDR and not a PSM-level FDR
    # (2) in contrast to MaxQuant, we set a global precursor-level FDR 
    #         instead of a per raw file PSM-level FDR
    postErrProbCutoff = fdr.calcPostErrProbCutoff([x[0] for x in postErrProbs if not helpers.isMbr(x[0])], psmQvalCutoff)
    logger.info(f"PEP-cutoff corresponding to {psmQvalCutoff*100:g}% PSM-level FDR: {postErrProbCutoff}")
    
    printNumPeptidesAtFDR(postErrProbs, postErrProbCutoff)
    
    logger.info("Filtering for identified precursors")
    # precursor = (peptide, charge) tuple
    # this filter also ensures that MBR precursors which were matched to 
    # unidentified precursors are removed
    for pgr in proteinGroupResults:
        pgr.precursorQuants = retainOnlyIdentifiedPrecursors(pgr.precursorQuants, postErrProbCutoff)
    
    columns: List[ProteinGroupColumns] = [
        UniquePeptideCountColumns(), 
        IdentificationTypeColumns(),
        SummedIntensityAndIbaqColumns(silacChannels, numIbaqPeptidesPerProtein),
        SequenceCoverageColumns(proteinSequences),
        EvidenceIdsColumns()
    ]
    
    if numTmtChannels > 0:
        columns.append(TMTIntensityColumns(numTmtChannels))
    else:
        columns.append(LFQIntensityColumns(silacChannels, minPeptideRatiosLFQ, stabilizeLargeRatiosLFQ, numThreads))
        # TODO: add SILAC functionality of Triqler
        if numSilacChannels == 0:
            columns.append(TriqlerIntensityColumns(params))
    
    for c in columns:
        c.append_headers(proteinGroupResults, experiments)
        c.append_columns(proteinGroupResults, experimentToIdxMap, postErrProbCutoff)


def parseEvidenceFiles(proteinGroupResults, mqEvidenceFiles, peptideToProteinMaps, 
                                            fileMapping, scoreType, discard_shared_peptides):    
    protein_groups = ProteinGroups.from_protein_group_results(proteinGroupResults)
    protein_groups.create_index()
    
    postErrProbs = list()
    shared_peptide_precursors, unique_peptide_precursors = 0, 0
    numTmtChannels, numSilacChannels = -1, -1
    parsedExperiments = set()
    missing_peptides_in_protein_groups = 0
    
    for peptide, proteins, charge, rawFile, experiment, fraction, intensity, postErrProb, tmtCols, silacCols, evidenceId in psm.parse_evidence_file_multiple(mqEvidenceFiles, peptide_to_protein_maps = peptideToProteinMaps, score_type = ProteinScoringStrategy("bestPEP"), for_quantification = True):
        if numTmtChannels == -1:
            # There are 3 columns per TMT channel: 
            #     Reporter intensity corrected, 
            #     Reporter intensity
            #     Reporter intensity count
            numTmtChannels = int(len(tmtCols) / 3) 
        if numSilacChannels == -1:
            numSilacChannels = len(silacCols)

        # override the parsed experiment and fraction if --file_list_file option is used
        if fileMapping:
            experiment, fraction = fileMapping[rawFile]
        elif experiment not in parsedExperiments:
            parsedExperiments.add(experiment)        
        
        protein_group_idxs = protein_groups.get_protein_group_idxs(proteins)
        
        # removes peptides not present in the proteinGroups.txt file
        if len(protein_group_idxs) == 0:
            logger.debug(f'Could not find any of the proteins {proteins} in proteinGroups.txt')
            missing_peptides_in_protein_groups += 1
            continue
        
        if discard_shared_peptides and helpers.is_shared_peptide(protein_group_idxs):
            shared_peptide_precursors += 1
            continue
        
        unique_peptide_precursors += 1
    
        if not helpers.isDecoy(proteins):
            postErrProbs.append((postErrProb, rawFile, experiment, peptide))
        
        if len(tmtCols) > 0:
            tmtCols = np.array(tmtCols, dtype = 'float64')
        else:
            tmtCols = None
        
        if len(silacCols) > 0:
            silacCols = np.array(silacCols, dtype = 'float64')
        else:
            silacCols = None

        for proteinGroupIdx in protein_group_idxs:
            precursorQuant = PrecursorQuant(peptide, 
                                            charge, 
                                            experiment, 
                                            fraction, 
                                            intensity, 
                                            postErrProb,
                                            tmtCols,
                                            silacCols,
                                            evidenceId)
            proteinGroupResults[proteinGroupIdx].precursorQuants.append(precursorQuant)
    
    # if missingPeptidesInFasta > 0:
    #     logger.warning(f"Skipped {missingPeptidesInFasta} precursors not present in the fasta file")
    
    if missing_peptides_in_protein_groups > 0:
        logger.debug(f"Skipped {missing_peptides_in_protein_groups} precursors from proteins not present in proteinGroups.txt file")
    
    logger.info(f"Found {unique_peptide_precursors} precursors from unique and {shared_peptide_precursors} precursors from shared peptides")
    
    return proteinGroupResults, postErrProbs, numTmtChannels, numSilacChannels, parsedExperiments

    
def printNumPeptidesAtFDR(postErrProbs, postErrProbCutoff):
    survivingModPeptides = set([x[3] for x in postErrProbs if x[0] <= postErrProbCutoff])
    
    peptidesPerRawFile = collections.defaultdict(list)
    peptidesPerExperiment = collections.defaultdict(list)
    peptidesPerRawFileMbr = collections.defaultdict(list)
    peptidesPerExperimentMbr = collections.defaultdict(list)
    for postErrProb, rawFile, experiment, peptide in postErrProbs:
        if postErrProb <= postErrProbCutoff:
            peptidesPerRawFile[rawFile].append(peptide)
            peptidesPerExperiment[experiment].append(peptide)
        elif helpers.isMbr(postErrProb) and peptide in survivingModPeptides:
            peptidesPerRawFileMbr[rawFile].append(peptide)
            peptidesPerExperimentMbr[experiment].append(peptide)

    logger.info("Precursor counts per rawfile (1% PSM-level FDR):")
    for rawFile, peptides in sorted(peptidesPerRawFile.items()):
        numPeptides = len(set(peptides))
        numPeptidesWithMbr = len(set(peptides + peptidesPerRawFileMbr[rawFile]))
        logger.info(f"    {rawFile}: {numPeptides} {'(' + str(numPeptidesWithMbr) + ' with MBR)' if numPeptidesWithMbr > numPeptides else ''}")
    
    logger.info("Precursor counts per experiment (1% PSM-level FDR):")
    for experiment, peptides in sorted(peptidesPerExperiment.items()):
        numPeptides = len(set(peptides))
        numPeptidesWithMbr = len(set(peptides + peptidesPerExperimentMbr[experiment]))
        logger.info(f"    {experiment}: {numPeptides} {'(' + str(numPeptidesWithMbr) + ' with MBR)' if numPeptidesWithMbr > numPeptides else ''}")


def getSilacChannels(numSilacChannels):
    silacChannels = list()
    if numSilacChannels == 3:
        silacChannels = ['L', 'M', 'H']
    elif numSilacChannels == 2:
        silacChannels = ['L', 'H']
    elif numSilacChannels != 0:
        sys.exit("ERROR: Found a number of SILAC channels not equal to 2 or 3")
    return silacChannels


def initTriqlerParams():
    params = dict()
    # TODO: make these parameters configurable from the command line
    params['decoyPattern'] = 'REV__'
    params["groups"] = []
    params["groupLabels"] = []
    params['numThreads'] = 4
    params['warningFilter'] = "ignore"
    params['foldChangeEval'] = 0.8
    params['returnPosteriors'] = False
    params["minSamples"] = 5
    return params


def parseFileList(fileListFile, params):
    fileInfoList = triqler.parsers.parseFileList(fileListFile)
    fileMapping = dict()
    experiments = list()
    for rawFile, condition, experiment, fraction in fileInfoList:
        if experiment not in experiments:
            experiments.append(experiment)
        fileMapping[rawFile] = (experiment, fraction)
        # Note that params["groupLabels"] and params["groups"] are only used by Triqler
        if condition not in params["groupLabels"]:
            params["groupLabels"].append(condition)
            params["groups"].append([])
        params["groups"][params["groupLabels"].index(condition)].append(experiments.index(experiment))
    return experiments, fileMapping, params


def retainOnlyIdentifiedPrecursors(peptideIntensityList, postErrProbCutoff):
    identifiedPrecursors = set()
    for precursor in peptideIntensityList:
        if precursor.postErrProb <= postErrProbCutoff:
            identifiedPrecursors.add((precursor.peptide, precursor.charge))
    return [precursorRow for precursorRow in peptideIntensityList if (precursorRow.peptide, precursorRow.charge) in identifiedPrecursors]


if __name__ == "__main__":
    main(sys.argv[1:])
