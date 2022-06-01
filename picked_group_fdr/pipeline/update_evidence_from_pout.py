import sys
import os
import collections
import logging
from typing import List, Dict, Set

from .. import parsers
from .. import helpers


# hacky way to get the package logger instead of just __main__ when running as python -m picked_group_fdr.pipeline.update_evidence_from_pout ...
logger = logging.getLogger(__package__ + "." + __file__)


def parseArgs(argv):
    import argparse
    apars = argparse.ArgumentParser(
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    apars.add_argument('--mq_evidence', default=None, metavar = "EV", nargs='+', required = True,
                                         help='''MaxQuant evidence.txt or msms.txt file(s). If you want to combine multiple evidence files, use spaces to separate the file paths.
                                                    ''')
                                                    
    apars.add_argument('--mq_evidence_out', default=None, metavar = "EV", required = True,
                                         help='''MaxQuant evidence.txt or msms.txt combined output file with updated percolator results.
                                                    ''')
    
    apars.add_argument('--perc_results', default=list(), metavar = "POUT", nargs='+',
                                         help='''Percolator output file(s) with PSMs separated by spaces. Also include the decoy results. If this flag is not set, the evidence files are simply concatenated without updating the PSMs.
                                                    ''')
                                                    
    apars.add_argument('--mq_msms', default=None, metavar = "M", nargs='+',
                                         help='''MaxQuant msms.txt file(s) (optional). Can help resolve MS1 features with multiple scans. If you want to combine multiple evidence files, use spaces to separate the file paths. Not implemented yet...
                                                    ''')
    
    apars.add_argument('--pout_input_type', default="andromeda", metavar = "SE",
                                         help='''Input type of percolator output file. Can be "andromeda" or "prosit".
                                                    ''')
    
    apars.add_argument('--mq_input_type', default="auto", metavar = "IT",
                                         help='''Input type of MaxQuant output file. Can be "peptides", "evidence", "msms" or "auto" (automatic detection).
                                                    ''')
                                                    
    # ------------------------------------------------
    args = apars.parse_args(argv)
    
    return args


def main(argv):
    args = parseArgs(argv)
    
    if os.path.isfile(args.mq_evidence_out):
        logger.info(f"Found updated evidence file {args.mq_evidence_out}, remove this file to rerun update_evidence.")
        return
    
    if args.mq_input_type == "peptides":
        updatePeptides(args.mq_evidence, args.perc_results, args.mq_evidence_out, args.mq_msms, args.pout_input_type)
    else:
        updateEvidence(args.mq_evidence, args.perc_results, args.mq_evidence_out, args.mq_msms, args.pout_input_type)

    
def updateEvidence(evidenceFiles, poutFiles, outEvidenceFile, msmsFiles, poutInputType):
    fixed_mods, resultsDict = get_percolator_results(poutFiles, poutInputType)
    
    logger.info("Writing updated combined evidence file")
    writer = parsers.getTsvWriter(outEvidenceFile + ".tmp")
    first = True
    mqPEPs, prositPEPs = list(), list()
    unexplainedMissingPSMs = 0
    unexplainedPeptides = list()
    for evidenceFile in evidenceFiles:
        logger.info(f"Processing {evidenceFile}")
        reader = parsers.getTsvReader(evidenceFile)
        headersOrig = next(reader) # save the header
        headers = list(map(lambda x : x.lower(), headersOrig))
        
        if first:
            first = False
            writer.writerow(headersOrig)
            firstHeaders = headers
        else:
            if headers != firstHeaders:
                logger.warning("Current column names are different from the first evidence file")
                logger.info('Column\tFirstHeaders\tCurrentHeaders')
                logger.info('\n'.join([str(i+1) + '\t' + x + '\t' + y for i, (x, y) in enumerate(zip(firstHeaders, headers)) if x != y]))
        
        # these columns will be updated
        scoreCol = headers.index('score')
        postErrProbCol = headers.index('pep')
        
        # these columns are needed to retrieve the PSM
        rawFileCol = headers.index('raw file')
        if 'ms/ms scan number' in headers:
            scanNrCol = headers.index('ms/ms scan number') # evidence.txt
        else:
            scanNrCol = headers.index('scan number') # msms.txt
        peptCol = headers.index('modified sequence')
        idTypeCol = headers.index('type') # MULTI-MSMS MULTI-MATCH MSMS MULTI-SECPEP MULTI-MATCH-MSMS
        reverseCol = headers.index('reverse')
        labelingStateCol = None
        if 'labeling state' in headers:
            labelingStateCol = headers.index('labeling state')

        missingRawFiles = set()
        for row in reader:
            if len(poutFiles) > 0 and len(row[scanNrCol]) > 0:
                rawFile = row[rawFileCol]
                scanNr = int(row[scanNrCol])
                peptideOriginal = row[peptCol]
                isDecoy = (row[reverseCol] == "+")
                idType = row[idTypeCol]

                if poutInputType == "prosit":
                    fixed_mods_tmp = fixed_mods
                    if is_heavy_labeled(row, labelingStateCol):
                        fixed_mods_tmp = parsers.SILAC_HEAVY_FIXED_MODS
                    peptide = parsers.maxquant_to_internal([peptideOriginal], fixed_mods=fixed_mods_tmp)[0]
                else:
                    peptide = peptideOriginal[1:-1]
                
                if len(resultsDict[rawFile]) == 0 and rawFile not in missingRawFiles:
                    logger.warning(f"Could not find any PSMs for raw file {rawFile} in the percolator result files")
                    missingRawFiles.add(rawFile)
                    continue

                percResult = resultsDict[rawFile].get((scanNr, peptide), None)
                if not percResult and poutInputType == "prosit" and has_unknown_silac_label(row, labelingStateCol): # try with heavy labeling for SILAC because labelingState column is not reliable when multiple MS/MS map to the precursor
                    fixed_mods_tmp = parsers.SILAC_HEAVY_FIXED_MODS
                    peptide = parsers.maxquant_to_internal([peptideOriginal], fixed_mods=fixed_mods_tmp)[0]
                    percResult = resultsDict[rawFile].get((scanNr, peptide), None)
                    
                if not isDecoy and rawFile in resultsDict:
                    mqPEPs.append((float(row[postErrProbCol]), idType))
                    if not percResult and (poutInputType != "prosit" or is_valid_prosit_peptide(peptide)): # mysterious missing PSMs
                        unexplainedMissingPSMs += 1
                        #logger.info(rawFile, peptide, scanNr)
                        if unexplainedMissingPSMs <= 10:
                            unexplainedPeptides.append(row[peptCol])
                
                if percResult:
                    if not isDecoy:
                        prositPEPs.append((percResult[1], idType))
                    row[scoreCol] = percResult[0]
                    row[postErrProbCol] = percResult[1]
                else:
                    #logger.info("Missing PSM in percolator output: " + str(row[rawFileCol]) + " " + str(row[scanNrCol]) + " " + str(row[peptCol]))
                    continue
            
            writer.writerow(row)
    
    logger.info(f"Unexplained missing PSMs in Percolator results: {unexplainedMissingPSMs}")
    if unexplainedMissingPSMs > 0:
        logger.info("\tFirst 10 missing peptides:")
        for peptide in unexplainedPeptides:
            logger.info("\t" + peptide)
    logger.info(f"Results written to {outEvidenceFile}")
    
    logger.info("#MQ identifications:")
    countBelowFDR(mqPEPs)
    
    logger.info(f"#{poutInputType} identifications:")
    countBelowFDR(prositPEPs)
    
    os.rename(outEvidenceFile + ".tmp", outEvidenceFile)


def updatePeptides(peptideFiles, poutFiles, outPeptideFile, msmsFiles, poutInputType):
    _, resultsDict = get_percolator_results(poutFiles, poutInputType)
        
    peptideResultsDict = convertPSMDictToPeptideDict(resultsDict)
    
    logger.info("Writing updated peptides file")
    writer = parsers.getTsvWriter(outPeptideFile + ".tmp")
    first = True
    mqPEPs, prositPEPs = list(), list()
    unexplainedMissingPeptides = 0
    unexplainedPeptides = list()
    for peptideFile in peptideFiles:
        logger.info(f"Processing {peptideFile}")
        reader = parsers.getTsvReader(peptideFile)
        headersOrig = next(reader) # save the header
        headers = list(map(lambda x : x.lower(), headersOrig))
        
        if first:
            first = False
            writer.writerow(headersOrig)
            firstHeaders = headers
        else:
            if headers != firstHeaders:
                logger.warning("Current column names are different from the first evidence file")
                logger.info('Column\tFirstHeaders\tCurrentHeaders')
                logger.info('\n'.join([str(i+1) + '\t' + x + '\t' + y for i, (x, y) in enumerate(zip(firstHeaders, headers)) if x != y]))
        
        scoreCol = headers.index('score')
        postErrProbCol = headers.index('pep')
        
        peptCol = headers.index('sequence')
        reverseCol = headers.index('reverse')
        
        for row in reader:
            if len(poutFiles) > 0:
                peptide = row[peptCol]
                isDecoy = (row[reverseCol] == "+")
                
                percResult = peptideResultsDict.get(peptide, None)
                    
                if not isDecoy:
                    mqPEPs.append((float(row[postErrProbCol]), "Unknown"))
                    if not percResult and (poutInputType != "prosit" or is_valid_prosit_peptide(peptide)): # mysterious missing PSMs
                        unexplainedMissingPeptides += 1
                        #logger.info(rawFile, peptide, scanNr)
                        if unexplainedMissingPeptides <= 10:
                            unexplainedPeptides.append(row[peptCol])
                
                if percResult:
                    if not isDecoy:
                        prositPEPs.append((percResult[1], "Unknown"))
                    row[scoreCol] = percResult[0]
                    row[postErrProbCol] = percResult[1]
                else:
                    #logger.info("Missing PSM in percolator output: " + str(row[rawFileCol]) + " " + str(row[scanNrCol]) + " " + str(row[peptCol]))
                    continue
            
            writer.writerow(row)
    
    logger.info("Unexplained missing peptides in Percolator results: {unexplainedMissingPeptides}")
    if unexplainedMissingPeptides > 0:
        logger.info("\tFirst 10 missing peptides:")
        for peptide in unexplainedPeptides:
            logger.info("\t" + peptide)
    logger.info(f"Results written to {outPeptideFile}")
    
    logger.info("MQ")
    countBelowFDR(mqPEPs)
    
    logger.info(poutInputType)
    countBelowFDR(prositPEPs)
    
    os.rename(outPeptideFile + ".tmp", outPeptideFile)


def get_percolator_results(pout_files: List[str], pout_input_type: str):
    results_dict = collections.defaultdict(dict)
    fixed_mods = None
    for pout_file in pout_files:
        logger.info(f"Processing {pout_file}")
        fixed_mods, results_dict = parsers.parsePercolatorOutFileToDict(pout_file, results_dict, pout_input_type)
    
    logger.info("Finished parsing percolator output files")
    logger.info("#PSMs per raw file:")
    for raw_file, psms in results_dict.items():
        logger.info(f"  {raw_file}: {len(psms)} PSMs")
    return fixed_mods, results_dict


def convertPSMDictToPeptideDict(resultsDict):
    peptideResultsDict = dict()
    for _, results in resultsDict.items():
        for (_, peptide), (score, PEP) in results.items():
            peptide = helpers.cleanPeptide(peptide, removeFlanks = False)
            currScore, currPEP = peptideResultsDict.get(peptide, (-1e10, 1e10))
            peptideResultsDict[peptide] = (max(currScore, score), min(currPEP, PEP))
    return peptideResultsDict


def is_valid_prosit_peptide(peptide):
    return (len(peptide) <= 30 and not "(ac)" in peptide)


def has_unknown_silac_label(row, labeling_state_col):
    return labeling_state_col is not None and len(row[labeling_state_col]) == 0


def is_heavy_labeled(row, labeling_state_col):
    return labeling_state_col is not None and len(row[labeling_state_col]) > 0 and int(row[labeling_state_col]) == 1


def countBelowFDR(pepsWithType, fdr = 0.01):
    pepsWithType = sorted(pepsWithType)
    counts = collections.defaultdict(int)
    sumPEP = 0.0
    for i, (pep, idType) in enumerate(pepsWithType):
        sumPEP += pep
        if sumPEP / (i+1) > 0.01:
            for k, v in sorted(counts.items()):
                logger.info(f"  {k}: {v}")
            break
        counts[idType] += 1


if __name__ == "__main__":
    main(sys.argv[1:])
