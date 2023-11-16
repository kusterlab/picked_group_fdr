import _csv
import sys
import os
import collections
import logging
from typing import List

from .. import __version__, __copyright__
from .. import helpers
from ..picked_group_fdr import ArgumentParserWithLogger
from ..parsers import tsv
from ..parsers import modifications
from ..parsers import percolator
from ..parsers import maxquant

# hacky way to get the package logger instead of just __main__ when running as python -m picked_group_fdr.pipeline.update_evidence_from_pout ...
logger = logging.getLogger(__package__ + "." + __file__)


def parseArgs(argv):
    import argparse
    apars = ArgumentParserWithLogger(
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    apars.add_argument('--mq_evidence', default=None, metavar = "EV", nargs='+', required = True,
                       help='''MaxQuant evidence.txt or msms.txt file(s). If you 
                               want to combine multiple evidence files, use 
                               spaces to separate the file paths.''')
                                                    
    apars.add_argument('--mq_evidence_out', default=None, metavar = "EV", required = True,
                       help='''MaxQuant evidence.txt or msms.txt combined output 
                               file with updated percolator results.''')
    
    apars.add_argument('--perc_results', default=list(), metavar = "POUT", nargs='+',
                       help='''Percolator output file(s) with PSMs separated by 
                               spaces. Also include the decoy results. If this 
                               flag is not set, the evidence files are simply 
                               concatenated without updating the PSMs.''')
                                                    
    apars.add_argument('--mq_msms', default=None, metavar = "M", nargs='+',
                       help='''MaxQuant msms.txt file(s) (optional). Can help 
                               resolve MS1 features with multiple scans. If you
                               want to combine multiple evidence files, use 
                               spaces to separate the file paths. Not 
                               implemented yet...''')
    
    apars.add_argument('--pout_input_type', default="andromeda", metavar = "SE",
                       help='''Input type of percolator output file. Can be 
                               "andromeda" or "prosit".''')
    
    apars.add_argument('--mq_input_type', default="auto", metavar = "IT",
                       help='''Input type of MaxQuant output file. Can be 
                               "peptides", "evidence", "msms" or "auto" 
                               (automatic detection).''')
                                                    
    # ------------------------------------------------
    args = apars.parse_args(argv)
    
    return args


def main(argv):
    logger.info(f'UpdateEvidence version {__version__}\n{__copyright__}')
    logger.info(f'Issued command: {os.path.basename(__file__)} {" ".join(map(str, argv))}')
    
    args = parseArgs(argv)
    
    if os.path.isfile(args.mq_evidence_out):
        logger.info(f"Found updated evidence file {args.mq_evidence_out}, remove this file to rerun update_evidence.")
        return
    
    if args.mq_input_type == "peptides":
        updatePeptides(args.mq_evidence, args.perc_results, args.mq_evidence_out, args.mq_msms, args.pout_input_type)
    else:
        updateEvidence(args.mq_evidence, args.perc_results, args.mq_evidence_out, args.mq_msms, args.pout_input_type)
    os.rename(args.mq_evidence_out + ".tmp", args.mq_evidence_out)

    
def updateEvidence(evidenceFiles, poutFiles, outEvidenceFile, msmsFiles, poutInputType):
    fixed_mods, resultsDict = get_percolator_results(poutFiles, poutInputType)
    
    logger.info("Writing updated combined evidence file")
    writer = tsv.get_tsv_writer(outEvidenceFile + ".tmp")
    firstHeaders = []
    for evidenceFile in evidenceFiles:
        firstHeaders = update_evidence_single(evidenceFile, writer, firstHeaders, fixed_mods, resultsDict, poutInputType)

    logger.info(f"Results written to {outEvidenceFile}")


def update_evidence_single(evidenceFile: str, writer: '_csv._writer', first_headers: List[str], fixed_mods, resultsDict, poutInputType: str):
    logger.info(f"Processing {evidenceFile}")
    reader = tsv.get_tsv_reader(evidenceFile)

    headersOrig = next(reader) # save the header
    headers = list(map(lambda x : x.lower(), headersOrig))
    if len(first_headers) == 0:
        writer.writerow(headersOrig)
    elif headers != first_headers:
        logger.warning("Current column names are different from the first evidence file")
        logger.info('Column\tFirstHeaders\tCurrentHeaders')
        logger.info('\n'.join([str(i+1) + '\t' + x + '\t' + y for i, (x, y) in enumerate(zip(first_headers, headers)) if x != y]))

    # these columns will be updated
    scoreCol = tsv.get_column_index(headers, 'score')
    postErrProbCol = tsv.get_column_index(headers, 'pep')

    mqPEPs, prositPEPs = list(), list()
    unexplainedMissingPSMs = 0
    unexplainedPeptides = list()
    writtenRows = 0
    missingRawFiles = set()
    for row, psm in maxquant.parse_evidence_file_for_percolator_matching(reader, headers):
        if len(resultsDict) > 0 and psm.scanNr >= 0:
            if len(resultsDict[psm.rawFile]) == 0 and psm.rawFile not in missingRawFiles:
                logger.warning(f"Could not find any PSMs for raw file {psm.rawFile} in the percolator result files")
                missingRawFiles.add(psm.rawFile)
                continue
            
            percResult, peptide = find_percolator_psm(psm, fixed_mods, resultsDict, poutInputType)
            
            if not psm.isDecoy and psm.rawFile in resultsDict:
                mqPEPs.append((psm.postErrProb, psm.idType))
            
            if not percResult:
                # check for unexplainable missing PSMs
                if not psm.isContaminant and psm.score > 0.0 and (poutInputType != "prosit" or is_valid_prosit_peptide(peptide)):
                    unexplainedMissingPSMs += 1
                    logger.debug(f"Unexplained missing PSM in percolator output: {psm.rawFile}, {peptide}, {psm.scanNr}")
                    if unexplainedMissingPSMs <= 10:
                        unexplainedPeptides.append(psm.peptideOriginal)
                logger.debug(f"Missing PSM in percolator output: {psm.rawFile}, {peptide}, {psm.scanNr}")
                continue
            
            if not psm.isDecoy:
                prositPEPs.append((percResult[1], psm.idType))

            row[scoreCol] = percResult[0]
            row[postErrProbCol] = percResult[1]
            
        if writtenRows % 500000 == 0:
            logger.info(f"    Writing line {writtenRows}")
        writtenRows += 1
        
        writer.writerow(row)
    
    unexplainedPercentage = int(unexplainedMissingPSMs/writtenRows*100)
    logger.info(f"Unexplained missing PSMs in Percolator results: {unexplainedMissingPSMs} out of {writtenRows} ({unexplainedPercentage}%)")
    if unexplainedMissingPSMs > 0:
        logger.info("If this percentage is low (<5%), it is probably the result of second peptide hits.")
        logger.info("\tFirst 10 missing peptides:")
        for peptide in unexplainedPeptides:
            logger.info("\t" + peptide)
    
    logger.info("#MQ identifications:")
    countBelowFDR(mqPEPs)
    
    logger.info(f"#{poutInputType} identifications:")
    countBelowFDR(prositPEPs)
    
    return headers


def find_percolator_psm(psm: maxquant.EvidenceRow, fixed_mods, resultsDict, poutInputType: str):
    peptide = psm.peptideOriginal[1:-1]
    if poutInputType == "prosit":
        fixed_mods_tmp = fixed_mods
        if maxquant.is_heavy_labeled(psm.labelingState):
            fixed_mods_tmp = modifications.SILAC_HEAVY_FIXED_MODS
        peptide = modifications.maxquant_mod_to_unimod([psm.peptideOriginal], fixed_mods=fixed_mods_tmp)[0]

    percResult = resultsDict[psm.rawFile].get((psm.scanNr, peptide), None)
    if poutInputType == "prosit" and not percResult and maxquant.has_unknown_silac_label(psm.labelingState):
        # if scanNr is not found, retry with heavy labeling for SILAC because 
        # labelingState column is not reliable when multiple MS/MS map to the precursor
        fixed_mods_tmp = modifications.SILAC_HEAVY_FIXED_MODS
        peptide = modifications.maxquant_mod_to_unimod([psm.peptideOriginal], fixed_mods=fixed_mods_tmp)[0]
        percResult = resultsDict[psm.rawFile].get((psm.scanNr, peptide), None)
    return percResult, peptide


def updatePeptides(peptideFiles, poutFiles, outPeptideFile, msmsFiles, poutInputType):
    _, resultsDict = get_percolator_results(poutFiles, poutInputType)
        
    peptideResultsDict = convertPSMDictToPeptideDict(resultsDict)
    
    logger.info("Writing updated peptides file")
    writer = tsv.get_tsv_writer(outPeptideFile + ".tmp")
    first = True
    mqPEPs, prositPEPs = list(), list()
    unexplainedMissingPeptides = 0
    unexplainedPeptides = list()
    for peptideFile in peptideFiles:
        logger.info(f"Processing {peptideFile}")
        reader = tsv.get_tsv_reader(peptideFile)
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
        
        scoreCol = tsv.get_column_index(headers, 'score')
        postErrProbCol = tsv.get_column_index(headers, 'pep')
        
        peptCol = tsv.get_column_index(headers, 'sequence')
        reverseCol = tsv.get_column_index(headers, 'reverse')
        
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


def get_percolator_results(pout_files: List[str], pout_input_type: str):
    results_dict = collections.defaultdict(dict)
    fixed_mods = None
    for pout_file in pout_files:
        logger.info(f"Processing {pout_file}")
        fixed_mods, results_dict = percolator.parse_percolator_out_file_to_dict(pout_file, results_dict, pout_input_type)
    
    logger.info("Finished parsing percolator output files")
    logger.info("#PSMs per raw file:")
    for raw_file, psms in results_dict.items():
        logger.info(f"  {raw_file}: {len(psms)} PSMs")
    return fixed_mods, results_dict


def convertPSMDictToPeptideDict(resultsDict):
    peptideResultsDict = dict()
    for _, results in resultsDict.items():
        for (_, peptide), (score, PEP) in results.items():
            peptide = helpers.clean_peptide(peptide, removeFlanks = False)
            currScore, currPEP = peptideResultsDict.get(peptide, (-1e10, 1e10))
            peptideResultsDict[peptide] = (max(currScore, score), min(currPEP, PEP))
    return peptideResultsDict


def is_valid_prosit_peptide(peptide):
    return (len(peptide) <= 30 and not "(ac)" in peptide)


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
