import sys
import os
import csv
import collections
import logging

import numpy as np
import triqler.qvality as qvality

from .. import __version__, __copyright__
from .. import helpers
from .. import digest
from .. import parsers
from ..picked_group_fdr import ArgumentParserWithLogger

# hacky way to get the package logger instead of just __main__ when running as python -m picked_group_fdr.pipeline.update_evidence_from_pout ...
logger = logging.getLogger(__package__ + "." + __file__)


def parseArgs(argv):
    import argparse
    apars = ArgumentParserWithLogger(
            description='''Merges percolator output files to a list of ranked 
                           peptides, again in percolator output format.
                           Peptides are remapped to their protein sequences
                           based on an in-silico digestion of the fasta database
                           Peptides are then ranked based on their q-values.
                           For target peptides, q-values are re-estimated from
                           the posterior error probabilities to obtain a
                           higher resolution in the high-confident region.''',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    apars.add_argument('--perc_results', default=list(), metavar = "POUT", nargs='+',
                       help='''Percolator output file(s), separated by 
                               spaces. Also include the decoy results!''')
    
    apars.add_argument('--fasta', default=None, metavar = "F", required = True,
                       help='''Fasta file with protein sequences''')

    apars.add_argument('--perc_merged', default="merged_pout.txt", metavar = "M",
                       help='''Output file with merged pout results.''')
    
    digest.addArguments(apars)
                                              
    # ------------------------------------------------
    args = apars.parse_args(argv)
    
    return args


def main(argv):
    logger.info(f'MergePout version {__version__}\n{__copyright__}')
    logger.info(f'Issued command: {os.path.basename(__file__)} {" ".join(map(str, argv))}')
    
    args = parseArgs(argv)
    
    if os.path.isfile(args.perc_merged):
        logger.info(f"Found merged pout file {args.perc_merged}, remove this file to rerun merge_pout.")
        return
    
    peptideToProteinMap = digest.getPeptideToProteinMapWithEnzyme(
        args.fasta, args.min_length, args.max_length, args.enzyme, 
        args.cleavages, list(args.special_aas), db = "concat")

    merge_pout(args.perc_results, peptideToProteinMap, args.perc_merged)


def merge_pout(perc_results, peptideToProteinMap, perc_merged):
    seenPeptides = dict()
    missingPeptides, matchedPeptides = 0, 0
    for poutFile in perc_results:
        poutReader = parsers.getTsvReader(poutFile)
        headers = next(poutReader)
        
        idCol, peptCol, scoreCol, qvalCol, postErrProbCol, _ = parsers.getPercolatorColumnIdxs(headers)

        sumPEP = 0.0
        matchedBefore = matchedPeptides
        logger.info(f"Parsing {poutFile}")
        for i, row in enumerate(poutReader):
            if i % 1000000 == 0:
                logger.info(f"Processing row {i}: #Missing peptides: {missingPeptides}, #Matched peptides: {matchedPeptides}")
            
            # removeFlanks=True only removes a single character (MaxQuant convention)
            # convert peptide string to upper case, since prosit converts modified amino acids to lower case
            peptide = helpers.cleanPeptide(row[peptCol][2:-2].upper(), removeFlanks=False)
            proteins = digest.getProteins(peptideToProteinMap, peptide)
            isDecoy = helpers.isDecoy(proteins)
            
            if isDecoy:
                qValue = float(row[qvalCol])
            else:
                # for targets, estimate "high resolution" q-value based on PEPs
                sumPEP += float(row[postErrProbCol])
                qValue = sumPEP / (i+1)

            if qValue < seenPeptides.get(peptide, (1.0, []))[0]:
                if len(proteins) > 0:
                    row = [row[idCol], row[scoreCol], row[qvalCol], row[postErrProbCol], "-." + peptide + ".-"] + proteins
                    matchedPeptides += 1                
                    seenPeptides[peptide] = (qValue, row, isDecoy)
                else:
                    if not helpers.isContaminant(proteins):
                        logger.debug(f"Could not find peptide {peptide} in fasta file, check your database and if the correct digestion parameters were specified")
                    missingPeptides += 1
        
        logger.info(f"Processing row {i}: #Missing peptides: {missingPeptides}, #Matched peptides: {matchedPeptides}")
        
        if matchedBefore == matchedPeptides:
            logger.warning(f"No new peptides added by {poutFile}")
    
    psm_infos = sorted(seenPeptides.values())
    peps = get_peptide_PEPs(psm_infos)
    
    write_updated_PSMs(perc_merged, psm_infos, peps, update_qvals=False)
    

def get_peptide_PEPs(psm_infos):
    targetScores, decoyScores = list(), list()
    for qval, _, isDecoy in psm_infos:
        if isDecoy:
            decoyScores.append(-1*np.log10(qval))
        else:
            targetScores.append(-1*np.log10(qval))
    
    _, peps = qvality.getQvaluesFromScores(targetScores, decoyScores, includeDecoys=True, plotRegressionCurve=False)
    return peps


def write_updated_PSMs(perc_merged, psm_infos, peps, update_qvals=False):
    writer = parsers.getTsvWriter(perc_merged + ".tmp")
    
    headers = parsers.getPercolatorNativeHeaders()
    writer.writerow(headers)
    _, _, _, qvalCol, _, _ = parsers.getPercolatorColumnIdxs(headers)
    
    sumPEP = 0.0
    counts = 0
    for (_, row, isDecoy), pep in zip(psm_infos, peps):
        if update_qvals:
            if not isDecoy:
                sumPEP += pep
                counts += 1
            row[qvalCol] = sumPEP / counts
        writer.writerow(row)
    
    os.rename(perc_merged + ".tmp", perc_merged)


if __name__ == "__main__":
    main(sys.argv[1:])
