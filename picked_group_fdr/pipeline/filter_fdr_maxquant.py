import sys
import collections
import re

import numpy as np

from .. import parsers
from .. import fdr


Scan = collections.namedtuple('Scan', 'postErrorProb rawFile scanNr sequence precCharge isDecoy')


def parseArgs(argv):
  import argparse
  apars = argparse.ArgumentParser(
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  apars.add_argument('--mq_msms', default=None, metavar="M", nargs='+', required=False,
                     help='MaxQuant msms.txt or evidence.txt file(s).')
  
  apars.add_argument('--mq_msms_out', default=None, metavar="O", required=False,
                     help='Output path for filtered msms.txt or evidence.txt file.')
  
  apars.add_argument('--mq_protein_groups', default=None, metavar="PG", nargs='+', required=False,
                     help='proteinGroups.txt filtered at 100% FDR file(s)')
  
  apars.add_argument('--mq_protein_groups_out', default=None, metavar="O", required=False,
                     help='Output path for filtered msms.txt or evidence.txt file.')
  
  apars.add_argument('--fdr_cutoff', default=0.01, metavar="C", type = float,
                     help='FDR threshold.')

  apars.add_argument('--psm_level_fdr',
                     help='Use PSM-level FDR cutoff instead of the default peptide-level FDR cutoff.',
                     action='store_true')
  
  apars.add_argument('--precursor_level_fdr',
                     help='Use precursor-level FDR cutoff instead of the default peptide-level FDR cutoff.',
                     action='store_true')

  apars.add_argument('--per_rawfile_fdr',
                     help='Apply the FDR cutoff per raw file, instead of globally.',
                     action='store_true')
   
  # ------------------------------------------------
  args = apars.parse_args(argv)
  
  return args

 
def main(argv):
  args = parseArgs(argv)
  
  if args.mq_msms:
    if args.mq_msms_out:
      filterAtFDR(args.mq_msms, args.mq_msms_out, args.fdr_cutoff, args.psm_level_fdr, args.precursor_level_fdr, args.per_rawfile_fdr)
    else:
      print("ERROR: if --mq_msms is set, --mq_msms_out also needs to be set")  
  
  if args.mq_protein_groups:
    if args.mq_protein_groups_out:
      filterProteinGroupsAtFDR(args.mq_protein_groups, args.mq_protein_groups_out, args.fdr_cutoff)
    else:
      print("ERROR: if --mq_protein_groups is set, --mq_protein_groups_out also needs to be set")


def filterProteinGroupsAtFDR(proteinGroupFiles, outputFile, fdrCutoff):
  for proteinGroupFile in proteinGroupFiles:
    if ".txt" not in proteinGroupFile:
      sys.exit("ERROR: could not detect file extension .txt in the input file")
    
    print("Writing filtered protein groups results to:", outputFile)
    writer = parsers.getTsvWriter(outputFile)
    reader = parsers.getTsvReader(proteinGroupFile)
    header = next(reader)
    writer.writerow(header)
    
    qvalCol = header.index('Q-value')
    for row in reader:
      if float(row[qvalCol]) <= fdrCutoff:
        writer.writerow(row)
    
    print("Finished writing")


def filterAtFDR(msmsFiles, msmsOutputFile, fdrCutoff, psmLevelFDR, precursorLevelFDR, perRawFileFDR):  
  level = 'peptide'
  if precursorLevelFDR:
    level = 'precursor'
  if psmLevelFDR:
    level = 'PSM'
  
  bestScans = collections.defaultdict(lambda: collections.defaultdict(lambda: Scan(np.inf, "", "", "", "", "")))
  rawFiles = set()
  for msmsFile in msmsFiles:
    print("Processing", msmsFile)
    reader = parsers.getTsvReader(msmsFile)
    headersOrig = next(reader) # save the header
    headers = list(map(lambda x : x.lower(), headersOrig))
    
    postErrorProbCol = headers.index('pep')
    seqCol = headers.index('sequence')
    reverseCol = headers.index('reverse')
    
    rawFileCol = None
    if 'raw file' in headers and perRawFileFDR:
      rawFileCol = headers.index('raw file')
      
    scanNrCol = headers.index('sequence')
    if 'scan number' in headers:
      scanNrCol = headers.index('scan number')
    
    chargeCol = None
    if 'charge' in headers:
      chargeCol = headers.index('charge')
    
    for i, row in enumerate(reader):
      if i % 500000 == 0:
        print("Processing line %d" % i)
      
      isDecoy = False
      if row[reverseCol] == '+':
        isDecoy = True
      
      postErrorProb = float(row[postErrorProbCol])
      sequence = row[seqCol]
      precCharge = int(row[chargeCol]) if chargeCol else 2
      rawFile = row[rawFileCol] if rawFileCol else ""
      scanNr = row[scanNrCol]
      
      key1 = sequence
      if precursorLevelFDR:
        key1 = (sequence, precCharge)
      elif psmLevelFDR:
        key1 = scanNr
      
      key = key1
      if perRawFileFDR:
        key = (key1, rawFile)
      rawFiles.add(rawFile)
      
      if not np.isnan(postErrorProb) and postErrorProb < bestScans[rawFile][key].postErrorProb:
        bestScans[rawFile][key] = Scan(postErrorProb, rawFile, scanNr, sequence, precCharge, isDecoy)
  
  print(f"Filtering best scan per {level}")
  
  survivingPeptides, postErrorProbCutoffs = getPeptidesBelowFdr(rawFiles, bestScans, fdrCutoff, level)
  writeMsmsOutput(msmsOutputFile, headersOrig, msmsFiles, postErrorProbCutoffs, survivingPeptides, perRawFileFDR)


def getPeptidesBelowFdr(rawFiles, bestScans, fdrCutoff, level):  
  survivingPeptides = list()
  postErrorProbCutoffs = dict()
  for rawFile in rawFiles:
    survivingScans = set(bestScans[rawFile].values())
    sortedPostErrorProbs = sorted([(x.postErrorProb, x.isDecoy) for x in survivingScans])
    postErrorProbCutoff = calculatePostErrProbCutoff(sortedPostErrorProbs, fdrCutoff)
    postErrorProbCutoffs[rawFile] = postErrorProbCutoff
    survivingPeptides.extend([(x.sequence, x.rawFile) for x in bestScans[rawFile].values() if x.postErrorProb <= postErrorProbCutoff])
    print(f"PEP at 1% {level}-level FDR for {rawFile}: {postErrorProbCutoff}")
  survivingPeptides = set(survivingPeptides)
  
  return survivingPeptides, postErrorProbCutoffs


def writeMsmsOutput(msmsOutputFile, headersOrig, msmsFiles, postErrorProbCutoffs, survivingPeptides, perRawFileFDR):
  print("Writing msms.txt output file")
  writer = parsers.getTsvWriter(msmsOutputFile)
  writer.writerow(headersOrig)
  
  survivingPSMs = 0
  for msmsFile in msmsFiles:
    print("Processing", msmsFile)
    reader = parsers.getTsvReader(msmsFile)
    headersOrig = next(reader)
    headers = list(map(lambda x : x.lower(), headersOrig))
    
    if 'raw file' in headers and perRawFileFDR:
      rawFileCol = headers.index('raw file')
    else:
      rawFileCol = None
    seqCol = headers.index('sequence')
    postErrorProbCol = headers.index('pep')
    for i, row in enumerate(reader):
      if i % 500000 == 0:
        print("Processing line %d" % i)
      
      postErrorProb = float(row[postErrorProbCol])
      sequence = row[seqCol]
      rawFile = row[rawFileCol] if rawFileCol else ""
      
      key = (sequence, rawFile)
      postErrorProbCutoff = postErrorProbCutoffs[rawFile]
      
      # this PEP comparison lets NaN values set by MBR pass on purpose. Therefore, we also need to check if the peptide actually has been identified in at least one run
      if postErrorProb > postErrorProbCutoff or key not in survivingPeptides:
        continue
      
      survivingPSMs += 1
      writer.writerow(row)
  
  print("Written {} rows".format(survivingPSMs))


def calculatePostErrProbCutoff(sortedPostErrorProbs, fdrCutoff):
  # both the PSM and peptide-level thresholds calculated here seem to be more conservative than MaxQuant's approach
  tps, fps = 1, 1
  fdrs, peps = list(), list()
  for postErrProb, isDecoy in sortedPostErrorProbs:
    if isDecoy:
      fps += 1
    else:
      tps += 1
    fdrs.append(fps / tps)
    peps.append(postErrProb)
  
  qvals = fdr.fdrsToQvals(fdrs)
  if qvals[-1] <= fdrCutoff:
    postErrorProbCutoff = peps[-1]
  else:
    postErrorProbCutoff = next(x[0] for x in zip(peps, qvals) if x[1] > fdrCutoff)
  
  return postErrorProbCutoff


if __name__ == "__main__":
  main(sys.argv[1:])
