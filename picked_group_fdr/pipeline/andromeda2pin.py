#!/usr/bin/python

'''
Converts Andromeda evidence.txt output to tab delimited percolator input file
'''

import sys
import os
import csv
import getopt
import numpy as np

from .. import digest
from .. import parsers

# TODO allow mqpar.xml as input
def main(argv):
  args = parseArgs(argv)
  
  andromedaTargetOutFNs = []
  if len(args.mq_evidence_file) > 0:
    with open(args.mq_evidence_file, 'r') as f:
      for line in f:
        if line.lower().startswith("sequence"):
          andromedaTargetOutFNs.append(args.mq_evidence_file)
          break
        else:
          if len(andromedaTargetOutFNs) == 0:
            print("Meta file detected, interpreting each line as a path to an evidence file")
          andromedaTargetOutFNs.append(line.rstrip())
  else:
    sys.exit("No input files provided")
  
  percInFN = args.outputTab
  decoyPattern = args.pattern
  #numHits = args.matches
  numHits = 1
  fastaFile = args.databases
  
  if len(percInFN) > 0:
    if os.path.isfile(percInFN):
      print(f"Found output file {percInFN}, remove this file to re-run andromeda2pin.")
      return
    else:
      print("Writing results to:", percInFN)
      writer = parsers.getTsvWriter(percInFN)
  else:
    print("Writing results to stdout")
    writer = parsers.getTsvWriter(sys.stdout)
  
  charges = list(range(2,7))
  peptideToProteinMap = getPeptideToProteinMap(fastaFile, args.min_length, args.max_length, args.enzyme, args.cleavages, list(args.special_aas), db = "concat")
  
  writeHeaders(writer, charges)
  
  for andromedaTargetOutFN in andromedaTargetOutFNs:
    convertAndromedaOutToPin(andromedaTargetOutFN, writer, charges, numHits, peptideToProteinMap, decoyPattern = decoyPattern)
  
  print("Finished writing percolator input")

def parseArgs(argv):
  import argparse
  apars = argparse.ArgumentParser(
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  apars.add_argument('mq_evidence_file', default=None, metavar = "evidence.txt",
   help='''MaxQuant evidence file, or a meta file if you have multiple evidence 
           files you want to combine. Meta files are text files containing the
           the paths of mzid-files, one path per line. For successful results,
           the different runs should be generated under similar conditions.
        ''')
                          
  apars.add_argument('-o', '--outputTab', default = None, metavar='pin.tab', 
                     help='''Save output in a tab delimited file
                          ''')
  
  #apars.add_argument('-m', '--matches', default = 1, metavar='M', type=int, 
  #                   help='''Maximal number of matches to take in consideration 
  #                           per spectrum
  #                        ''')
  apars.add_argument('-P', '--pattern', default = "REV__", metavar='P',
                     help='''Pattern used to identify the decoy PSMs.
                          ''')
  
  apars.add_argument('-F', '--databases', default = None, metavar='F',
                     help='''Fasta database used in the search 
                             against the spectra file.
                          ''')
  
  digest.addArguments(apars)
                          
  # ------------------------------------------------
  args = apars.parse_args(argv)
  
  return args
  
def writeHeaders(writer, charges):
  writer.writerow(["SpecId", "Label", "ScanNr", "ExpMass", "AndromedaScore", "DeltaScore", "PepLen"] + ["Charge" + str(i) for i in charges] + ["Mass", "enzN", "enzC", "enzInt", "numMods", "dM", "absdM", "Peptide", "Proteins"])
  #writer.writerow(["DefaultDirection", "-", "-", "-", 1, 0.5, 0] + [0 for i in charges] + [0, 0, 0, -1.5, -2, 0, -1])

def parseMqEvidenceFile(mqEvidenceFile, razor = False):
  if mqEvidenceFile.endswith('.csv'):
    delimiter = ','
  else:
    delimiter = '\t'
  reader = parsers.getTsvReader(mqEvidenceFile, delimiter)
  headers = next(reader) # save the header
  headers = list(map(lambda x : x.lower(), headers))
  
  if mqEvidenceFile.endswith('.csv'):
    headers = [x.replace(".", " ") for x in headers]
  
  peptCol = headers.index('modified sequence')
  idCol = headers.index('ms/ms scan number')
  fileCol = headers.index('raw file')
  chargeCol = headers.index('charge')
  massCol = headers.index('mass')
  if 'mass error [ppm]' in headers:
    deltaMassCol = headers.index('mass error [ppm]')
  else:
    deltaMassCol = -1
  #intensityCol = headers.index('intensity')
  if 'leading proteins' in headers:
    proteinCol = headers.index('leading proteins') # all proteins the peptide matches to, does not return REV__ proteins
  else:
    proteinCol = headers.index('proteins')
  scoreCol = headers.index('score')
  deltaScoreCol = headers.index('delta score')
  #postErrCol = headers.index('pep')
  #rtCol = headers.index('retention time')
  
  #fractionCol = headers.index('fraction') if 'Fraction' in headers else -1
  if 'experiment' in headers:
    experimentCol = headers.index('experiment')
  else:
    experimentCol = -1
  
  print("Parsing MaxQuant evidence.txt file")
  quants = list()
  for lineIdx, row in enumerate(reader):
    if lineIdx % 500000 == 0:
      print("  Reading line", lineIdx)
    
    if len(row[idCol]) == 0:
      continue
    
    scanNr = int(row[idCol])
    charge = int(row[chargeCol])
    fileName = row[fileCol]
    peptide = "-." + row[peptCol][1:-1].replace('pS', 'S[80]').replace('pT', 'T[80]').replace('pY', 'Y[80]').replace('(ph)', '[80]').replace('(ox)', '[16]').replace('(ac)', '[42]').replace('(Acetyl (Protein N-term))', '[42]').replace('(Oxidation (M))', '[16]') + ".-"
    proteins = row[proteinCol].split(";")
    if experimentCol >= 0:
      experiment = row[experimentCol]
    else:
      experiment = "Experiment1"

    score = float(row[scoreCol])
    deltaScore = float(row[deltaScoreCol])
    mass = float(row[massCol])
    deltaMass = float(row[deltaMassCol])
    
    if np.isnan(deltaMass):
      deltaMass = 0.0
    
    if not np.isnan(score) and score > 0.0:
      yield scanNr, charge, fileName, peptide, proteins, experiment, score, deltaScore, mass, deltaMass

def getPeptideToProteinMap(fastaFile, min_len, max_len, enzyme, miscleavages, specialAAs, db):
  pre, not_post = digest.getCleavageSites(enzyme)
  
  if len(fastaFile) > 0:
    return digest.getPeptideToProteinMap(fastaFile, db, digestion = 'full', min_len = min_len, max_len = max_len, pre = pre, not_post = not_post, miscleavages = miscleavages, methionineCleavage = True, specialAAs = specialAAs)
  else:
    return dict()
  
def convertAndromedaOutToPin(andromedaOutFN, writer, charges, numHits, peptideToProteinMap, decoyPattern = ""):
  print("Reading", andromedaOutFN)
  
  rows = list()
  for scanNr, charge, fileName, peptide, tmp_proteins, experiment, andromedaScore, deltaScore, expMass, deltaMass in parseMqEvidenceFile(andromedaOutFN):
    rank = 1
    psmId = fileName + "_" + str(scanNr) + "_" + str(charge) + "_" + str(rank)
    modPeptide, cleanPeptide, pepLen, enzN, enzC, enzInt, numMods = getPeptideStats(peptide, deltaMass)
    absDeltaMass = abs(deltaMass)
    
    if pepLen >= 6:
      if len(peptideToProteinMap) > 0:
        proteins = digest.getProteins(peptideToProteinMap, cleanPeptide[2:-2])
        if len(proteins) == 0:
          if not "CON__" in tmp_proteins[0]:
            print("WARNING: Could not find peptide " + peptide + " " + str(tmp_proteins) + " in fasta database, skipping PSM")
          continue
      
      if len(decoyPattern) > 0:
        if sum(1 for p in proteins if p.startswith(decoyPattern)) == len(proteins):
          label = -1
        else:
          label = 1
      
      r = [psmId, label, scanNr, expMass, andromedaScore, deltaScore, pepLen] + [1 if charge == i else 0 for i in charges] + [expMass, enzN, enzC, enzInt, numMods, deltaMass, absDeltaMass, modPeptide] + proteins
      writer.writerow(r)

def getPeptideStats(peptide, deltaMass, accurateModMasses = False):
  cleanPeptide = peptide[:2]
  modPeptide = peptide[:2]
  numMods, enzInt = 0, 0
  aaIdx = 2
  while aaIdx < len(peptide):
    if peptide[aaIdx] == "[":
      modStart = aaIdx
      while peptide[aaIdx+1].isdigit():
        aaIdx += 1
      mod = "["
      
      multFactor = 1
      if peptide[modStart] == "-":
        multFactor = -1
      modMass = multFactor * float(peptide[modStart+1:aaIdx+1])
          
      # modification masses are typically not integers, add the deltamass to the first mod to ensure calcmass = expmass
      if accurateModMasses and numMods == 0:
        modMass += deltaMass
        mod += str(round(modMass, 4))
      else:
        mod += str(int(modMass))
      
      mod += "]"
      aaIdx += 1
      modPeptide += mod
      numMods += 1
    elif peptide[aaIdx] == ".":
      modPeptide += peptide[aaIdx:]
      cleanPeptide += peptide[aaIdx:]
      break
    else:
      modPeptide += peptide[aaIdx]
      cleanPeptide += peptide[aaIdx]
      if digest.isEnzymatic(cleanPeptide[-2], cleanPeptide[-1], not_post = [], methionineCleavage = False):
        enzInt += 1
    aaIdx += 1
  enzN = int(digest.isEnzymatic(cleanPeptide[0], cleanPeptide[2], not_post = [])) # Andromeda uses Trypsin/P
  enzC = int(digest.isEnzymatic(cleanPeptide[-3], cleanPeptide[-1], not_post = []))
  pepLen = len(cleanPeptide) - 4
  return modPeptide, cleanPeptide, pepLen, enzN, enzC, enzInt, numMods

if __name__ == "__main__":
   main(sys.argv[1:])

