import sys
import csv
import itertools
import collections

import numpy as np
import scipy.stats


def main(argv):  
  proteotypicityFile = argv[0]
  peptideToProteinMapFile = argv[1]
  simulatedEvidenceFile = argv[2]
  simulatedProteinGroups = simulatedEvidenceFile.replace('/evidence.', '/proteinGroups.')
  seed = int(argv[3])
  np.random.seed(seed)
  
  proteinToGeneMapFile = ""
  if len(argv) > 4:
    proteinToGeneMapFile = argv[4]
  
  numExperiments = 10
  numProteinsMean = 10000
  numProteinsStdev = 1000
  truePosPeptideScoresMean = 2.5
  truePosPeptideScoresStdev = 0.7
  falsePosPeptideScoresMean = 0.0
  falsePosPeptideScoresStdev = 0.7
  peptideFdr = 0.01
  incorrectRatio = 0.6
  recalculatePeptFDR = False
  
  with open(simulatedEvidenceFile + '.params.txt', 'w') as f:
    f.write('numExperiments = {}\n'.format(numExperiments))
    f.write('numProteinsMean = {}\n'.format(numProteinsMean))
    f.write('numProteinsStdev = {}\n'.format(numProteinsStdev))
    f.write('truePosPeptideScoresMean = {}\n'.format(truePosPeptideScoresMean))
    f.write('truePosPeptideScoresStdev = {}\n'.format(truePosPeptideScoresStdev))
    f.write('falsePosPeptideScoresMean = {}\n'.format(falsePosPeptideScoresMean))
    f.write('falsePosPeptideScoresStdev = {}\n'.format(falsePosPeptideScoresStdev))
    f.write('peptideFdr = {}\n'.format(peptideFdr))
    f.write('incorrectRatio = {}\n'.format(incorrectRatio))
    f.write('recalculatePeptFDR = {}\n'.format(recalculatePeptFDR))
    f.write('proteotypicityFile = {}\n'.format(proteotypicityFile))
    f.write('peptideToProteinMapFile = {}\n'.format(peptideToProteinMapFile))
    f.write('seed = {}\n'.format(seed))
  
  pi0 = 1.0 - incorrectRatio
  
  minScore = scipy.stats.norm.ppf(1 - pi0 * peptideFdr / incorrectRatio, falsePosPeptideScoresMean, falsePosPeptideScoresStdev)
  print("Score threshold:", minScore)
  
  #plotScoreDistributions(truePosPeptideScoresMean, truePosPeptideScoresStdev, falsePosPeptideScoresMean, falsePosPeptideScoresStdev, minScore, incorrectRatio)
  
  minScoreTruePos = (minScore - truePosPeptideScoresMean) / truePosPeptideScoresStdev
  minScoreFalsePos = (minScore - falsePosPeptideScoresMean) / falsePosPeptideScoresStdev
  
  peptideToProteinMap = getPeptideToProteinMap(peptideToProteinMapFile)
  
  proteinToGeneMap = dict()
  if len(proteinToGeneMapFile) > 0:
    proteinToGeneMap, geneToProteinMap = getProteinToGeneMap(proteinToGeneMapFile)
  
  proteinProbs, peptideProbsPerProtein, geneProbs = getProbabilities(proteotypicityFile, peptideToProteinMap, proteinToGeneMap)
  
  allPeptides = list(peptideToProteinMap.keys())
  
  scoresCombined = list()
  identifiedPeptidesCombined = list()
  presentProteinsCombined = list()
  peptideLabelsCombined = list()
  for experimentId in range(numExperiments): 
    print("Simulating experiment:", experimentId)
       
    numProteins = np.random.normal(numProteinsMean, numProteinsStdev, 1)
    if len(proteinToGeneMapFile) > 0:
      presentProteins = list()
      presentGenes = np.random.choice(list(geneProbs.keys()), max(200, int(numProteins[0])), p = list(geneProbs.values()))
      for gene in presentGenes:
        if len(geneToProteinMap[gene]) == 1:
          presentProteins.append(geneToProteinMap[gene][0])
        else:
          proteinsForGene, proteinProbsForGene = zip(*[(p, proteinProbs.get(p,0)) for p in geneToProteinMap[gene]])
          proteinProbsForGene = [x / sum(proteinProbsForGene) for x in proteinProbsForGene]
          chosenProtein = np.random.choice(proteinsForGene, 1, p = proteinProbsForGene)[0]
          presentProteins.append(chosenProtein)
    else:
      presentProteins = np.random.choice(list(proteinProbs.keys()), max(200, int(numProteins[0])), replace = False, p = list(proteinProbs.values()))
      
    presentProteinsCombined.extend(presentProteins)
    
    identifiedPeptides = list()
    for protein in presentProteins:
      peptides, peptideProbs = zip(*peptideProbsPerProtein[protein])

      identifiedPeptidesForProtein = []
      while len(identifiedPeptidesForProtein) == 0:
        uniformRandomDraws = np.random.rand(len(peptides))
        identifiedPeptidesForProtein = [p for p, prob, draw in zip(peptides, peptideProbs, uniformRandomDraws) if draw < prob]
      identifiedPeptides.extend(identifiedPeptidesForProtein)
    
    scores = list(scipy.stats.truncnorm.rvs(minScoreTruePos, np.inf, truePosPeptideScoresMean, truePosPeptideScoresStdev, len(identifiedPeptides)))
    peptideLabels = [[True, experimentId]]*len(identifiedPeptides)
    
    numTruePosPeptides = len(identifiedPeptides)
    numFalsePosPeptides = int(2 * numTruePosPeptides*peptideFdr / (1-peptideFdr))
    drawIndices = np.random.choice(len(peptideToProteinMap), numFalsePosPeptides)
    scores.extend(list(scipy.stats.truncnorm.rvs(minScoreFalsePos, np.inf, falsePosPeptideScoresMean, falsePosPeptideScoresStdev, numFalsePosPeptides)))
    
    falsePositivePeptides = [allPeptides[i] for i in drawIndices]
    identifiedPeptides.extend(falsePositivePeptides)
    peptideLabels.extend([[False, experimentId]]*numFalsePosPeptides)
    
    scoresCombined.extend(scores)
    identifiedPeptidesCombined.extend(identifiedPeptides)
    peptideLabelsCombined.extend(peptideLabels)
  
  presentProteinsCombined = set(presentProteinsCombined)
  
  if recalculatePeptFDR:
    seenPeptides = set()
    filteredPeptides = list()
    tp, fp = 0, 0
    for score, peptide, label in sorted(zip(scoresCombined, identifiedPeptidesCombined, peptideLabelsCombined), reverse = True):
      if peptide not in seenPeptides:
        seenPeptides.add(peptide)
        isTarget = not (sum([1 for x in peptideToProteinMap[peptide] if "REV__" in x]) == len(peptideToProteinMap[peptide]))
        if isTarget:
          tp += 1
        else:
          fp += 1
        if float(fp) / tp >= peptideFdr:
          print("recalculated peptide FDR score threshold:", score)
          break
        else:
          filteredPeptides.append((score, peptide, label))
    scoresCombined, identifiedPeptidesCombined, peptideLabelsCombined = zip(*filteredPeptides)  
  
  print("# present proteins:", len(presentProteinsCombined))
  print("# peptides:", len(identifiedPeptidesCombined))
  
  print("Calculating PEPs")
  scoresCombined = np.array(scoresCombined)
  PEPsCombined = scipy.stats.norm.pdf(scoresCombined, falsePosPeptideScoresMean, falsePosPeptideScoresStdev) * incorrectRatio / (scipy.stats.norm.pdf(scoresCombined, truePosPeptideScoresMean, truePosPeptideScoresStdev) * (1 - incorrectRatio) + scipy.stats.norm.pdf(scoresCombined, falsePosPeptideScoresMean, falsePosPeptideScoresStdev) * incorrectRatio)
  
  print("Protein grouping")  
  proteinGroups = replicateMaxQuantGrouping(identifiedPeptidesCombined, peptideToProteinMap)
  
  proteinToGroupMap = getProteinToGroupMap(proteinGroups)
  
  print("Calculating number of peptides per protein")
  proteinToObservedPeptidesMap = collections.defaultdict(list)
  for peptideIdx, (peptide, score, PEP) in enumerate(zip(identifiedPeptidesCombined, scoresCombined, PEPsCombined)):
    for protein in peptideToProteinMap[peptide]:
      proteinToObservedPeptidesMap[protein].append(peptide)
  
  proteinToObservedPeptidesMap = { protein : len(set(peptides)) for protein, peptides in proteinToObservedPeptidesMap.items() }
  
  entrapmentMarkedProteinGroups = list()
  for proteinGroup in proteinGroups:    
    proteinGroup = sorted(proteinGroup, key = lambda x : proteinToObservedPeptidesMap[x], reverse = True)
    proteinGroup = markEntrapment(proteinGroup, presentProteinsCombined)
    entrapmentMarkedProteinGroups.append(proteinGroup)
    #writerProteinGroups.writerow([";".join(proteinGroup), proteinScore, truePositive])
  
  print("Writing evidence.txt")
  writer = csv.writer(open(simulatedEvidenceFile, 'w'), delimiter = '\t')
  writer.writerow(['Modified sequence', 'Proteins', 'Leading proteins', 'Leading razor protein', 'Score', 'PEP', 'Experiment', 'True positive', 'Target'])
  proteinScoresList = collections.defaultdict(list)
  for peptideIdx, (peptide, score, PEP, (truePositive, experimentId)) in enumerate(zip(identifiedPeptidesCombined, scoresCombined, PEPsCombined, peptideLabelsCombined)):
    leadingProteins = list(set([entrapmentMarkedProteinGroups[proteinToGroupMap[protein]][0] for protein in peptideToProteinMap[peptide]]))
    leadingRazorProtein = sorted([(proteinToObservedPeptidesMap.get(p.replace("_entrapment", ""), 0), p) for p in leadingProteins], reverse = True)[0][1]
    isTarget = not (sum([1 for x in peptideToProteinMap[peptide] if "REV__" in x]) == len(peptideToProteinMap[peptide]))
    
    proteinScoresList[leadingRazorProtein].append((PEP, peptide))
    
    writer.writerow(["_" + peptide + "_", ";".join(peptideToProteinMap[peptide]), ";".join(leadingProteins), leadingRazorProtein, score, PEP, experimentId, truePositive, isTarget])
  
  proteinScores = calculateProteinScores(proteinScoresList)
  
  print("Writing proteinGroups.txt")
  writerProteinGroups = csv.writer(open(simulatedProteinGroups, 'w'), delimiter = '\t')
  writerProteinGroups.writerow(['Protein IDs', 'Score', 'True positive'])
  for proteinGroup in entrapmentMarkedProteinGroups:
    truePositive = False
    for protein in proteinGroup:
      if protein in presentProteinsCombined:
        truePositive = True
        break
    proteinScore = proteinScores.get(proteinGroup[0], -100)
    writerProteinGroups.writerow([";".join(proteinGroup), proteinScore, truePositive])

def calculateProteinScores(proteinScoresList):  
  proteinScores = dict()
  for leadingProtein in proteinScoresList:
    multPEP = 0.0
    seenPeptides = set()
    for PEP, peptide in proteinScoresList[leadingProtein]:
      if peptide not in seenPeptides:
        seenPeptides.add(peptide)
        multPEP -= np.log10(PEP / 1e-1)
    proteinScores[leadingProtein] = multPEP
  
  return proteinScores
    
def markEntrapment(proteins, presentProteins):
  return [p if p.replace("REV__", "") in presentProteins else p + "_entrapment" for p in proteins]

def getPeptideToProteinMap(peptideToProteinMapFile):
  peptideToProteinMap = dict()
  reader = csv.reader(open(peptideToProteinMapFile, 'r'), delimiter = '\t')
  for row in reader:
    peptideToProteinMap[row[0]] = row[1].split(";")
  return peptideToProteinMap

def getProteinToGeneMap(proteinToGeneMapFile):
  proteinToGeneMap = dict()
  geneToProteinMap = collections.defaultdict(list)
  reader = csv.reader(open(proteinToGeneMapFile, 'r'), delimiter = '\t')
  for row in reader:
    proteinToGeneMap[row[0]] = row[1]
    geneToProteinMap[row[1]].append(row[0])
  return proteinToGeneMap, geneToProteinMap
  
def getProbabilities(proteotypicityFile, peptideToProteinMap, proteinToGeneMap):
  # JPT input file: 
  # delimiter, proteinCol, peptideCol, psmsCol, occurenceCol, proteotypicityCol = ',', 2, 5, 6, 7, 8
  
  # ProteomicsDB API input file: 
  delimiter, proteinCol, peptideCol, psmsCol, occurenceCol, proteotypicityCol = '\t', 0, 1, 4, 3, 2
  
  reader = csv.reader(open(proteotypicityFile, 'r'), delimiter = delimiter)
  
  seenPeptides = set()
  psmsPerProtein = collections.defaultdict(int)
  psmsPerGene = collections.defaultdict(int)
  peptideProbsPerProtein = collections.defaultdict(list)
  for row in reader:
    protein = row[proteinCol]
    peptide = row[peptideCol]
    psms = int(row[psmsCol])
    occurrence = int(float(row[occurenceCol]))
    proteotypicity = float(row[proteotypicityCol])
    #print(protein, peptide, psms, occurrence, proteotypicity)
    if peptide in peptideToProteinMap:
      if len(proteinToGeneMap) > 0:
        if peptide not in seenPeptides:
          seenPeptides.add(peptide)
        psmsPerGene[proteinToGeneMap[protein]] += psms
      
      psmsPerProtein[protein] += psms
      peptideProbsPerProtein[protein].append([peptide, proteotypicity])
  
  proteinProbs = getProteinProbs(psmsPerProtein)
  geneProbs = getProteinProbs(psmsPerGene)
  
  return proteinProbs, peptideProbsPerProtein, geneProbs
  
def getProteinProbs(psmsPerProtein):  
  totalPsms = float(sum(psmsPerProtein.values()))
  proteinProbs = {key: value / totalPsms for (key, value) in psmsPerProtein.items()}
  #plt.hist(proteinProbs.values(), bins = np.linspace(0,0.0005,50))
  #plt.show()
  
  return proteinProbs

def replicateMaxQuantGrouping(observedPeptides, peptideToProteinMap):
  observedPeptideToProteinMap = dict()
  observedProteinToPeptideMap = collections.defaultdict(list)
  
  for peptide in observedPeptides:
    proteins = peptideToProteinMap[peptide]
    
    observedPeptideToProteinMap[peptide] = proteins
    for protein in proteins:
      observedProteinToPeptideMap[protein].append(peptide)
  
  proteinGroups = list()
  proteinToGroupMap = dict()
  proteinGroupIdx = 0
  
  for protein in observedProteinToPeptideMap:
    proteinGroups.append([protein])
    proteinToGroupMap[protein] = proteinGroupIdx
    proteinGroupIdx += 1
  
  subsetProteins = collections.defaultdict(list)
  for protein, peptides in observedProteinToPeptideMap.items():
    candidateProteins = observedPeptideToProteinMap[peptides[0]]
    for peptide in peptides[1:]:
      newCandidateProteins = list()
      for candidateProtein in candidateProteins:
        if candidateProtein in observedPeptideToProteinMap[peptide]:
          newCandidateProteins.append(candidateProtein)
      candidateProteins = newCandidateProteins
    
    moved = False
    # peptides from "protein" form a subset of each candidateProtein's peptides
    if len(candidateProteins) > 1:
      candidateProteins = sorted(candidateProteins, key = lambda x : len(observedProteinToPeptideMap[x]), reverse = True)
      for candidateProtein in candidateProteins:
        if candidateProtein != protein:
          if len(proteinGroups[proteinToGroupMap[candidateProtein]]) > 0:
            proteinGroups[proteinToGroupMap[candidateProtein]].extend(proteinGroups[proteinToGroupMap[protein]])
            moved = True
            break
    
    if moved:
      proteinGroups[proteinToGroupMap[protein]] = []
  
  proteinGroups = list(filter(lambda x : len(x) > 0, proteinGroups))
  
  #for pg in proteinGroups:
  #  if len(pg) > 1:
  #    print(";".join(pg))
  
  return proteinGroups

def getProteinToGroupMap(proteinGroups):
  proteinToGroupMap = dict()
  for proteinGroupIdx, proteinGroup in enumerate(proteinGroups):
    for protein in proteinGroup:
      proteinToGroupMap[protein] = proteinGroupIdx
  
  return proteinToGroupMap


def plotScoreDistributions(truePosPeptideScoresMean, truePosPeptideScoresStdev, falsePosPeptideScoresMean, falsePosPeptideScoresStdev, minScore, incorrectRatio):  
  import matplotlib.pyplot as plt
  
  xs = np.linspace(-2,5,100)
  plt.plot(xs, (1.0 - incorrectRatio) * scipy.stats.norm.pdf(xs, truePosPeptideScoresMean, truePosPeptideScoresStdev), label = 'Correct peptides')
  plt.plot(xs, incorrectRatio * scipy.stats.norm.pdf(xs, falsePosPeptideScoresMean, falsePosPeptideScoresStdev), label = 'Incorrect peptides')
  print((1.0 - incorrectRatio) * (1.0 - scipy.stats.norm.cdf(minScore, truePosPeptideScoresMean, truePosPeptideScoresStdev)))
  print(incorrectRatio * (1.0 - scipy.stats.norm.cdf(minScore, falsePosPeptideScoresMean, falsePosPeptideScoresStdev)))
  
  #PEPsCombined = scipy.stats.norm.pdf(xs, falsePosPeptideScoresMean, falsePosPeptideScoresStdev) * incorrectRatio / (scipy.stats.norm.pdf(xs, truePosPeptideScoresMean, truePosPeptideScoresStdev) * (1 - incorrectRatio) + scipy.stats.norm.pdf(xs, falsePosPeptideScoresMean, falsePosPeptideScoresStdev) * incorrectRatio)
  #plt.plot(xs, PEPsCombined, label = "PEP")
  
  plt.plot([minScore, minScore], [0, 0.5], 'k--', label = '1% FDR')
  plt.xlabel('Score')
  plt.ylabel('Probability')
  plt.legend()
  plt.show()
    
if __name__ == "__main__":
  main(sys.argv[1:])
