import sys
import csv
import os
import collections

def main(argv):
    proteinGroupsFile = argv[0]
    fastaFile = argv[1]
    proteinGroupsFileOut = argv[2]
    
    proteinToGeneMap = getProteinToGeneMap(fastaFile)
    
    reader = csv.reader(open(proteinGroupsFile, 'r'), delimiter = '\t')
    
    headers = next(reader)
    
    proteinIdsCol = headers.index('Protein IDs')
    qvalCol = headers.index('Q-value')
    numUniquePeptidesCol = headers.index('Peptide counts (unique)')
    bestPeptideCol = headers.index('Best peptide')
    scoreCol = headers.index('Score')
    reverseCol = headers.index('Reverse')
    potentialContaminantCol = headers.index('Potential contaminant')
    
    proteinsWithoutGeneName = 0
    proteinGroupsWithGenes = list()
    geneIdsCount = collections.defaultdict(int)
    for row in reader:
        if row[reverseCol] == "+" or row[potentialContaminantCol] == "+":
            continue

        proteinIds = row[proteinIdsCol].split(';')
        geneIds = makeListUnique(filter(lambda x : x != "None", map(lambda x : getGeneId(x, proteinToGeneMap), proteinIds)))
        if len(geneIds) > 0:
            leadingGene = geneIds[0]
            geneIds = ";".join(geneIds)
            proteinIds = row[proteinIdsCol]
            qval = float(row[qvalCol])
            numUniquePeptides = row[numUniquePeptidesCol]
            bestPeptide = row[bestPeptideCol]
            bestPeptidePEP = 10**(-1*float(row[scoreCol]))
            
            if qval < 0.01:
                proteinGroupsWithGenes.append([leadingGene, geneIds, proteinIds, qval, numUniquePeptides, bestPeptide, bestPeptidePEP])
                geneIdsCount[geneIds] += 1
        else:
            proteinsWithoutGeneName += 1

    print("#Protein groups without gene name:", proteinsWithoutGeneName)
    
    writer = csv.writer(open(proteinGroupsFileOut, 'w'), delimiter = '\t')
    writer.writerow(["Leading gene", "Gene group",  "Proteing Group", "Protein group q-value", "Number of unique peptides", "Best scoring peptide", "Best scoring peptide's PEP"])
    
    geneIdsWithMultipleGroups = [geneId for geneId, count in geneIdsCount.items() if count > 1]
    print("#Genes with multiple protein groups:", len(geneIdsWithMultipleGroups))
    
    for proteinGroup in sorted(proteinGroupsWithGenes):
        if proteinGroup[0] in geneIdsWithMultipleGroups:
            writer.writerow(proteinGroup)


def getGeneId(proteinId, proteinToGeneMap):
    if proteinId.replace("REV__", "") not in proteinToGeneMap:
        return "None"
        
    if proteinId.startswith("REV__"):
        return "REV__" + proteinToGeneMap[proteinId.replace("REV__", "")]
    else:
        return proteinToGeneMap[proteinId]


def makeListUnique(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


def getProteinToGeneMap(fastaFile):
    proteinToGeneMap = dict()
    with open(fastaFile, 'r') as f:
        for line in f:
            if line.startswith('>'):
                line = line[:-1]
                proteinId = line.split(" ")[0][1:]
                if "GN=" in line:
                    geneId = line.split("GN=")[1].split(" ")[0]
                else:
                    continue
                    #geneId = proteinId
                proteinToGeneMap[proteinId] = geneId
    return proteinToGeneMap


if __name__ == "__main__":
    main(sys.argv[1:])
