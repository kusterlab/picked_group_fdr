from timeit import Timer
import random
import pytest
import numpy as np

from picked_group_fdr.quant.lfq import LFQIntensityColumns
from picked_group_fdr.quant.precursor_quant import PrecursorQuant


def test_performanceLFQ():
    num_experiments = 600 # most experiments = 600 for prosit job 964 (Chengdong plasma samples)
    num_peptides = 1250 # most peptides = 1243 for Q8WZ42 in prosit job 776
    
    peptideIntensityList = getPeptideIntensityList(num_experiments, num_peptides)
    pq = getLFQIntensityColumns()
    experimentToIdxMap = getExperimentToIdxMap(num_experiments)
    
    random.shuffle(peptideIntensityList)
    
    t = Timer(lambda: pq._getLFQIntensities(peptideIntensityList, experimentToIdxMap, 0.1))
    num_proteins = 1
    print(f"execution time for {num_proteins} runs: {t.timeit(number = num_proteins)}")


def getPeptideIntensityList(num_experiments, num_peptides):
    peptideIntensityList = list()
    for experiment in range(num_experiments):
        for peptide in range(num_peptides):
            x = random.uniform(1.0, 100.0)
            if x > 80.0: # test performance on sparse matrix
                continue
            peptideIntensityList.append(PrecursorQuant(f'_APEPTIDE{peptide}_', 2, f'file{experiment}', 1, x, 0.001, [], [], 1))

    return peptideIntensityList


def getLFQIntensityColumns():
    return LFQIntensityColumns(silacChannels=[], minPeptideRatiosLFQ=2, stabilizeLargeRatiosLFQ=False)
      

def getExperimentToIdxMap(num_experiments):
    experiments = [f'file{i}' for i in range(num_experiments)]
    experimentToIdxMap = dict([(v,k) for k, v in enumerate(experiments)])
    return experimentToIdxMap


if __name__ == "__main__":
    #import cProfile
    #cProfile.run("test_performanceLFQ()")
    test_performanceLFQ()
