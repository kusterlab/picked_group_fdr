from timeit import Timer
import random
import pytest
import numpy as np

from picked_group_fdr.quant.lfq import LFQIntensityColumns
from picked_group_fdr.quant.precursor_quant import PrecursorQuant
from picked_group_fdr.results import ProteinGroupResult


def test_performanceLFQ():
    num_experiments = 600 # most experiments = 600 for prosit job 964 (Chengdong plasma samples)
    num_peptides = 600 # most peptides = 22500 for TTN in WP3 pipeline
    num_proteins = 50
    
    proteinGroupResults = []
    for i in range(num_proteins):
        peptideIntensityList = getPeptideIntensityList(num_experiments, num_peptides)
        
        random.shuffle(peptideIntensityList)
        
        pgr = ProteinGroupResult()
        pgr.precursorQuants = peptideIntensityList
        #pgr.precursorQuants = np.ones(100000000)
        
        proteinGroupResults.append(pgr)
    
    columns = [LFQIntensityColumns(silacChannels=[], minPeptideRatiosLFQ=1, stabilizeLargeRatiosLFQ=False, numThreads=1)]
    
    experimentToIdxMap = getExperimentToIdxMap(num_experiments)
    for c in columns:
        #c.append_headers(proteinGroupResults, experiments)
        c.append_columns(proteinGroupResults, experimentToIdxMap, postErrProbCutoff=1.0)

    #t = Timer(lambda: _getLFQIntensities(peptideIntensityList, experimentToIdxMap, 0.1, 1))
    #print(f"execution time for {num_proteins} runs: {t.timeit(number = num_proteins)}")


def getPeptideIntensityList(num_experiments, num_peptides):
    peptideIntensityList = list()
    for experiment in range(num_experiments):
        for peptide in range(num_peptides):
            x = random.uniform(1.0, 100.0)
            if x > 80.0: # test performance on sparse matrix
                continue
            peptideIntensityList.append(PrecursorQuant(f'_APEPTIDE{peptide}_', 2, f'file{experiment}', 1, x, 0.001, None, None, 1))

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
    print("Starting LFQ performance test")
    test_performanceLFQ()
