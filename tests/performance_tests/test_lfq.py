import random
import numpy as np
from timeit import default_timer as timer

import picked_group_fdr.columns as columns
from picked_group_fdr.precursor_quant import PrecursorQuant
from picked_group_fdr.results import ProteinGroupResults, ProteinGroupResult


def test_performanceLFQ():
    num_experiments = 100 # most experiments = 600 for prosit job 964 (Chengdong plasma samples)
    num_peptides = 600 # most peptides = 22500 for TTN in WP3 pipeline
    num_proteins = 200
    
    proteinGroupResults = []
    for i in range(num_proteins):
        peptideIntensityList = getPeptideIntensityList(num_experiments, num_peptides)
        
        random.shuffle(peptideIntensityList)
        
        pgr = ProteinGroupResult()
        pgr.precursorQuants = peptideIntensityList
        #pgr.precursorQuants = np.ones(100000000)
        
        proteinGroupResults.append(pgr)
    
    lfq_columns = [columns.LFQIntensityColumns(minPeptideRatiosLFQ=1, stabilizeLargeRatiosLFQ=False, numThreads=2)]
    
    start = timer()

    proteinGroupResults = ProteinGroupResults(proteinGroupResults)
    for experiment in range(num_experiments):
        proteinGroupResults.experiments.append(f'file{experiment}')
    
    for c in lfq_columns:
        c.append(proteinGroupResults, post_err_prob_cutoff=1.0)
    
    end = timer()
    print(f"execution took {'%.1f' % (end - start)} seconds wall clock time")


def getPeptideIntensityList(num_experiments, num_peptides):
    peptideIntensityList = list()
    for experiment in range(num_experiments):
        for peptide in range(num_peptides):
            x = random.uniform(1.0, 100.0)
            if x > 80.0: # test performance on sparse matrix
                continue
            peptideIntensityList.append(PrecursorQuant(f'APEPTIDE{peptide}', 2, f'file{experiment}', 1, x, 0.001, None, None, 1))

    return peptideIntensityList


def getLFQIntensityColumns():
    return columns.LFQIntensityColumns(minPeptideRatiosLFQ=2, stabilizeLargeRatiosLFQ=False)
      

def getExperimentToIdxMap(num_experiments):
    experiments = [f'file{i}' for i in range(num_experiments)]
    experimentToIdxMap = dict([(v,k) for k, v in enumerate(experiments)])
    return experimentToIdxMap



# Profile with cProfile:
if __name__ == "__main__":
    # Perform cProfile profiling (no overhead but harder to find exact bottleneck):
    # 1. run "python -m cProfile -o program.prof tests/performance_tests/test_lfq.py"
    # 2. run "snakeviz -H 0.0.0.0 -s program.prof"
    
    # Perform line-by-line profiling (significant overhead but easy to interpret):
    # 1. add @profile decorators in lfq.py
    # 2. run "kernprof -lv tests/performance_tests/test_lfq.py"

    # Perform memory profiling:
    # 1. mprof run --include-children --backend psutil_pss --python python tests/performance_tests/test_lfq.py
    # 2. mprof plot -o mprof_plot.png --backend agg
    
    print("Starting LFQ performance test")
    test_performanceLFQ()
