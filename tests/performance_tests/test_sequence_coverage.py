from timeit import Timer
import random
from string import ascii_uppercase

import pytest
import numpy as np

from picked_group_fdr.columns.sequence_coverage import SequenceCoverageColumns
from picked_group_fdr.precursor_quant import PrecursorQuant


def test_performanceSequenceCoverage():
    num_experiments = 150
    num_peptides = 1250
    
    protein_sequence = getLongSequence()
    peptideIntensityList = getPeptideIntensityList(num_experiments, num_peptides, protein_sequence)
    pq = getSequenceCoverageColumns({'proteinA' : protein_sequence})
    experimentToIdxMap = getExperimentToIdxMap(num_experiments)
    
    random.shuffle(peptideIntensityList)
    
    t = Timer(lambda: pq.get_sequence_coverages(peptideIntensityList, experimentToIdxMap, 0.1, 'proteinA'))
    num_proteins = 1
    print(f"execution time for {num_proteins} runs: {t.timeit(number = num_proteins)}")


def getLongSequence():
    return ''.join(random.choice(ascii_uppercase) for i in range(38138)) # length is taken from TITIN_HUMAN


def getPeptideIntensityList(num_experiments, num_peptides, protein_sequence):
    peptideIntensityList = list()
    for experiment in range(num_experiments):
        for peptide in range(num_peptides):
            start = random.randint(0, len(protein_sequence) - 60)
            length = random.randint(6, 50)
            peptide_sequence = protein_sequence[start:start+length]
            peptideIntensityList.append(PrecursorQuant(peptide_sequence, 2, f'file{experiment}', 1, random.uniform(1.0, 100.0), 0.001, [], [], 1))

    return peptideIntensityList 


def getSequenceCoverageColumns(proteinSequences):
    return SequenceCoverageColumns(proteinSequences)
      

def getExperimentToIdxMap(num_experiments):
    experiments = [f'file{i}' for i in range(num_experiments)]
    experimentToIdxMap = dict([(v,k) for k, v in enumerate(experiments)])
    return experimentToIdxMap


if __name__ == "__main__":
    #import cProfile
    #cProfile.run("test_performanceSequenceCoverage()")
    test_performanceSequenceCoverage()
