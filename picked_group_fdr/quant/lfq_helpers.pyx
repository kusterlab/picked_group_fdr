# cython: linetrace=True
# cython: binding=True
# distutils: language = c++
# distutils: define_macros=CYTHON_TRACE=1
cdef extern from "<algorithm>" namespace "std":
    Iter remove_if[Iter, UnaryPred](Iter first, Iter last, UnaryPred pred)
    void for_each[Iter, UnaryFunction](Iter first, Iter last, UnaryFunction f) except +  # actually returns f
    void sort[Iter](Iter first, Iter last) except +
    Iter unique[Iter](Iter first, Iter last) except +
    

# cythonize -a -i lfq_helpers.pyx
cimport cython
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.set cimport set as c_set
from libcpp cimport bool

import itertools
import sys
import collections
from typing import List, Dict, Tuple

import numpy as np
cimport numpy as np
import bottleneck as bn

from .precursor_quant import PrecursorQuant
from .. import helpers

cdef extern from "precursor_quant.h":
    cdef cppclass PrecursorQuantC:
        int peptide
        int charge
        int experiment
        int fraction
        double intensity
        double postErrProb
        vector[double] silacIntensities
        int evidenceId
        
        PrecursorQuantC()
        PrecursorQuantC(int, int, int, int, double, double, int)

    cdef cppclass IsMissingOrUnidentified:
        IsMissingOrUnidentified(double postErrProbCutoff)
    
    double get_peptide_intensity_matrix(
        vector[PrecursorQuantC]& peptideIntensityList,
        int numSilacChannels,
        int numExperiments,
        vector[vector[double]]& peptideIntensities)
    
    cdef cppclass FloatWithNanSort:
        FloatWithNanSort()

    double nanmedian(vector[double]& v)
    
    void divide(double* v1, 
                double* v2, 
                size_t n,
                vector[double]& result)
    void divide_all(double* v1, 
                size_t m, 
                size_t n,
                vector[double]& result)


def _getPeptideIntensities(
        peptideIntensityListPython: List[PrecursorQuant],
        experimentToIdxMap: Dict[str, int],
        postErrProbCutoff: float,
        numSilacChannels: int,
        numExperiments: int) \
        -> Tuple[Dict[Tuple[str, int], List[float]], float]:
    """
    Collects all precursor intensities per experiment
    """
    cdef vector[PrecursorQuantC] peptideIntensityList
    
    # according to line_profiler, this for loop is ~15% of the total runtime of 
    # _getLFQIntensities(). However, without the profiler, it is only <8%
    for p in peptideIntensityListPython:
        peptideIntensityList.push_back(PrecursorQuantC(hash(p.peptide), p.charge, experimentToIdxMap[p.experiment], p.fraction, p.intensity, p.postErrProb, p.evidenceId))
        if not p.silacIntensities is None:
            peptideIntensityList.back().silacIntensities = p.silacIntensities
    
    # filter missing and unidentified precursors
    peptideIntensityList.erase(
        remove_if(peptideIntensityList.begin(), 
                  peptideIntensityList.end(), 
                  IsMissingOrUnidentified(postErrProbCutoff)), 
        peptideIntensityList.end())
    
    # for each (peptide, charge, experiment, fraction) tuple, only keep the one 
    # with the highest intensity, in case of a tie, select the lowest PEP (= most confident PSM) 
    sort(peptideIntensityList.begin(), peptideIntensityList.end())
    peptideIntensityList.erase(unique(peptideIntensityList.begin(), peptideIntensityList.end()), peptideIntensityList.end())

    cdef vector[vector[double]] peptideIntensities
    totalIntensity = get_peptide_intensity_matrix(peptideIntensityList, numSilacChannels, numExperiments, peptideIntensities)
    intensityMatrix = np.array(peptideIntensities)
    
    return intensityMatrix, totalIntensity


def _getLogMedianPeptideRatiosLoop(columns, valid_vals, minPeptideRatiosLFQ, intensityMatrix):    
    if len(columns) == 0:
        return dict()
    
    cdef vector[double] ratios
    ratios = np.ones(len(columns[0][-1]))
    
    peptideRatios, experimentPairs = list(), list()
    for idx_i, i, col_i in columns:
        for idx_j, j, col_j in columns[idx_i+1:]:
            if valid_vals[idx_i, idx_j] < minPeptideRatiosLFQ:
                continue
            
            divide_wrapper(col_i, col_j, ratios)
            
            peptideRatios.append(nanmedian(ratios)) # this is the slow part: ~60% of the total runtime of _getLFQIntensities()
            experimentPairs.append((i,j))
    
    return dict(zip(experimentPairs, np.log(peptideRatios)))


cdef divide_wrapper(np.ndarray[np.double_t, ndim=1] col_i, np.ndarray[np.double_t, ndim=1] col_j, vector[double]& ratios):
    divide(&col_i[0], &col_j[0], len(col_i), ratios)


cdef divide_all_wrapper(np.ndarray[np.double_t, ndim=2] columns, vector[double]& ratios):
    """
    tried doing it all in C++ but it's not faster...
    cdef vector[double] peptideRatios
    divide_all_wrapper(intensityMatrix.T, peptideRatios)
    """
    divide_all(&columns[0, 0], columns.shape[0], columns.shape[1], ratios)
