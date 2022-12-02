#ifndef PRECURSOR_QUANT_H
#define PRECURSOR_QUANT_H

#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <iostream>

class PrecursorQuantC {
  public:
    int peptide;
    int charge;
    int experiment;
    int fraction;
    double intensity;
    double postErrProb;
    std::vector<double> silacIntensities;
    int evidenceId;
    
    PrecursorQuantC() {};
    PrecursorQuantC(
        int peptide, 
        int charge, 
        int experiment, 
        int fraction, 
        double intensity, 
        double postErrProb, 
        int evidenceId) : 
        peptide(peptide), 
        charge(charge),
        experiment(experiment),
        fraction(fraction),
        intensity(intensity),
        postErrProb(postErrProb),
        evidenceId(evidenceId) {};
    
    friend bool operator< (const PrecursorQuantC& p1, const PrecursorQuantC& p2);
    friend bool operator== (const PrecursorQuantC& p1, const PrecursorQuantC& p2);
};

bool operator< (const PrecursorQuantC& p1, const PrecursorQuantC& p2) {
    return (p1.peptide < p2.peptide) ||
           (p1.peptide == p2.peptide && p1.charge < p2.charge) ||
           (p1.peptide == p2.peptide && p1.charge == p2.charge && p1.experiment < p2.experiment) ||
           (p1.peptide == p2.peptide && p1.charge == p2.charge && p1.experiment == p2.experiment && p1.fraction < p2.fraction) ||
           (p1.peptide == p2.peptide && p1.charge == p2.charge && p1.experiment == p2.experiment && p1.fraction == p2.fraction && p1.intensity > p2.intensity) ||
           (p1.peptide == p2.peptide && p1.charge == p2.charge && p1.experiment == p2.experiment && p1.fraction == p2.fraction && p1.intensity == p2.intensity && p1.postErrProb < p2.postErrProb);
}

bool operator== (const PrecursorQuantC& p1, const PrecursorQuantC& p2) {
    return p1.peptide == p2.peptide && p1.charge == p2.charge && p1.experiment == p2.experiment && p1.fraction == p2.fraction;
}

class IsMissingOrUnidentified {
    double postErrProbCutoff;

  public:
    IsMissingOrUnidentified(double postErrProbCutoff) :
        postErrProbCutoff(postErrProbCutoff)
    {}

    bool operator()(const PrecursorQuantC& p) const {
        return p.intensity <= 0.0 || (!isnan(p.postErrProb) && p.postErrProb > postErrProbCutoff);
    }
};


inline size_t key(int i,int j) {return (size_t) i << 32 | (unsigned int) j;}

double get_peptide_intensity_matrix(
        std::vector<PrecursorQuantC>& peptideIntensityList,
        int numSilacChannels,
        int numExperiments,
        std::vector<std::vector<double> >& peptideIntensities) {
    std::unordered_map<size_t, size_t> precursorMap;
    double totalIntensity = 0.0;
    for (auto& p : peptideIntensityList) {
        /* for each (peptide, charge) tuple, only sum the intensities of the PSMs
           with the highest intensity per fraction */
        size_t precursor = key(p.peptide, p.charge);
        if (precursorMap.find(precursor) == precursorMap.end()) {
            std::vector<double> v(numExperiments, 0.0);
            precursorMap[precursor] = peptideIntensities.size();
            peptideIntensities.push_back(v);
        }
        
        if (numSilacChannels > 0) {
            for (size_t silacIdx = 0; silacIdx < numSilacChannels; ++silacIdx) {
                double silacIntensity = p.silacIntensities.at(silacIdx);
                size_t silacExpIdx = p.experiment * numSilacChannels + silacIdx;
                peptideIntensities[precursorMap[precursor]][silacExpIdx] += (
                        silacIntensity);
                totalIntensity += silacIntensity;
            }
        } else {
            peptideIntensities[precursorMap[precursor]][p.experiment] += p.intensity;
            totalIntensity += p.intensity;
        }
    }
    return totalIntensity;
}

class FloatWithNanSort {
  public:
    FloatWithNanSort() {};

    bool operator()(double x, double y) const {
        if (isnan(x)) return false;
        if (isnan(y)) return true;
        return x < y;
    }
};

double nanmedian(std::vector<double>& v) {
    size_t m = v.size() / 2;
    std::nth_element(v.begin(), v.begin()+m, v.end(), FloatWithNanSort());
    if (v.size() % 2 == 0) {
        return (v[m] + v[m-1]) / 2;
    }
    return v[m];
}

void divide(double* v1, double* v2, size_t n, std::vector<double>& result) {
    for (size_t i = 0; i < n; ++i) {
        result[i] = v1[i] / v2[i];
    }
}

void divide_all(double* v1, size_t m, size_t n, std::vector<double>& results) {
    std::vector<double> tmp(n);
    for (size_t j = 0; j < m; ++j) {
        for (size_t k = j+1; k < m; ++k) {
            for (size_t i = 0; i < n; ++i) {
                tmp[i] = v1[i+m*j] / v1[i+m*k];
            }
            results.push_back(nanmedian(tmp));
        }
    }
}

#endif
