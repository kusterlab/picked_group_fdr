from dataclasses import dataclass, field
from typing import List

import numpy as np


@dataclass
class PrecursorQuant:
    peptide: str
    charge: int
    experiment: str
    fraction: int
    intensity: float
    postErrProb: float
    tmtIntensities: np.array
    silacIntensities: np.array
    evidenceId: int
    
    def _totuple(self):
        return (hash(self.peptide), self.charge, hash(self.experiment), self.fraction, self.intensity, self.postErrProb, self.evidenceId)
    
    def __array__(self):
        return np.array(self._totuple())
    
    def __len__(self):
        return self._totuple().__len__()
    
    def __getitem__(self, item):
        return self._totuple().__getitem__(item)        
