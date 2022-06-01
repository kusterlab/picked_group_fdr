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
