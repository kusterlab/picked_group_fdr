from dataclasses import dataclass
from typing import Optional

import numpy as np


@dataclass
class PrecursorQuant:
    peptide: str
    charge: int
    experiment: str
    fraction: int
    intensity: float
    post_err_prob: float
    tmt_intensities: Optional[np.array]
    silac_intensities: Optional[np.array]
    evidence_id: int
    assigned_mods: Optional[str] = None
    observed_mods: Optional[str] = None
