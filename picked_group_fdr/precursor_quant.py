from dataclasses import dataclass
from typing import Optional

import numpy as np


@dataclass
class PrecursorQuant:
    """Represents a quantified precursor.

    The peptide sequence can contain modifications in parentheses or square
    brackets but has no flanking characters, e.g. AP(ox)EPTIDE.

    assigned_mods and observed_mods fields are only used for FragPipe input.
    """    
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
