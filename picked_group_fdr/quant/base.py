from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Dict, List

from .. import results


class ProteinGroupColumns(ABC):
    @abstractmethod
    def append_headers(
        self,
        proteinGroupResults: results.ProteinGroupResults,
        experiments: List[str],
    ) -> None:
        pass

    @abstractmethod
    def append_columns(
        self,
        proteinGroupResults: results.ProteinGroupResults,
        experimentToIdxMap: Dict[str, int],
        postErrProbCutoff: float,
    ) -> None:
        pass
