from abc import ABC, abstractmethod
from typing import Dict, List

from ..results import ProteinGroupResults


class ProteinGroupColumns(ABC):
    @abstractmethod
    def append_headers(
        self,
        proteinGroupResults: ProteinGroupResults,
        experiments: List[str],
    ) -> None:
        pass

    @abstractmethod
    def append_columns(
        self,
        proteinGroupResults: ProteinGroupResults,
        experimentToIdxMap: Dict[str, int],
        postErrProbCutoff: float,
    ) -> None:
        pass
