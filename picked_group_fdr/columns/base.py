from __future__ import annotations

from abc import ABC, abstractmethod

# for type hints only
from .. import results


class ProteinGroupColumns(ABC):
    def append(
        self,
        protein_group_results: results.ProteinGroupResults,
        post_err_prob_cutoff: float,
    ):
        if not self.is_valid(protein_group_results):
            return
        self.append_headers(protein_group_results)
        self.append_columns(protein_group_results, post_err_prob_cutoff)

    def is_valid(self, protein_group_results: results.ProteinGroupResults) -> bool:
        return True

    @abstractmethod
    def append_headers(
        self,
        proteinGroupResults: results.ProteinGroupResults,
    ) -> None:
        pass

    @abstractmethod
    def append_columns(
        self,
        proteinGroupResults: results.ProteinGroupResults,
        postErrProbCutoff: float,
    ) -> None:
        pass
