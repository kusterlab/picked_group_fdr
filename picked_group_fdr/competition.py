from __future__ import annotations

from abc import ABC, abstractmethod
from typing import List, Tuple

import numpy as np

from .protein_groups import ProteinGroups
from .results import ProteinGroupResult
from . import helpers

# for type hints only
from . import scoring_strategy
from . import peptide_info


def ProteinCompetitionStrategyFactory(
    method="picked_group",
) -> ProteinCompetitionStrategy:
    methods = {}
    methods["picked"] = PickedStrategy
    methods["picked_group"] = PickedGroupStrategy
    methods["classic"] = ClassicStrategy
    if method not in methods:
        raise ValueError(
            f"Unknown pickedStrategy {method['pickedStrategy']}, should be one of 'picked', 'picked_group' or 'classic'"
        )
    return methods[method]()


class ProteinCompetitionStrategy(ABC):
    @abstractmethod
    def _add_seen_proteins(self, proteins: List[str]) -> None:
        pass

    @abstractmethod
    def _is_protein_seen(self, proteins: List[str]) -> bool:
        pass

    @abstractmethod
    def short_description(self) -> str:
        pass

    @abstractmethod
    def long_description(self) -> str:
        pass

    @abstractmethod
    def reset(self) -> None:
        pass

    def _select_proteins_for_picked(
        self,
        proteins: List[str],
        protein_group_score_list: List[Tuple[float, str, List[str]]],
    ) -> List[str]:
        return proteins

    def do_competition(
        self,
        protein_groups: ProteinGroups,
        protein_group_peptide_infos: peptide_info.ProteinGroupPeptideInfos,
        score_type: scoring_strategy.ProteinScoringStrategy,
    ) -> Tuple[ProteinGroups, List, List[float]]:
        protein_scores = map(score_type.calculate_score, protein_group_peptide_infos)
        is_obsolete = map(helpers.is_obsolete, protein_groups)
        score_group_tuples = list(
            zip(
                protein_groups, protein_group_peptide_infos, protein_scores, is_obsolete
            )
        )
        score_group_tuples = list(
            filter(lambda x: len(x[1]) > 0, score_group_tuples)
        )  # remove protein groups without peptides

        # shuffle protein_groups before sorting by protein score. This avoids
        # biases between groups with equal score, e.g. when all target protein groups
        # are listed above decoy protein groups.
        np.random.shuffle(score_group_tuples)
        score_group_tuples = sorted(
            score_group_tuples, key=lambda x: (x[2], not x[3]), reverse=True
        )

        filtered_score_group_tuples = []
        for (
            protein_group,
            protein_group_score_list,
            protein_score,
            _,
        ) in score_group_tuples:
            if self._is_protein_seen(protein_group) or helpers.is_contaminant(
                protein_group
            ):
                continue

            picking_proteins = self._select_proteins_for_picked(
                protein_group, protein_group_score_list
            )
            self._add_seen_proteins(picking_proteins)

            filtered_score_group_tuples.append(
                (protein_group, protein_group_score_list, protein_score)
            )

        self.reset()  # clears list of seen proteins

        # shuffle and sort again to randomly order obsolete and regular protein groups with the same score
        np.random.shuffle(filtered_score_group_tuples)
        filtered_score_group_tuples = sorted(
            filtered_score_group_tuples, key=lambda x: x[2], reverse=True
        )

        (
            filtered_protein_groups,
            filtered_protein_group_peptide_infos,
            filtered_protein_scores,
        ) = zip(*filtered_score_group_tuples)
        return (
            ProteinGroups(filtered_protein_groups),
            list(filtered_protein_group_peptide_infos),
            filtered_protein_scores,
        )


class PickedStrategy(ProteinCompetitionStrategy):
    seen_proteins: set

    def __init__(self):
        self.seen_proteins = set()

    def _add_seen_proteins(self, proteins: List[str]) -> None:
        self.seen_proteins.add(self._get_protein_group_string(proteins))

    def _is_protein_seen(self, proteins: List[str]) -> bool:
        return self._get_protein_group_string(proteins) in self.seen_proteins

    def _get_protein_group_string(self, proteins: List[str]) -> str:
        return ";".join(map(_clean_protein_id, proteins))

    def short_description(self) -> str:
        return "pT"

    def long_description(self) -> str:
        return "picked target-decoy strategy"

    def reset(self) -> None:
        self.seen_proteins = set()


class PickedGroupStrategy(ProteinCompetitionStrategy):
    seen_proteins: set

    def __init__(self, picking_strategy: str = "leading"):
        self.seen_proteins = set()
        self.picking_strategy = picking_strategy

    def _add_seen_proteins(self, proteins: List[str]) -> None:
        for protein in proteins:
            self.seen_proteins.add(_clean_protein_id(protein))

    def _is_protein_seen(self, proteins: List[str]) -> bool:
        return any(_clean_protein_id(p) in self.seen_proteins for p in proteins)

    def _select_proteins_for_picked(
        self,
        proteins: List[str],
        protein_group_score_list: List[Tuple[float, str, List[str]]],
    ) -> List[str]:
        num_unique_peptides_per_protein = ProteinGroupResult._get_peptide_counts(
            protein_group_score_list, 1.01
        )
        if self.picking_strategy == "all":
            return proteins
        elif self.picking_strategy == "majority":
            return [
                p
                for p in proteins
                if num_unique_peptides_per_protein[p]
                >= max(num_unique_peptides_per_protein.values()) / 2
            ]
        elif self.picking_strategy == "leading":
            return [
                p
                for p in proteins
                if num_unique_peptides_per_protein[p]
                == max(num_unique_peptides_per_protein.values())
            ]
        else:
            raise ValueError(f"Unknown picking strategy {self.picking_strategy}")

    def short_description(self) -> str:
        return "pgT"

    def long_description(self) -> str:
        return "picked group target-decoy strategy"

    def reset(self) -> None:
        self.seen_proteins = set()


class ClassicStrategy(ProteinCompetitionStrategy):
    def _add_seen_proteins(self, proteins: List[str]) -> None:
        pass

    def _is_protein_seen(self, proteins: List[str]) -> bool:
        return False

    def short_description(self) -> str:
        return "cT"

    def long_description(self) -> str:
        return "classic target-decoy strategy"

    def reset(self) -> None:
        pass


def _clean_protein_id(protein_id):
    return protein_id.replace("REV__", "").replace("OBSOLETE__", "").replace("rev_", "")
