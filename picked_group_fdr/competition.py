from abc import ABC, abstractmethod
from typing import List, Tuple

import numpy as np

from .protein_groups import ProteinGroups
from .results import ProteinGroupResult
from .scoring import ProteinScoringStrategy, ProteinGroupPeptideInfos
from . import helpers


def ProteinCompetitionStrategyFactory(method='picked_group'):
    methods = {}
    methods['picked'] = PickedStrategy
    methods['picked_group'] = PickedGroupStrategy
    methods['classic'] = ClassicStrategy
    if method not in methods:
        raise ValueError(f"Unknown pickedStrategy {method['pickedStrategy']}, should be one of 'picked', 'picked_group' or 'classic'")
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
    
    def _select_proteins_for_picked(self, 
            proteins: List[str], 
            protein_group_score_list: List[Tuple[float,str,List[str]]]) -> List[str]:
        return proteins

    def do_competition(self, 
            proteinGroups: ProteinGroups, 
            proteinGroupPeptideInfos: ProteinGroupPeptideInfos,
            scoreType: ProteinScoringStrategy) -> Tuple[ProteinGroups, List, List[float]]:
        proteinScores = map(scoreType.calculate_score, proteinGroupPeptideInfos)
        isObsolete = map(helpers.isObsolete, proteinGroups)
        scoreGroupTuples = list(zip(proteinGroups, proteinGroupPeptideInfos, proteinScores, isObsolete))
        scoreGroupTuples = list(filter(lambda x : len(x[1]) > 0, scoreGroupTuples)) # remove protein groups without peptides
        
        # shuffle proteinGroups before sorting by protein score. This avoids 
        # biases between groups with equal score, e.g. when all target protein groups 
        # are listed above decoy protein groups.
        np.random.shuffle(scoreGroupTuples)
        scoreGroupTuples = sorted(scoreGroupTuples, key = lambda x : (x[2], not x[3]), reverse = True)
        
        filteredScoreGroupTuples = []
        for proteinGroup, proteinGroupScoreList, protein_score, _ in scoreGroupTuples:
            if self._is_protein_seen(proteinGroup) or helpers.isContaminant(proteinGroup):
                continue
            
            pickingProteins = self._select_proteins_for_picked(proteinGroup, proteinGroupScoreList)
            self._add_seen_proteins(pickingProteins)
            
            filteredScoreGroupTuples.append((proteinGroup, proteinGroupScoreList, protein_score))
        
        self.reset() # clears list of seen proteins
        
        # shuffle and sort again to randomly order obsolete and regular protein groups with the same score
        np.random.shuffle(filteredScoreGroupTuples)
        filteredScoreGroupTuples = sorted(filteredScoreGroupTuples, key = lambda x : x[2], reverse = True)
        
        filteredProteinGroups, filteredProteinGroupPeptideInfos, filteredProteinScores = zip(*filteredScoreGroupTuples)        
        return ProteinGroups(filteredProteinGroups), list(filteredProteinGroupPeptideInfos), filteredProteinScores


class PickedStrategy(ProteinCompetitionStrategy):
    seenProteins: set
    
    def __init__(self):
        self.seenProteins = set()
    
    def _add_seen_proteins(self, proteins: List[str]) -> None:
        self.seenProteins.add(self._get_protein_group_string(proteins))

    def _is_protein_seen(self, proteins: List[str]) -> bool:
        return self._get_protein_group_string(proteins) in self.seenProteins
    
    def _get_protein_group_string(self, proteins: List[str]) -> str:
        return ";".join(map(_clean_protein_id, proteins))

    def short_description(self) -> str:
        return "pT"
    
    def long_description(self) -> str:
        return "picked target-decoy strategy"
    
    def reset(self) -> None:
        self.seenProteins = set()


class PickedGroupStrategy(ProteinCompetitionStrategy):
    seenProteins: set
    
    def __init__(self, picking_strategy: str="leading"):
        self.seenProteins = set()
        self.picking_strategy = picking_strategy
    
    def _add_seen_proteins(self, proteins: List[str]) -> None:
        for protein in proteins:
            self.seenProteins.add(_clean_protein_id(protein))

    def _is_protein_seen(self, proteins: List[str]) -> bool:
        return any(_clean_protein_id(p) in self.seenProteins for p in proteins)
    
    def _select_proteins_for_picked(self, 
            proteins: List[str], 
            protein_group_score_list: List[Tuple[float,str,List[str]]]) -> List[str]:
        numUniquePeptidesPerProtein = ProteinGroupResult._get_peptide_counts(protein_group_score_list, 1.01)
        if self.picking_strategy == "all":
            return proteins
        elif self.picking_strategy == "majority":
            return [p for p in proteins if numUniquePeptidesPerProtein[p] >= max(numUniquePeptidesPerProtein.values()) / 2]
        elif self.picking_strategy == "leading":
            return [p for p in proteins if numUniquePeptidesPerProtein[p] == max(numUniquePeptidesPerProtein.values())]
        else:
            raise ValueError(f"Unknown picking strategy {self.picking_strategy}")
        
    def short_description(self) -> str:
        return "pgT"
    
    def long_description(self) -> str:
        return "picked group target-decoy strategy"
    
    def reset(self) -> None:
        self.seenProteins = set()


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
    return protein_id.replace("REV__", "").replace("OBSOLETE__", "")
