from __future__ import annotations
from typing import List, Set, Dict
import logging

from . import helpers
from . import entrapment
from .parsers import protein_groups as pgp

# for type hints only
from .peptide_info import ProteinGroupPeptideInfos


logger = logging.getLogger(__name__)


class ProteinGroups:
    protein_groups: List[List[str]]
    protein_to_group_idx_map: Dict[str, int]
    valid_idx: bool

    def __init__(self, protein_groups: List[List[str]] = None):
        if protein_groups is None:  # https://docs.python-guide.org/writing/gotchas/
            self.protein_groups = []
        else:
            self.protein_groups = protein_groups
        self.protein_to_group_idx_map = dict()
        self.valid_idx = False

    @staticmethod
    def from_mq_protein_groups_file(mq_protein_groups_file: str):
        protein_groups = [
            protein_group
            for protein_group, _ in pgp.parse_protein_groups_file_single(
                mq_protein_groups_file
            )
        ]
        return ProteinGroups.init_from_list(protein_groups)

    @staticmethod
    def from_observed_peptide_map(protein_to_peptides_dict):
        protein_groups = [[protein] for protein in protein_to_peptides_dict.keys()]
        return ProteinGroups.init_from_list(protein_groups)

    @staticmethod
    def from_protein_group_results(protein_groups_results):
        protein_groups = [pgr.proteinIds.split(";") for pgr in protein_groups_results]
        return ProteinGroups.init_from_list(protein_groups)

    @classmethod
    def init_from_list(cls, protein_groups):
        protein_groups = cls(protein_groups)
        protein_groups.create_index()
        return protein_groups

    def __iter__(self):
        return iter(self.protein_groups)

    def __next__(self):
        return next(self.protein_groups)

    def __len__(self):
        return len(self.protein_groups)

    def __eq__(self, other):
        return sorted(self.protein_groups) == sorted(other.protein_groups)

    def append(self, protein_group: List[str]):
        self.protein_groups.append(protein_group)
        self.valid_idx = False

    def extend(self, protein_groups: ProteinGroups):
        self.protein_groups.extend(protein_groups)
        self.valid_idx = False

    def size(self) -> int:
        return len(self.protein_groups)

    def get_all_proteins(self):
        return set(
            [
                protein
                for protein_group in self.protein_groups
                for protein in protein_group
            ]
        )

    def get_leading_proteins(self, proteins: List[str]) -> Set[str]:
        """Returns the leading protein for each of the protein groups that the proteins belong to."""
        get_leading_protein = lambda protein: self.get_protein_group(protein)[0]
        return set(map(get_leading_protein, proteins))

    def create_index(self):
        """Creates a map of proteins to their index in self.protein_groups."""
        self.protein_to_group_idx_map = dict()
        for protein_group_idx, protein_group in enumerate(self.protein_groups):
            for protein in protein_group:
                self.protein_to_group_idx_map[protein] = protein_group_idx
        self.valid_idx = True

    def get_protein_group(self, protein: str, check_idx_valid: bool = True):
        return self.protein_groups[
            self._get_protein_group_idx(protein, check_idx_valid)
        ]

    def _get_protein_group_idx(self, protein: str, check_idx_valid: bool = True):
        if not self.valid_idx and check_idx_valid:
            raise Exception("Trying to get group index while index is invalid")
        return self.protein_to_group_idx_map[protein]

    def get_protein_group_idxs(
        self, proteins: List[str], check_idx_valid: bool = True
    ) -> Set[int]:
        if not self.valid_idx and check_idx_valid:
            raise Exception("Trying to get group index while index is invalid")

        protein_group_idxs = set()
        for protein in proteins:
            if protein in self.protein_to_group_idx_map:
                protein_group_idxs.add(
                    self._get_protein_group_idx(protein, check_idx_valid)
                )
            else:
                protein_group_idxs.add(-1)
        return protein_group_idxs

    def get_protein_groups(
        self, proteins: List[str], check_idx_valid: bool = True
    ) -> List[List[str]]:
        protein_group_idxs = self.get_protein_group_idxs(proteins, check_idx_valid)
        return [self.protein_groups[idx] for idx in protein_group_idxs]

    def merge_groups(self, superset_protein: str, protein: str):
        """Merge two protein groups"""
        superset_group_idx = self._get_protein_group_idx(
            superset_protein, check_idx_valid=False
        )
        self.protein_groups[superset_group_idx].extend(
            self.get_protein_group(protein, check_idx_valid=False)
        )
        self.protein_groups[
            self._get_protein_group_idx(protein, check_idx_valid=False)
        ] = []
        self.valid_idx = False

    def add_unseen_protein_groups(
        self,
        protein_groups_other: ProteinGroups,
        protein_group_peptide_infos_other: ProteinGroupPeptideInfos,
    ):
        """Adds missing protein groups from another ProteinGroups object.

        Adds protein groups to the current ProteinGroups object from
        protein_groups_other that have proteins which are not represented in a
        protein group of the current ProteinGroups object.

        Additionally, this function returns obsolete protein groups and
        corresponding peptides needed for picked protein group competition.
        Obsolete protein groups are protein groups from the original protein
        grouping, i.e. before the rescuing step, that were merged in the
        rescuing step due to removal of low-confident unique peptides. These are
        needed for picked protein group competition using only leading proteins
        (the default setting) and deal with the following situation:

        protein groups before rescuing step:
        - proteinA, #peptides: 10, bestPEP: 1e-5
        - proteinB, #peptides: 3, bestPEP: 0.08
        - REV__proteinA, #peptides: 2, bestPEP: 0.13
        - REV__proteinB, #peptides: 2, bestPEP: 0.11

        protein groups after rescuing step:
        - proteinA;proteinB, #peptides: 5;1, bestPEP: 1e-5
        - REV__proteinA, #peptides: 1, bestPEP: 0.13
        - REV__proteinB, #peptides: 1, bestPEP: 0.11
        - OBSOLETE__proteinA, #peptides: 10, bestPEP: 1e-5
        - OBSOLETE__proteinB, #peptides: 3, bestPEP: 0.08

        For the picked protein group competition, only the leading protein
        (proteinA) is considered for the protein group proteinA;proteinB. The
        regular picked target-decoy competition is replaced by a
        target-decoy-obsolete competition. Without the obsolete protein group
        OBSOLETE__proteinB, rev_proteinB would never receive competition and
        will therefore always be retained. This leads to an overconservative FDR
        estimate.

        Obsolete protein groups are used for FDR calculation but not reported
        in the proteinGroups.txt (see ProteinGroupResults.from_protein_groups).
        """
        identified_proteins = self.get_all_proteins()
        obsolete_protein_groups = ProteinGroups()
        obsolete_protein_group_peptide_infos = []
        for protein_group_other, protein_group_peptide_info_other in zip(
            protein_groups_other, protein_group_peptide_infos_other
        ):
            unseen_protein_group = list(
                filter(
                    lambda protein: protein not in identified_proteins,
                    protein_group_other,
                )
            )
            if len(unseen_protein_group) == 0:
                # all proteins in protein_group_other are already represented
                # in a protein group in the current protein groups object.
                obsolete_protein_group = list(
                    map(lambda x: f"OBSOLETE__{x}", protein_group_other)
                )
                obsolete_protein_groups.append(obsolete_protein_group)
                obsolete_protein_group_peptide_infos.append(
                    protein_group_peptide_info_other
                )
            else:
                self.append(unseen_protein_group)

        self.create_index()

        return obsolete_protein_groups, obsolete_protein_group_peptide_infos

    def remove_empty_groups(self):
        self.protein_groups = list(filter(lambda x: len(x) > 0, self.protein_groups))
        self.create_index()

    def print(self, min_size: int = 0):
        for pg in self.protein_groups:
            if len(pg) >= min_size:
                logger.info(
                    f'{";".join(pg)} {"entrapmentGroup" if entrapment.is_entrapment(pg) else ""} {"decoyGroup" if helpers.is_decoy(pg) else ""}'
                )
