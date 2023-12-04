from __future__ import annotations

import collections
import logging
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional

import numpy as np

from . import writers
from . import helpers
from .parsers import tsv

# for type hints only
from .peptide_info import ProteinGroupPeptideInfos
from .protein_groups import ProteinGroups
from . import columns

logger = logging.getLogger(__name__)


@dataclass
class ProteinGroupResult:
    proteinIds: str = ""
    majorityProteinIds: str = ""
    peptideCountsUnique: str = ""
    bestPeptide: str = ""
    numberOfProteins: int = 0
    qValue: float = np.nan
    score: float = np.nan
    reverse: str = ""
    potentialContaminant: str = ""
    precursorQuants: List[columns.PrecursorQuant] = field(default_factory=list)
    extraColumns: List[float] = field(default_factory=list)

    def extend(self, e):
        self.extraColumns.extend(e)

    def append(self, a):
        self.extraColumns.append(a)

    @classmethod
    def from_protein_group(
        cls,
        protein_group,
        peptide_scores,
        reported_fdr,
        protein_score,
        score_cutoff,
        keep_all_proteins,
    ):
        num_unique_peptides_per_protein = cls._get_peptide_counts(
            peptide_scores, score_cutoff
        )
        peptide_counts_unique = [
            num_unique_peptides_per_protein[p] for p in protein_group
        ]
        if sum(peptide_counts_unique) == 0 and not keep_all_proteins:
            return None

        protein_group, peptide_counts_unique = zip(
            *[
                (p, num_peptides)
                for p, num_peptides in zip(protein_group, peptide_counts_unique)
                if num_peptides > 0 or keep_all_proteins
            ]
        )

        best_peptide = sorted([(p[0], p[1]) for p in peptide_scores])[0][1]
        majority_protein_ids = ";".join(
            [
                p
                for p, num_peptides in zip(protein_group, peptide_counts_unique)
                if num_peptides >= max(peptide_counts_unique) / 2
            ]
        )
        number_of_proteins = len(protein_group)

        peptide_counts_unique = ";".join(map(str, peptide_counts_unique))
        protein_ids = ";".join(protein_group)

        qval = reported_fdr
        score = protein_score
        reverse = "+" if helpers.is_decoy(protein_group) else ""
        potential_contaminant = "+" if helpers.is_contaminant(protein_group) else ""
        return cls(
            protein_ids,
            majority_protein_ids,
            peptide_counts_unique,
            best_peptide,
            number_of_proteins,
            qval,
            score,
            reverse,
            potential_contaminant,
        )

    @staticmethod
    def _get_peptide_counts(score_peptide_pairs, score_cutoff):
        protein_peptide_count = collections.defaultdict(int)
        seen_peptides = set()
        for post_err_prob, peptide, proteins in sorted(score_peptide_pairs):
            if post_err_prob > score_cutoff:
                break
            if peptide not in seen_peptides:
                seen_peptides.add(peptide)
                for protein in proteins:
                    protein_peptide_count[protein] += 1
        return protein_peptide_count

    def to_list(self, format_extra_columns=None):
        if format_extra_columns is None:
            format_extra_columns = writers.format_extra_columns
        return [
            self.proteinIds,
            self.majorityProteinIds,
            self.peptideCountsUnique,
            self.bestPeptide,
            self.numberOfProteins,
            self.qValue,
            self.score,
            self.reverse,
            self.potentialContaminant,
        ] + [format_extra_columns(x) for x in self.extraColumns]


class ProteinGroupResults:
    headers: List[str]
    protein_group_results: List[ProteinGroupResult]
    experiments: List[str]
    num_tmt_channels: int = -1
    num_silac_channels: int = -1

    def __init__(self, protein_group_results: List[ProteinGroupResult] = None):
        """
        NOTE: cannot use empty list as default argument: https://docs.python-guide.org/writing/gotchas/
        """
        self.protein_group_results = []
        if protein_group_results is not None:
            self.protein_group_results = protein_group_results
        self.experiments = []
        self.headers = writers.PROTEIN_GROUP_HEADERS.copy()

    def __iter__(self):
        return iter(self.protein_group_results)

    def __next__(self) -> ProteinGroupResult:
        return next(self.protein_group_results)

    def __len__(self) -> int:
        return len(self.protein_group_results)

    def __getitem__(self, indices) -> ProteinGroupResult:
        return self.protein_group_results[indices]

    def append_header(self, header: str) -> None:
        if header in self.headers:
            raise ValueError(f"Trying to add duplicate column name: {header}")
        self.headers.append(header)

    def append_headers(self, headers: List[str]) -> None:
        for header in headers:
            self.append_header(header)

    def remove_column(self, header: str) -> None:
        if header not in self.headers:
            logger.warning(f"Attempted to remove non-existing column {header}")
            return

        column_idx = self.headers.index(header)
        del self.headers[column_idx]

        column_idx -= len(writers.PROTEIN_GROUP_HEADERS)
        for pgr in self.protein_group_results:
            del pgr.extraColumns[column_idx]

    def get_experiment_to_idx_map(self):
        return {experiment: idx for idx, experiment in enumerate(self.experiments)}

    def write(
        self,
        output_file: str,
        header_dict: Optional[Dict[str, str]] = None,
        format_extra_columns: Callable = None,
    ) -> None:
        writer = tsv.get_tsv_writer(output_file)
        if header_dict is None:
            writer.writerow(self.headers)
        else:
            writer.writerow(header_dict.keys())

        for protein_row in self.protein_group_results:
            out_row = protein_row.to_list(format_extra_columns)
            if header_dict is not None:
                out_row = [
                    out_row[self.headers.index(original_col_name)]
                    for original_col_name in header_dict.values()
                ]
            writer.writerow(out_row)

    def remove_protein_groups_without_precursors(self) -> None:
        filtered_protein_group_results = list()
        for pgr in self.protein_group_results:
            if len(pgr.precursorQuants) > 0:
                filtered_protein_group_results.append(pgr)

        self.protein_group_results = filtered_protein_group_results

    @classmethod
    def from_protein_groups(
        cls,
        protein_groups: ProteinGroups,
        protein_group_peptide_infos: ProteinGroupPeptideInfos,
        protein_scores: List[float],
        reported_qvals: List[float],
        score_cutoff: float,
        keep_all_proteins: bool,
    ) -> ProteinGroupResults:
        protein_group_results = list()
        for protein_group, peptide_scores, protein_score, reported_fdr in zip(
            protein_groups, protein_group_peptide_infos, protein_scores, reported_qvals
        ):
            if helpers.is_obsolete(protein_group):
                continue
            pgr = ProteinGroupResult.from_protein_group(
                protein_group,
                peptide_scores,
                reported_fdr,
                protein_score,
                score_cutoff,
                keep_all_proteins,
            )
            if (
                pgr is not None
            ):  # protein groups can get filtered out when they do not have any PSMs below the PSM FDR cutoff
                protein_group_results.append(pgr)

        return cls(protein_group_results)
