from __future__ import annotations

import hashlib
import logging
from typing import Callable, Dict, List

from . import fdr, helpers
from .observed_peptides import ObservedPeptides
from .scoring import (
    BestAndromedaScore,
    BestPEPScore,
    MQProteinScore,
    MultPEPScore,
    ProteinScore,
)
from .score_origin import (
    FragPipeInput,
    MaxQuantInput,
    MaxQuantInputNoRemap,
    PercolatorInput,
    PercolatorInputRemapped,
    SageInput,
    ScoreOrigin,
)

# for type hints only
from .protein_groups import ProteinGroups
from .peptide_info import PeptideInfoList, ProteinGroupPeptideInfos

logger = logging.getLogger(__name__)


class ProteinScoringStrategy:
    use_razor: bool
    use_shared_peptides: bool
    protein_score: ProteinScore
    peptide_score_cutoff: float
    score_origin: ScoreOrigin
    peptide_counts_per_protein: Dict[str, int]
    best_peptide_score_per_protein: Dict[str, float]

    def __init__(
        self, score_description, mq_protein_groups_file=""
    ):
        if "multPEP" in score_description:
            self.protein_score = MultPEPScore()
        elif "bestPEP" in score_description:
            self.protein_score = BestPEPScore()
        elif "Andromeda" in score_description:
            self.protein_score = BestAndromedaScore()
        elif "MQ_protein" in score_description:
            self.protein_score = MQProteinScore(mq_protein_groups_file)
        else:
            raise NotImplementedError

        self.use_razor = "razor" in score_description
        self.use_shared_peptides = "with_shared" in score_description

        if "Perc" in score_description:
            # Add "remap" to the score_type if the fasta database used for protein grouping is different from the one used by Percolator
            if "remap" in score_description:
                self.score_origin = PercolatorInputRemapped()
            else:
                self.score_origin = PercolatorInput()
        elif "FragPipe" in score_description:
            self.score_origin = FragPipeInput()
        elif "Sage" in score_description:
            self.score_origin = SageInput()
        else:
            if "no_remap" in score_description:
                self.score_origin = MaxQuantInputNoRemap()
            else:
                self.score_origin = MaxQuantInput()

    def get_evidence_file(self, args) -> str:
        return self.score_origin.get_evidence_file(args)

    def get_evidence_parser(self) -> Callable:
        return self.score_origin.get_evidence_parser()

    def get_quantification_file(self, args) -> str:
        return self.score_origin.get_quantification_file(args)

    def get_quantification_parser(self) -> Callable:
        return self.score_origin.get_quantification_parser()

    def remaps_peptides_to_proteins(self) -> bool:
        return self.score_origin.remaps_peptides_to_proteins()

    def can_do_quantification(self) -> bool:
        return self.score_origin.can_do_quantification()

    def optimize_hyperparameters(
        self, protein_groups, protein_group_peptide_infos, protein_group_fdr_threshold: float
    ) -> float:
        return self.protein_score.optimize_hyperparameters(
            protein_groups, protein_group_peptide_infos, protein_group_fdr_threshold
        )

    def calculate_score(self, score_peptide_pairs) -> float:
        return self.protein_score.calculate_score(score_peptide_pairs)

    def get_score_column(self) -> str:
        return self.protein_score.get_score_column(
            self.score_origin.short_description() == "p"
        )

    def can_do_protein_group_rescue(self) -> bool:
        return self.protein_score.can_do_protein_group_rescue()

    def short_description(self) -> str:
        return (
            self.protein_score.short_description()
            + self.score_origin.short_description()
            + "P"
        )

    def short_description_razor(self) -> str:
        return "rS" if self.use_razor else "dS"

    def long_description(self) -> str:
        return f"{self.protein_score.long_description()} {self.score_origin.long_description()} PEP"

    def long_description_razor(self) -> str:
        return "razor peptides" if self.use_razor else "discard shared peptides"

    def filter_proteins(self, proteins) -> List[str]:
        if self.use_razor:
            return self._retain_protein_with_most_observed_peptides(proteins)
        else:
            return proteins

    def set_peptide_counts_per_protein(
        self, peptide_info_list: PeptideInfoList
    ) -> None:
        if self.use_razor:
            observed_peptides = ObservedPeptides()
            observed_peptides.create(peptide_info_list)
            self.peptide_counts_per_protein = (
                observed_peptides.get_peptide_counts_per_protein()
            )
            self.best_peptide_score_per_protein = (
                observed_peptides.get_best_peptide_score_per_protein()
            )

    def _retain_protein_with_most_observed_peptides(
        self, proteins: List[str]
    ) -> List[str]:
        """Retains only the protein with the most observed peptides.
        Ties are first broken on best scoring (potentially shared) peptide and
        otherwise by the md5 hash of the protein identifier. The latter
        ensures that we randomly select a protein, but that this happens
        consistently across different peptides."""
        num_peptides_per_protein_pairs = [
            (
                self.peptide_counts_per_protein.get(protein, 0),
                -1 * self.best_peptide_score_per_protein.get(protein, 1.0),
                hashlib.md5(protein.encode("utf-8")).hexdigest(),
                protein,
            )
            for protein in proteins
        ]
        pair_with_most_observed_peptides = sorted(
            num_peptides_per_protein_pairs, reverse=True
        )[0]
        protein_with_most_observed_peptides = pair_with_most_observed_peptides[-1]
        return [protein_with_most_observed_peptides]

    def collect_peptide_scores_per_protein(
        self,
        protein_groups: ProteinGroups,
        peptide_info_list: PeptideInfoList,
        peptide_qval_cutoff: float,
        suppress_missing_protein_warning: bool = False,
    ) -> ProteinGroupPeptideInfos:
        """Groups peptides with associated scores by protein

        :param protein_groups: ProteinGroups object
        :param peptide_info_list: Dict of peptide -> (score, proteins)
        :param suppress_missing_protein_warning: suppresses the warning for missing proteins
            in the proteinGroups object. This is set during the rescuing grouping procedure
            since some protein groups will have been filtered out in the rescuing step.
        :returns: lists of (score, peptide, proteins) tuples per protein group
        """
        if not self.get_score_column():
            return self.protein_score.get_protein_scores_from_file()

        logger.info("Assigning peptides to protein groups")
        shared_peptides, unique_peptides = 0, 0
        protein_group_peptide_infos = [list() for _ in range(len(protein_groups))]
        post_err_probs = list()
        for peptide, (score, proteins) in peptide_info_list.items():
            proteins = self.filter_proteins(
                proteins
            )  # filtering for razor peptide approach

            protein_group_idxs = protein_groups.get_protein_group_idxs(proteins)
            if (
                helpers.is_missing_in_protein_groups(protein_group_idxs)
                and not suppress_missing_protein_warning
            ):
                raise Exception(
                    f"Could not find any of the proteins {proteins} in the ProteinGroups object, check if the identifier format is the same. \
                                  1st protein group in ProteinGroups object: {protein_groups.protein_groups[0]}"
                )

            if not self.use_shared_peptides and helpers.is_shared_peptide(
                protein_group_idxs
            ):  # ignore shared peptides
                shared_peptides += 1
                continue

            unique_peptides += 1
            for protein_group_idx in protein_group_idxs:
                protein_group_peptide_infos[protein_group_idx].append(
                    (score, peptide, proteins)
                )

            if not helpers.is_decoy(proteins) and not helpers.is_mbr(score):
                post_err_probs.append(score)

        self.peptide_score_cutoff = fdr.calc_post_err_prob_cutoff(
            post_err_probs, peptide_qval_cutoff
        )
        logger.info(
            f"#Precursors: Shared peptides = {shared_peptides}; Unique peptides = {unique_peptides}"
        )
        return protein_group_peptide_infos
