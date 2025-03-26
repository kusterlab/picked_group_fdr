from abc import ABC, abstractmethod
import logging

import numpy as np


from . import helpers
from . import fdr
from .parsers import protein_groups as pgp

logger = logging.getLogger(__name__)


"""
N.B.: Higher protein score indicates more confident identification

scoreType:  MQ_protein = MaxQuant protein score from proteinGroups.txt (see below)
            multPEP = emulate MaxQuant protein score with multPEP with dividing of PEPs by fixed value
            bestPEP = best scoring peptide based on PEP
            Andromeda = raw Andromeda score from evidence.txt

razor: True = use Occam's razor to decide which protein to assign to a shared peptide

useSharedPeptides: True = use peptides shared by multiple proteins / protein groups
                   False = discard peptides shared by multiple proteins / protein groups
"""


class ProteinScore(ABC):
    @abstractmethod
    def calculate_score(self, scorePeptidePairs):
        pass

    @abstractmethod
    def can_do_protein_group_rescue(self):
        pass

    @abstractmethod
    def get_score_column(self, percolator_input):
        pass

    @abstractmethod
    def short_description(self):
        pass

    def optimize_hyperparameters(
        self,
        protein_groups,
        protein_group_peptide_infos,
        protein_group_fdr_threshold: float,
    ):
        pass


class MQProteinScore(ProteinScore):
    mq_protein_groups_file: str

    def __init__(self, mq_protein_groups_file):
        self.mq_protein_groups_file = mq_protein_groups_file

    def calculate_score(self, scorePeptidePairs):
        raise NotImplementedError

    def can_do_protein_group_rescue(self):
        return False

    def get_score_column(self, percolator_input):
        return None

    def short_description(self):
        return "m"

    def long_description(self):
        return "multiplication of"

    def get_protein_scores_from_file(self):
        protein_group_peptide_infos = list()
        for protein_group, protein_score in pgp.parse_protein_groups_file_single(
            self.mq_protein_groups_file
        ):
            protein_group_peptide_infos.append([(protein_score, "NA", protein_group)])
        return protein_group_peptide_infos


class BestAndromedaScore(ProteinScore):
    def calculate_score(self, score_peptide_pairs):
        return (
            max([y[0] for y in score_peptide_pairs])
            if len(score_peptide_pairs) > 0
            else -100.0
        )

    def can_do_protein_group_rescue(self):
        return False

    def get_score_column(self, percolator_input):
        return "score"

    def short_description(self):
        return "b"

    def long_description(self):
        return "best"


class BestPEPScore(ProteinScore):
    def calculate_score(self, score_peptide_pairs):
        # return max([-1*y[0] for y in scorePeptidePairs]) if len(scorePeptidePairs) > 0 else -100.0
        return (
            max([-1 * np.log10(y[0] + np.nextafter(0, 1)) for y in score_peptide_pairs])
            if len(score_peptide_pairs) > 0
            else -100.0
        )

    def can_do_protein_group_rescue(self):
        return True

    def get_score_column(self, percolator_input):
        if percolator_input:
            return "posterior_error_prob"
        else:
            return "pep"

    def short_description(self):
        return "b"

    def long_description(self):
        return "best"


class MultPEPScore(ProteinScore):
    div: float

    def __init__(self):
        self.div = 1.0

    def optimize_hyperparameters(
        self, protein_groups, protein_group_peptide_infos, protein_group_fdr_threshold
    ):
        protein_score_tuples = list()
        for protein_group, protein_group_score_list in zip(
            protein_groups, protein_group_peptide_infos
        ):
            mult_pep, num_peptides = self._get_score_and_num_peptides(
                protein_group_score_list
            )
            if num_peptides > 0 and not np.isnan(mult_pep):
                protein_score_tuples.append(
                    [mult_pep, num_peptides, helpers.is_decoy(protein_group)]
                )

        protein_score_tuples = np.array(protein_score_tuples)
        min_range, max_range = 0.1, 1.0
        for step_size in [1e-1, 1e-2, 1e-3, 1e-4]:
            num_identified_targets, div = self._get_optimal_div(
                protein_score_tuples,
                np.arange(min_range, max_range, step_size),
                protein_group_fdr_threshold,
            )
            min_range, max_range = div - step_size * 0.9, div + step_size * 0.9

        self.div = div
        logger.info(
            f"Optimal division factor for multPEP score: {self.div:.4f} (#targets = {num_identified_targets})"
        )

    def _get_optimal_div(
        self,
        protein_score_tuples,
        div_range: np.array,
        protein_group_fdr_threshold: float,
    ):
        most_identified_targets = (-np.inf, 1.0)
        for div in div_range:
            protein_scores = protein_score_tuples[:, 0] + protein_score_tuples[
                :, 1
            ] * np.log10(div)

            sorted_idxs = np.argsort(protein_scores)[::-1]
            num_decoys = protein_score_tuples[:, 2][sorted_idxs].cumsum()
            num_targets = (-protein_score_tuples[:, 2] + 1)[sorted_idxs].cumsum()

            fdrs = np.divide(num_decoys + 1, num_targets + 1)
            qvals = fdr.fdrs_to_qvals(fdrs)
            num_identified_targets = fdr.count_below_threshold(
                qvals,
                protein_group_fdr_threshold,
                protein_score_tuples[:, 2][sorted_idxs],
            )
            if num_identified_targets > most_identified_targets[0]:
                most_identified_targets = (num_identified_targets, div)
        return most_identified_targets

    def calculate_score(self, score_peptide_pairs):
        mult_pep, num_peptides = self._get_score_and_num_peptides(score_peptide_pairs)
        mult_pep += np.log10(self.div) * num_peptides
        if num_peptides == 0 or np.isnan(mult_pep):
            mult_pep = -100.0
        return mult_pep

    def _get_score_and_num_peptides(self, score_peptide_pairs):
        mult_pep = 0.0
        seen_peptides = set()
        for PEP, peptide, _ in sorted(score_peptide_pairs):
            if peptide not in seen_peptides:
                seen_peptides.add(peptide)
                mult_pep -= np.log10(PEP + np.nextafter(0, 1))
        return mult_pep, len(seen_peptides)

    def can_do_protein_group_rescue(self):
        return True

    def get_score_column(self, percolator_input):
        if percolator_input:
            return "posterior_error_prob"
        else:
            return "pep"

    def short_description(self):
        return "m"

    def long_description(self):
        return "multiplication of"
