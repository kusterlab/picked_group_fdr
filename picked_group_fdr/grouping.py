from abc import ABC, abstractmethod
from typing import List
import logging

import numpy as np

from .protein_groups import ProteinGroups
from .observed_peptides import ObservedPeptides

# for type hints only
from .results import ProteinGroupResults
from .peptide_info import PeptideInfoList, ProteinGroupPeptideInfos

logger = logging.getLogger(__name__)


def ProteinGroupingStrategyFactory(method="rescued_subset"):
    methods = {}
    methods["no"] = NoGrouping
    methods["subset"] = SubsetGrouping
    methods["rescued_subset"] = RescuedSubsetGrouping
    methods["mq_native"] = MQNativeGrouping
    methods["rescued_mq_native"] = RescuedMQNativeGrouping
    methods["pseudo_gene"] = PseudoGeneGrouping
    if method not in methods:
        raise ValueError(
            f"Unknown pickedStrategy {method['pickedStrategy']}, should be one of 'no', 'subset' or 'rescued_subset'"
        )
    return methods[method]()


class ProteinGroupingStrategy(ABC):
    @abstractmethod
    def needs_peptide_to_protein_map(self):
        pass

    @abstractmethod
    def group_proteins(
        self, peptide_info_list: PeptideInfoList, mq_protein_groups_file: str
    ) -> ProteinGroups:
        pass

    def get_rescue_steps(self):
        return [False]

    def rescue_protein_groups(
        self,
        peptide_info_list: PeptideInfoList,
        protein_group_results: ProteinGroupResults,
        protein_group_fdr_threshold: float,
        old_protein_groups: ProteinGroups,
        old_protein_group_peptide_infos,
    ) -> ProteinGroups:
        pass

    def update_protein_groups(
        self, protein_groups: ProteinGroups, protein_group_peptide_infos
    ):
        return protein_groups, protein_group_peptide_infos

    @abstractmethod
    def short_description(self, rescue_step):
        pass

    @abstractmethod
    def long_description(self, rescue_step):
        pass


class NoGrouping(ProteinGroupingStrategy):
    def needs_peptide_to_protein_map(self):
        return False

    def group_proteins(
        self, peptide_info_list: PeptideInfoList, mq_protein_groups_file: str
    ) -> ProteinGroups:
        """
        Creates a single protein group for each protein in peptide_info_list

        :param peptide_info_list: dictionary of peptide -> (score, proteins)
        :param mq_protein_groups_file: not used
        :returns: protein groups object
        """
        protein_groups = ProteinGroups()
        seen_proteins = set()
        for _, (_, proteins) in peptide_info_list.items():
            for protein in proteins:
                if protein not in seen_proteins:
                    seen_proteins.add(protein)
                    protein_groups.append([protein])
        protein_groups.create_index()
        return protein_groups

    def short_description(self, rescue_step):
        return "nG"

    def long_description(self, rescue_step):
        return "no protein grouping"


class SubsetGrouping(ProteinGroupingStrategy):
    def needs_peptide_to_protein_map(self):
        return False

    def group_proteins(
        self, peptide_info_list: PeptideInfoList, mq_protein_groups_file: str
    ) -> ProteinGroups:
        logger.info("Grouping proteins (subset strategy):")

        observed_peptides = ObservedPeptides()
        observed_peptides.create(peptide_info_list)

        protein_groups = observed_peptides.generate_protein_groups()
        return protein_groups

    def short_description(self, rescue_step):
        return "sG"

    def long_description(self, rescue_step):
        return "subset protein grouping"


class MQNativeGrouping(ProteinGroupingStrategy):
    """Use grouping provided by a MaxQuant proteinGroups.txt file"""

    def needs_peptide_to_protein_map(self):
        return False

    def group_proteins(
        self, peptide_info_list: PeptideInfoList, mq_protein_groups_file: str
    ) -> ProteinGroups:
        if not mq_protein_groups_file:
            raise ValueError("Missing MQ protein groups file input --mq_protein_groups")

        protein_groups = ProteinGroups.from_mq_protein_groups_file(
            mq_protein_groups_file
        )
        return protein_groups

    def short_description(self, rescue_step):
        return "sG"

    def long_description(self, rescue_step):
        return "subset protein grouping"


class RescuedGrouping:
    score_cutoff: float
    obsolete_protein_groups: ProteinGroups
    obsolete_protein_group_peptide_infos: ProteinGroupPeptideInfos

    def rescue_protein_groups(
        self,
        peptide_info_list: PeptideInfoList,
        protein_group_results: ProteinGroupResults,
        protein_group_fdr_threshold: float,
        old_protein_groups: ProteinGroups,
        old_protein_group_peptide_infos,
    ) -> ProteinGroups:
        self._calculate_rescue_score_cutoff(
            protein_group_results, protein_group_fdr_threshold
        )
        peptide_info_list_filtered = self._filter_peptide_list_by_score_cutoff(
            peptide_info_list
        )
        return self.merge_with_rescued_protein_groups(
            peptide_info_list_filtered,
            old_protein_groups,
            old_protein_group_peptide_infos,
        )

    def get_rescue_steps(self) -> List[bool]:
        return [False, True]

    def short_description(self, rescue_step: bool) -> str:
        if rescue_step:
            return "rsG"
        else:
            return "sG"

    def long_description(self, rescue_step: bool) -> str:
        if rescue_step:
            return "rescued subset protein grouping"
        else:
            return "subset protein grouping"

    def get_rescued_protein_groups(
        self, peptide_info_list: PeptideInfoList
    ) -> ProteinGroups:
        """Generates new protein grouping by "rescuing" protein groups that were split
        due to low-scoring PSMs uniquely mapping to particular isoforms

        :param peptide_info_list: dictionary of peptide -> (score, proteins)
        :returns: updated list of protein groups
        """

        logger.info(
            "Redoing protein grouping using peptides below equivalent protein FDR threshold"
        )

        observed_peptides = ObservedPeptides()
        observed_peptides.create(peptide_info_list)

        updated_protein_groups = observed_peptides.generate_protein_groups()

        connected_protein_graphs = observed_peptides.get_connected_proteins(
            updated_protein_groups
        )
        updated_protein_groups = connected_protein_graphs.decouple_connected_proteins(
            updated_protein_groups
        )

        return updated_protein_groups

    def update_protein_groups(
        self,
        protein_groups: ProteinGroups,
        protein_group_peptide_infos: ProteinGroupPeptideInfos,
    ):
        protein_groups.extend(self.obsolete_protein_groups)
        protein_group_peptide_infos.extend(self.obsolete_protein_group_peptide_infos)
        return protein_groups, protein_group_peptide_infos

    def merge_with_rescued_protein_groups(
        self,
        peptide_info_list: PeptideInfoList,
        protein_groups: ProteinGroups,
        protein_group_peptide_infos,
    ) -> ProteinGroups:
        new_protein_groups = self.get_rescued_protein_groups(peptide_info_list)
        (
            self.obsolete_protein_groups,
            self.obsolete_protein_group_peptide_infos,
        ) = new_protein_groups.add_unseen_protein_groups(
            protein_groups, protein_group_peptide_infos
        )

        logger.info(f"#protein groups before rescue: {len(protein_groups)}")
        logger.info(f"#protein groups after rescue: {len(new_protein_groups)}")

        return new_protein_groups

    def _filter_peptide_list_by_score_cutoff(self, peptide_info_list: PeptideInfoList):
        return {
            peptide: (score, proteins)
            for peptide, (score, proteins) in peptide_info_list.items()
            if score < self.score_cutoff
        }

    def _calculate_rescue_score_cutoff(
        self,
        protein_group_results: ProteinGroupResults,
        protein_group_fdr_threshold: float,
    ):
        """Calculate PEP corresponding to protein_group_fdr_threshold.
        N.B. this only works if the protein score is bestPEP!"""
        identified_protein_scores = [
            pfr.score
            for pfr in protein_group_results
            if pfr.qValue < protein_group_fdr_threshold
        ]

        if len(identified_protein_scores) == 0:
            logger.warning(
                f"Could not calculate rescuing threshold as no proteins were found below {protein_group_fdr_threshold*100:g}% FDR. Setting rescuing threshold to worst scoring protein's score."
            )
            identified_protein_scores = [pfr.score for pfr in protein_group_results]

        self.score_cutoff = np.power(10, min(identified_protein_scores) * -1)
        logger.info(
            f"Rescuing threshold: protein score = {'{0:.3g}'.format(min(identified_protein_scores))}, peptide PEP = {'{0:.3g}'.format(self.score_cutoff)}"
        )


class RescuedSubsetGrouping(RescuedGrouping, SubsetGrouping):
    def long_description(self, rescue_step):
        if rescue_step:
            return "rescued subset protein grouping"
        else:
            return "subset protein grouping"


class RescuedMQNativeGrouping(RescuedGrouping, MQNativeGrouping):
    def long_description(self, rescue_step):
        if rescue_step:
            return "rescued subset protein grouping"
        else:
            return "subset protein grouping"


class PseudoGeneGrouping(ProteinGroupingStrategy):
    """
    Groups proteins that have at least one peptide in common, forming pseudo-genes.
    This is necessary for consistent quantification across runs in the presence of
    isoform specific peptides only detected in a subset of the runs.

    For example:

    Run 1:
    peptide1 -> (isoformA1, isoformA2)
    peptide2 -> (isoformA1)
    peptide3 -> (isoformA2)

    Run 2:
    peptide1 -> (isoformA1, isoformA2)

    The regular protein grouping would create 2 protein groups, one with isoformA1
    and one with isoformA2. However, since Run 2 only has shared peptides for both
    protein groups, it would report neither as detected/quantified even though
    evidence exists for the corresponding gene. By grouping both isoforms in one
    group, it can be detected/quantified in both runs.
    """

    def needs_peptide_to_protein_map(self):
        return True

    def group_proteins(
        self, peptide_info_list: PeptideInfoList, mq_protein_groups_file: str
    ) -> ProteinGroups:
        observed_peptides = ObservedPeptides()
        observed_peptides.create(peptide_info_list)

        new_protein_groups = observed_peptides.generate_protein_groups()

        connected_protein_graphs = observed_peptides.get_connected_proteins(
            new_protein_groups, exclude_identified_proteins=False
        )
        new_protein_groups = connected_protein_graphs.get_connected_proteins(
            new_protein_groups
        )

        return new_protein_groups

    def short_description(self, rescue_step):
        return "gG"

    def long_description(self, rescue_step):
        return "pseudo-gene grouping"
