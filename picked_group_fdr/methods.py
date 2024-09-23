import os
from typing import List
import toml
from dataclasses import dataclass

from .competition import ProteinCompetitionStrategyFactory, ProteinCompetitionStrategy
from .scoring_strategy import ProteinScoringStrategy
from .grouping import ProteinGroupingStrategyFactory, ProteinGroupingStrategy


@dataclass
class MethodConfig:
    picked_strategy: ProteinCompetitionStrategy
    score_type: ProteinScoringStrategy
    grouping_strategy: ProteinGroupingStrategy
    label: str

    def short_description(
        self,
        rescue_step: bool,
        sep: str = "_",
    ) -> str:
        score_label = self.score_type.short_description()
        grouping_label = self.grouping_strategy.short_description(rescue_step=rescue_step)
        razor_label = self.score_type.short_description_razor()
        fdr_label = self.picked_strategy.short_description()
        return f"{score_label}{sep}{grouping_label}{sep}{razor_label}{sep}{fdr_label}"

    def long_description(
        self,
        rescue_step: bool = True,
        sep: str = ", ",
    ) -> str:
        score_label = self.score_type.long_description()
        grouping_label = self.grouping_strategy.long_description(
            rescue_step=rescue_step
        )
        razor_label = self.score_type.long_description_razor()
        fdr_label = self.picked_strategy.long_description()
        return f"{score_label}{sep}{grouping_label}{sep}{razor_label}{sep}{fdr_label}"


def get_methods(method_names: str, use_pseudo_genes: bool) -> List[MethodConfig]:
    """Get the picking, grouping and scoring strategy from TOML files.

    If no method is specified, use the picked protein group method, which was
    the most sensitive well-calibrated method in a benchmark of methods.
    """
    methods = ["picked_protein_group"]
    if method_names:
        methods = method_names.split(",")

    return [parse_method_toml(method, use_pseudo_genes) for method in methods]


def parse_method_toml(method_name: str, use_pseudo_genes: bool) -> MethodConfig:
    """Parse a method configuration from a TOML file.

    The function reads a TOML file containing the configuration of a method
    used for protein group-level estimation. It retrieves parameters such as
    the protein picking strategy, scoring type, grouping strategy, and label
    for the method.

    Picked Strategy Options:
        - PickedStrategy(): Picked FDR
        - ClassicStrategy(): Classic FDR

    Grouping Options:
        - MQNativeGrouping(): MaxQuant grouping from proteinGroups.txt
        - SubsetGrouping(): Emulate protein grouping of MaxQuant based on evidence.txt
                             (currently does not work with simulated datasets since
                             peptideToProteinMap does not contain entrapment labels)
        - NoGrouping(): No protein grouping, each protein is in its own group
        - +Rescued: Rescue protein groups by only considering peptides below, e.g. 1%, protein FDR threshold

    Args:
        method_name (str): The name of the method or path to the TOML file containing the method configuration.
        use_pseudo_genes (bool): Indicates whether pseudo genes should be used in grouping.

    Raises:
        FileNotFoundError: If the method TOML file cannot be found.

    Returns:
        dict: A dictionary containing the parsed method configuration parameters.
            Keys:
                - "pickedStrategy": The picked strategy instance.
                - "scoreType": The scoring strategy instance.
                - "grouping": The grouping strategy instance.
                - "label": The label for the method.
    """
    dir_path = os.path.dirname(os.path.realpath(__file__))
    method_toml_file = os.path.join(dir_path, "methods", f"{method_name}.toml")
    if not os.path.isfile(method_toml_file):
        if method_name.endswith(".toml") and os.path.isfile(method_name):
            method_toml_file = method_name
        else:
            raise FileNotFoundError(
                f"Could not find method {method_name}. Please ensure that it is one of the builtin methods or that it is a path to a TOML file with a .toml file extension."
            )

    method_name = toml.load(method_toml_file)

    if use_pseudo_genes:
        method_name["grouping"] = "pseudo_gene"

    picked_strategy = ProteinCompetitionStrategyFactory(method_name["pickedStrategy"])
    score_type = method_name["scoreType"]
    if method_name["sharedPeptides"] == "razor":
        score_type += " razor"
    scoring_strategy = ProteinScoringStrategy(score_type)
    grouping_strategy = ProteinGroupingStrategyFactory(method_name["grouping"])

    return MethodConfig(
        picked_strategy, scoring_strategy, grouping_strategy, method_name["label"]
    )


def requires_peptide_to_protein_map(method_configs: List[MethodConfig]) -> bool:
    """Check if any configuration requires a peptide-to-protein map."""
    for method_config in method_configs:
        if (
            method_config.grouping_strategy.needs_peptide_to_protein_map()
            or method_config.score_type.remaps_peptides_to_proteins()
        ):
            return True
    return False
