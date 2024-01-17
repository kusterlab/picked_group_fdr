import os
import toml

from picked_group_fdr import methods

from .competition import ProteinCompetitionStrategyFactory
from .scoring_strategy import ProteinScoringStrategy
from .grouping import ProteinGroupingStrategyFactory


def get_methods(args):
    """Get the picking, grouping and scoring strategy from toml files.

    If no method is specified, use the picked protein group method, which was
    the most sensitive well-calibrated method in a benchmark of methods.

    pickedStrategy: PickedStrategy() = picked FDR
                    ClassicStrategy() = classic FDR

    grouping:       MQNativeGrouping() = MaxQuant grouping from proteinGroups.txt
                    SubsetGrouping() = Emulate protein grouping of MaxQuant based on evidence.txt
                                       (currently does not work with simulated datasets since peptideToProteinMap does not contain entrapment labels)
                    NoGrouping() = No protein grouping, each protein is in its own group
                    +Rescued = Rescue protein groups by only considering peptides below 1% protein FDR threshold
    """
    methods = ["picked_protein_group"]
    if args.methods:
        methods = args.methods.split(",")

    return [parse_method_toml(method) for method in methods]


def parse_method_toml(method: str):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    method_toml_file = os.path.join(dir_path, "methods", f"{method}.toml")
    if not os.path.isfile(method_toml_file):
        if method.endswith(".toml") and os.path.isfile(method):
            method_toml_file = method
        else:
            raise FileNotFoundError(
                f"Could not find method {method}. Please ensure that it is one of the builtin methods or that it is a path to a TOML file with a .toml file extension."
            )

    method = toml.load(method_toml_file)

    picked_strategy = ProteinCompetitionStrategyFactory(method["pickedStrategy"])
    score_type = method["scoreType"]
    if method["sharedPeptides"] == "razor":
        score_type += " razor"
    scoring_strategy = ProteinScoringStrategy(score_type)
    grouping_strategy = ProteinGroupingStrategyFactory(method["grouping"])

    return {
        "pickedStrategy": picked_strategy,
        "scoreType": scoring_strategy,
        "grouping": grouping_strategy,
        "label": method["label"],
    }


def short_description(
    score_type, grouping_strategy, picked_strategy, rescue_step, sep="_"
):
    score_label = score_type.short_description()
    grouping_label = grouping_strategy.short_description(rescue_step=rescue_step)
    razor_label = score_type.short_description_razor()
    fdr_label = picked_strategy.short_description()
    return f"{score_label}{sep}{grouping_label}{sep}{razor_label}{sep}{fdr_label}"


def long_description(
    score_type, grouping_strategy, picked_strategy, rescue_step, sep=", "
):
    score_label = score_type.long_description()
    grouping_label = grouping_strategy.long_description(rescue_step=rescue_step)
    razor_label = score_type.long_description_razor()
    fdr_label = picked_strategy.long_description()
    return f"{score_label}{sep}{grouping_label}{sep}{razor_label}{sep}{fdr_label}"


def get_method_description(config):
    method_description_long = methods.long_description(
        config["scoreType"],
        config["grouping"],
        config["pickedStrategy"],
        rescue_step=True,
    )
    label = config.get("label", "")
    return label, method_description_long
