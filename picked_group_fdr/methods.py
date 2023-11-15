import os
import toml

from .competition import ProteinCompetitionStrategyFactory
from .scoring import ProteinScoringStrategy
from .grouping import ProteinGroupingStrategyFactory
    

def get_methods(args):
    """
    pickedStrategy: PickedStrategy() = picked FDR
                    ClassicStrategy() = classic FDR
    
    grouping:       MQNativeGrouping() = MaxQuant grouping from proteinGroups.txt
                    SubsetGrouping() = Emulate protein grouping of MaxQuant based on evidence.txt
                                       (currently does not work with simulated datasets since peptideToProteinMap does not contain entrapment labels)
                    NoGrouping() = No protein grouping, each protein is in its own group
                    +Rescued = Rescue protein groups by only considering peptides below 1% protein FDR threshold
    """
    configs = list()
    
    ### WELL-CALIBRATED METHODS
    
    if args.methods:
        for method in args.methods.split(","):
            configs.append(parse_method_toml(method))
        return configs
    
    # final method
    configs.append(parse_method_toml('picked_protein_group'))    
    
    return configs


def parse_method_toml(method: str):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    method_toml_file = os.path.join(dir_path, 'methods', f'{method}.toml')
    if not os.path.isfile(method_toml_file):
        if method.endswith(".toml") and os.path.isfile(method):
            method_toml_file = method
        else:
            raise FileNotFoundError(f"Could not find method {method}. Please ensure that it is one of the builtin methods or that it is a path to a TOML file with a .toml file extension.")

    method = toml.load(method_toml_file)
    
    picked_strategy = ProteinCompetitionStrategyFactory(method['pickedStrategy'])
    score_type = method['scoreType']
    if method['sharedPeptides'] == "razor":
        score_type += " razor"
    scoring_strategy = ProteinScoringStrategy(score_type)
    grouping_strategy = ProteinGroupingStrategyFactory(method['grouping'])
        
    return {'pickedStrategy': picked_strategy, 'scoreType': scoring_strategy, 'grouping': grouping_strategy, 'label': method['label'] }


def short_description(scoreType, groupingStrategy, pickedStrategy, rescue_step, sep='_'):
    scoreLabel = scoreType.short_description()
    groupingLabel = groupingStrategy.short_description(rescue_step = rescue_step)
    razorLabel = scoreType.short_description_razor()
    fdrLabel = pickedStrategy.short_description()
    return f"{scoreLabel}{sep}{groupingLabel}{sep}{razorLabel}{sep}{fdrLabel}"


def long_description(scoreType, groupingStrategy, pickedStrategy, rescue_step, sep=', '):
    scoreLabel = scoreType.long_description()
    groupingLabel = groupingStrategy.long_description(rescue_step = rescue_step)
    razorLabel = scoreType.long_description_razor()
    fdrLabel = pickedStrategy.long_description()
    return f"{scoreLabel}{sep}{groupingLabel}{sep}{razorLabel}{sep}{fdrLabel}"

