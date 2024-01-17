from typing import Dict, List, Tuple

# TODO: transform into proper class
PeptideInfoList = Dict[str, Tuple[float, List[str]]]

# for each protein group, a list of tuples (posterior_error_probability, peptide_sequence, proteins)
# TODO: transform into proper class
ProteinGroupPeptideInfos = List[List[Tuple[float,str,List[str]]]]
