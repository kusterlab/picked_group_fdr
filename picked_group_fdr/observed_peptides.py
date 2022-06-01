import collections
from typing import List, Set, Dict, Tuple
import logging

from . import graphs
from .protein_groups import ProteinGroups


logger = logging.getLogger(__name__)


class ObservedPeptides:
    peptide_to_proteins_dict: Dict[str, List[str]]
    protein_to_peptides_dict: Dict[str, List[str]]
    
    def __init__(self) -> None:
        self.peptide_to_proteins_dict = dict()
        self.protein_to_peptides_dict = collections.defaultdict(list) 
    
    def create(self, peptideInfoList: Dict[str, Tuple[float, List[str]]]) -> None:
        """Populates two dictionaries, the first mapping peptides to proteins and the
        second mapping proteins to peptides
        
        :param peptideInfoList: dictionary of peptide -> (score, proteins)
        """
        for peptide, (score, proteins) in peptideInfoList.items():
            self.peptide_to_proteins_dict[peptide] = proteins
            for protein in proteins:
                self.protein_to_peptides_dict[protein].append(peptide)
    
    def get_proteins(self, peptide: str) -> List[str]:
        """Returns list of proteins that a peptide matches to"""
        return self.peptide_to_proteins_dict.get(peptide, [])
    
    def get_peptides(self, protein: str) -> List[str]:
        """Returns list of peptides that a protein matches to"""
        return self.protein_to_peptides_dict[protein]
    
    def get_peptide_counts_per_protein(self) -> Dict[str, int]:
        return { protein : len(set(peptides)) for protein, peptides in self.protein_to_peptides_dict.items() }

    def generate_protein_groups(self) -> ProteinGroups:
        """Generating protein groups from a dictionary of observed_peptide -> protein"""
        logger.info("Generating protein groups from observed peptides")
        
        proteinGroups = ProteinGroups.from_observed_peptide_map(self.protein_to_peptides_dict)
        proteinGroups.create_index()
        
        for protein, peptides in self.protein_to_peptides_dict.items():
            superset_proteins = self._get_superset_proteins(peptides)
            
            # add to the protein with the most peptides out of the proteins that this protein is a subset of
            superset_proteins = sorted(superset_proteins, key = lambda x : len(self.protein_to_peptides_dict[x]), reverse = True)
            for superset_protein in superset_proteins:
                if len(proteinGroups.get_protein_group(superset_protein, check_idx_valid = False)) > 0 and protein != superset_protein:
                    proteinGroups.merge_groups(superset_protein, protein)
                    break
        
        proteinGroups.remove_empty_groups()
        
        return proteinGroups

    def _get_superset_proteins(self, peptides: List[str]) -> List[str]:
        """Get all proteins whose set of peptides forms a (non-strict) superset of the peptides for this protein"""
        candidateProteins = self.get_proteins(peptides[0])
        for peptide in peptides[1:]:
            proteins = set(self.get_proteins(peptide))
            candidateProteins = [p for p in candidateProteins if p in proteins] # this is 15% faster than set.intersection()
            
            # if there is only 1 protein candidate left, the only protein that forms a 
            # superset is the protein itself. No need to look further, so break and return.
            if len(candidateProteins) == 1:
                break
        return candidateProteins

    def get_connected_proteins(self, proteinGroups, exclude_identified_proteins = True):
        """ Combines protein groups for which no identified peptide is left after
        the removal of shared peptides. 
        
        This helps with cases such as:
        
        peptide1 => proteinA, proteinB
        peptide2 => proteinB, proteinC
        peptide3 => proteinA, proteinC
        
        Since none of the proteins form a subset of another protein, they will each
        form a separate protein group, none of which has a unique peptide. This
        function looks for such connected components and combines proteinA, proteinB
        and proteinC into a protein group.
        
        Protein groups with identified peptides are excluded from this process. For
        example, if proteinC had a unique peptide, only proteinA and proteinB would
        be grouped.
        
        For larger connected components, we iteratively try to remove shared peptides
        such that we have two or more connected components each with identified
        peptides.
        
        This helps with cases such as:
        
        peptide1 => proteinA, proteinB
        peptide2 => proteinB, proteinC
        peptide3 => proteinC, proteinD
        
        Here, we would remove peptide2 and still have 2 identified protein groups
        (proteinA;proteinB and proteinC;proteinD).

        :param proteinGroups: list of protein groups, where each entry is itself a 
            list of proteins in that protein group
        :param exclude_identified_proteins: do not add proteins with unique peptides
            to the graph. Set this to False for pseudo-gene grouping.
        :returns: updated list of protein groups where connected components are added
            as protein groups
        """
        identifiedProteinGroupIdxs = set()
        if exclude_identified_proteins:
            identifiedProteinGroupIdxs = self._get_protein_group_idxs_with_unique_peptides(proteinGroups)
        
        bipartiteGraph = graphs.PeptideProteinGraph()
        bipartiteGraph.create_graph(
                proteinGroups, identifiedProteinGroupIdxs, self)
        
        subgraphs = bipartiteGraph.get_connected_components()
        return graphs.ConnectedProteinGraphs(subgraphs)

    def _get_protein_group_idxs_with_unique_peptides(self, proteinGroups: ProteinGroups) -> Set[int]:
        """
        returns protein group idxs of protein groups with >=1 unique peptide(s)
        """
        identifiedProteinGroupIdxs = set()
        for peptide, proteins in self.peptide_to_proteins_dict.items():
            proteinGroupIdxsForPeptide = proteinGroups.get_protein_group_idxs(proteins)
            if len(proteinGroupIdxsForPeptide) == 1:
                identifiedProteinGroupIdxs.add(list(proteinGroupIdxsForPeptide)[0])
        return identifiedProteinGroupIdxs
