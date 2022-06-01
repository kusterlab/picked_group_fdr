import collections
import itertools
from typing import List, Set, Dict
import logging

from networkx.algorithms.connectivity import build_auxiliary_node_connectivity, minimum_st_node_cut
from networkx.algorithms.flow import build_residual_network
import networkx as nx


logger = logging.getLogger(__name__)


class PeptideProteinGraph:
  """Creates a bipartite graph between proteins and pseudo peptides.
    
  Instead of using the actual peptides, we represent each peptide by all the
  proteins it maps to. In this way, peptides shared between the same proteins
  only create a single pseudo-peptide node.
  """
  
  G: nx.Graph
  
  def __init__(self):
    self.G = nx.Graph()
    
  def create_graph(self, proteinGroups, identifiedProteinGroupIdxs,
                   observedPeptideProteinMaps):
    """Creates the bipartite graph.
    
    :param proteinGroups: ProteinGrouping object
    :param identifiedProteinGroupIdxs: set of protein group idxs that have a
      unique peptide already
    :param observedProteinPeptideMaps: dictionary of protein => observed peptides 
      below the PSM-level cutoff
    :returns: updated list of protein groups where connected components are added
      as protein groups
    """
    for idx, pg in enumerate(proteinGroups):
      # no need to add protein groups to graph that already have a unique peptide
      if idx in identifiedProteinGroupIdxs:
        continue
      
      protein = pg[0]
      self.G.add_node(protein, nodeType = "protein")
      for peptide in observedPeptideProteinMaps.get_peptides(protein):
        proteinsForPeptide = observedPeptideProteinMaps.get_proteins(peptide)
        leadingProteins = proteinGroups.get_leading_proteins(proteinsForPeptide)
        # add_edge creates nodes for the pseudo peptides automatically
        self.G.add_edge(protein, "peptide:" + ";".join(sorted(leadingProteins)))
  
  def get_connected_components(self):
    subgraphs = collections.deque()
    for c in nx.connected_components(self.G):
      subgraphs.append(self.G.subgraph(c).copy())
    return subgraphs
  
  def print_summary(self):
    numProteins = len([x for x,y in self.G.nodes(data=True) if 'nodeType' in y and y['nodeType'] == "protein"])
    logger.info(f"#Proteins: {numProteins}, #Nodes: {len(self.G)}, #Edges: {self.G.number_of_edges()}")


class ConnectedProteinGraphs:
  subgraphs: List[nx.Graph]
  
  def __init__(self, subgraphs):
    self.subgraphs = subgraphs
  
  def decouple_connected_proteins(self, proteinGroups):
    while len(self.subgraphs) > 0:
      c = self.subgraphs.popleft()
      proteins = self._get_protein_nodes(c)
      peptides = self._get_peptide_nodes(c)
      
      subgraphsLocal = []
      if len(proteins) > 1:
        subgraphsLocal = self._split_single_connected_component(c, proteins, peptides)
      
      if len(subgraphsLocal) == 0:
        # the connected component couldn't be split into multiple connected components by removing peptides
        leadingProtein = proteins[0]
        for protein in proteins[1:]:
          proteinGroups.merge_groups(leadingProtein, protein)
      else:
        for subgraph in subgraphsLocal:
          self.subgraphs.append(subgraph)
    
    proteinGroups.remove_empty_groups()
    
    return proteinGroups
  
  def get_connected_proteins(self, proteinGroups):
    while len(self.subgraphs) > 0:
      c = self.subgraphs.popleft()
      proteins = self._get_protein_nodes(c)
      peptides = self._get_peptide_nodes(c)
      
      leadingProtein = proteins[0]
      for protein in proteins[1:]:
        proteinGroups.merge_groups(leadingProtein, protein)
      
    proteinGroups.remove_empty_groups()
    
    return proteinGroups

  # adapted from https://stackoverflow.com/questions/47146755/networkx-find-all-minimal-cuts-consisting-of-only-nodes-from-one-set-in-a-bipar
  def _split_single_connected_component(self, G, A, B):  
    # build auxiliary networks
    H = build_auxiliary_node_connectivity(G)
    R = build_residual_network(H, 'capacity')

    # get all cuts that consist of nodes exclusively from B which disconnet
    # nodes from A
    best_cut = []
    seen = []
    
    cur_min = float('inf')
    for s, t in itertools.combinations(A,2):
      cut = minimum_st_node_cut(G, s, t, auxiliary=H, residual=R)
      if len(cut) < cur_min and set(cut).issubset(B) and cut not in seen:
        subGraphs = self._get_subgraphs_after_cut(G, cut)
        # check if each component after the cut has at least 2 proteins
        disconnectedComponents = list(filter(lambda x : len(x) == 1, subGraphs))
        if len(disconnectedComponents) == 0:
          best_cut = subGraphs
          if len(cut) == 1:
            break
          seen.append(cut)
    
    return best_cut

  def _get_protein_nodes(self, G):
    return [x for x,y in G.nodes(data=True) if 'nodeType' in y and y['nodeType'] == "protein"]

  def _get_peptide_nodes(self, G):
    return [x for x,y in G.nodes(data=True) if 'nodeType' not in y or y['nodeType'] != "protein"]

  def _get_subgraphs_after_cut(self, G, cut):
    G2 = G.copy()
    G2.remove_nodes_from(cut)
    return [G2.subgraph(c).copy() for c in nx.connected_components(G2)]
    
