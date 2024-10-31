import collections
import itertools
from typing import List
import logging

from networkx.algorithms.connectivity import (
    build_auxiliary_node_connectivity,
    minimum_st_node_cut,
)
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

    def create_graph(
        self,
        protein_groups,
        identified_protein_group_idxs,
        observed_peptide_protein_maps,
    ):
        """Creates the bipartite graph.

        :param protein_groups: ProteinGrouping object
        :param identified_protein_group_idxs: set of protein group idxs that have a
            unique peptide already
        :param observed_peptide_protein_maps: dictionary of protein => observed peptides
            below the PSM-level cutoff
        :returns: updated list of protein groups where connected components are added
            as protein groups
        """
        for idx, pg in enumerate(protein_groups):
            # no need to add protein groups to graph that already have a unique peptide
            if idx in identified_protein_group_idxs:
                continue

            protein = pg[0]
            self.G.add_node(protein, node_type="protein")
            for peptide in observed_peptide_protein_maps.get_peptides(protein):
                proteins_for_peptide = observed_peptide_protein_maps.get_proteins(
                    peptide
                )
                leading_proteins = protein_groups.get_leading_proteins(
                    proteins_for_peptide
                )
                # add_edge creates nodes for the pseudo peptides automatically
                self.G.add_edge(
                    protein, "peptide:" + ";".join(sorted(leading_proteins))
                )

    def get_connected_components(self):
        subgraphs = collections.deque()
        for c in nx.connected_components(self.G):
            subgraphs.append(self.G.subgraph(c).copy())
        return subgraphs

    def print_summary(self):
        num_proteins = len(
            [
                x
                for x, y in self.G.nodes(data=True)
                if "node_type" in y and y["node_type"] == "protein"
            ]
        )
        logger.info(
            f"#Proteins: {num_proteins}, #Nodes: {len(self.G)}, #Edges: {self.G.number_of_edges()}"
        )


class ConnectedProteinGraphs:
    subgraphs: List[nx.Graph]

    def __init__(self, subgraphs):
        self.subgraphs = subgraphs

    def decouple_connected_proteins(self, protein_groups):
        while len(self.subgraphs) > 0:
            c = self.subgraphs.popleft()
            proteins = self._get_protein_nodes(c)
            peptides = self._get_peptide_nodes(c)

            subgraphs_local = []
            if len(proteins) > 1:
                subgraphs_local = self._split_single_connected_component(
                    c, proteins, peptides
                )

            if len(subgraphs_local) == 0:
                # the connected component couldn't be split into multiple connected components by removing peptides
                leading_protein = proteins[0]
                for protein in proteins[1:]:
                    protein_groups.merge_groups(leading_protein, protein)
            else:
                for subgraph in subgraphs_local:
                    self.subgraphs.append(subgraph)

        protein_groups.remove_empty_groups()

        return protein_groups

    def get_connected_proteins(self, protein_groups):
        while len(self.subgraphs) > 0:
            c = self.subgraphs.popleft()
            proteins = self._get_protein_nodes(c)

            leading_protein = proteins[0]
            for protein in proteins[1:]:
                protein_groups.merge_groups(leading_protein, protein)

        protein_groups.remove_empty_groups()

        return protein_groups

    # adapted from https://stackoverflow.com/questions/47146755/networkx-find-all-minimal-cuts-consisting-of-only-nodes-from-one-set-in-a-bipar
    def _split_single_connected_component(
        self, G: nx.Graph, A: List[str], B: List[str]
    ):
        # build auxiliary networks
        H = build_auxiliary_node_connectivity(G)
        R = build_residual_network(H, "capacity")

        # get all cuts that consist of nodes exclusively from B which disconnect
        # nodes from A
        best_cut = []
        seen = []

        cur_min = float("inf")
        for s, t in itertools.combinations(A, 2):
            cut = minimum_st_node_cut(G, s, t, auxiliary=H, residual=R)
            if len(cut) < cur_min and set(cut).issubset(B) and cut not in seen:
                sub_graphs = self._get_subgraphs_after_cut(G, cut)
                # check if each component after the cut has at least 2 proteins
                disconnected_components = list(
                    filter(lambda x: len(x) == 1, sub_graphs)
                )
                if len(disconnected_components) == 0:
                    best_cut = sub_graphs
                    if len(cut) == 1:
                        break
                    seen.append(cut)

        return best_cut

    def _get_protein_nodes(self, G: nx.Graph):
        return sorted(
            x
            for x, y in G.nodes(data=True)
            if "node_type" in y and y["node_type"] == "protein"
        )

    def _get_peptide_nodes(self, G: nx.Graph):
        return sorted(
            x
            for x, y in G.nodes(data=True)
            if "node_type" not in y or y["node_type"] != "protein"
        )

    def _get_subgraphs_after_cut(self, G: nx.Graph, cut):
        G2 = G.copy()
        G2.remove_nodes_from(cut)
        return [G2.subgraph(c).copy() for c in nx.connected_components(G2)]
