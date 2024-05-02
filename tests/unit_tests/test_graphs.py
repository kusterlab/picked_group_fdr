import collections

import pytest
import networkx as nx

import picked_group_fdr.graphs as graphs
import picked_group_fdr.protein_groups as protein_groups


class TestConnectedProteinGraphs:
    def test_splitConnectedComponent_noCutPossible(self):
        G = nx.Graph()

        G.add_node("proteinA", node_type="protein")
        G.add_node("proteinB", node_type="protein")
        G.add_node("proteinC", node_type="protein")

        G.add_edge("proteinA", "peptide:proteinA;proteinB")
        G.add_edge("proteinB", "peptide:proteinA;proteinB")
        G.add_edge("proteinA", "peptide:proteinA;proteinC")
        G.add_edge("proteinC", "peptide:proteinA;proteinC")
        G.add_edge("proteinB", "peptide:proteinB;proteinC")
        G.add_edge("proteinC", "peptide:proteinB;proteinC")

        proteins = [
            x
            for x, y in G.nodes(data=True)
            if "node_type" in y and y["node_type"] == "protein"
        ]
        peptides = [x for x, y in G.nodes(data=True) if "node_type" not in y]

        connectedProteinGraphs = graphs.ConnectedProteinGraphs([])
        assert (
            connectedProteinGraphs._split_single_connected_component(
                G, proteins, peptides
            )
            == list()
        )

    def test_splitConnectedComponent_oneCutPossible(self):
        G = nx.Graph()

        G.add_node("proteinA", node_type="protein")
        G.add_node("proteinB", node_type="protein")
        G.add_node("proteinC", node_type="protein")
        G.add_node("proteinD", node_type="protein")

        G.add_edge("proteinA", "peptide:proteinA;proteinB")
        G.add_edge("proteinB", "peptide:proteinA;proteinB")
        G.add_edge("proteinC", "peptide:proteinC;proteinD")
        G.add_edge("proteinD", "peptide:proteinC;proteinD")

        # this link can be removed to create two groups: proteinA+B and proteinC+D
        G.add_edge("proteinB", "peptide:proteinB;proteinC")
        G.add_edge("proteinC", "peptide:proteinB;proteinC")

        proteins = [
            x
            for x, y in G.nodes(data=True)
            if "node_type" in y and y["node_type"] == "protein"
        ]
        peptides = [x for x, y in G.nodes(data=True) if "node_type" not in y]

        connectedProteinGraphs = graphs.ConnectedProteinGraphs([])
        assert set(
            [
                ";".join(sorted(connectedProteinGraphs._get_protein_nodes(x)))
                for x in connectedProteinGraphs._split_single_connected_component(
                    G, proteins, peptides
                )
            ]
        ) == {"proteinA;proteinB", "proteinC;proteinD"}

    def test_decouple_connected_proteins_oneCutPossible(self):
        G = nx.Graph()

        G.add_node("proteinA", node_type="protein")
        G.add_node("proteinB", node_type="protein")
        G.add_node("proteinC", node_type="protein")
        G.add_node("proteinD", node_type="protein")

        G.add_edge("proteinA", "peptide:proteinA;proteinB")
        G.add_edge("proteinB", "peptide:proteinA;proteinB")
        G.add_edge("proteinC", "peptide:proteinC;proteinD")
        G.add_edge("proteinD", "peptide:proteinC;proteinD")

        # this link can be removed to create two groups: proteinA+B and proteinC+D
        G.add_edge("proteinB", "peptide:proteinB;proteinC")
        G.add_edge("proteinC", "peptide:proteinB;proteinC")

        proteins = [
            x
            for x, y in G.nodes(data=True)
            if "node_type" in y and y["node_type"] == "protein"
        ]
        peptides = [x for x, y in G.nodes(data=True) if "node_type" not in y]

        proteinGroups = protein_groups.ProteinGroups(
            [
                ["proteinA", "proteinE"],
                ["proteinB"],
                ["proteinC"],
                ["proteinD", "proteinF"],
            ]
        )
        proteinGroups.create_index()

        subgraphs = collections.deque()
        subgraphs.append(G)

        connectedProteinGraphs = graphs.ConnectedProteinGraphs(subgraphs)
        assert connectedProteinGraphs.decouple_connected_proteins(
            proteinGroups
        ).protein_groups == [
            ["proteinA", "proteinE", "proteinB"],
            ["proteinC", "proteinD", "proteinF"],
        ]


class TestPeptideProteinGraph:
    def test_getConnectedComponents(self):
        G = nx.Graph()

        G.add_node("proteinA", node_type="protein")
        G.add_node("proteinB", node_type="protein")
        G.add_node("proteinC", node_type="protein")
        G.add_node("proteinD", node_type="protein")

        G.add_edge("proteinA", "peptide:proteinA;proteinB")
        G.add_edge("proteinB", "peptide:proteinA;proteinB")
        G.add_edge("proteinC", "peptide:proteinC;proteinD")
        G.add_edge("proteinD", "peptide:proteinC;proteinD")

        connectedProteinGraphs = graphs.ConnectedProteinGraphs([])
        peptideProteinGraph = graphs.PeptideProteinGraph()
        peptideProteinGraph.G = G
        assert set(
            [
                ";".join(sorted(connectedProteinGraphs._get_protein_nodes(x)))
                for x in peptideProteinGraph.get_connected_components()
            ]
        ) == {"proteinA;proteinB", "proteinC;proteinD"}
