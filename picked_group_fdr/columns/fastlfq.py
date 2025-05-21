import networkx as nx
import numpy as np


def build_binary_matrix(samples, peptide_index):
    """
    Create a binary matrix where each row represents a sample and each column represents a peptide.
    The value is 1 if the peptide is present in the sample, else 0.
    """
    n_samples = len(samples)
    n_peptides = len(peptide_index)
    binary_matrix = np.zeros((n_samples, n_peptides), dtype=int)

    for i, sample in enumerate(samples):
        for peptide in sample:
            peptide_id = peptide_index[peptide]  # Get index for peptide
            binary_matrix[i, peptide_id] = 1

    return binary_matrix


def compute_distance_matrix(binary_matrix):
    """
    Compute pairwise distances based on the overlap between samples using binary matrix.
    """
    # Dot product of the binary matrix will give the number of overlapping peptides (intersection)
    overlap_matrix = binary_matrix @ binary_matrix.T  # Matrix multiplication

    # Normalize to get distances (1 / (overlap + small_value))
    np.fill_diagonal(
        overlap_matrix, 0
    )  # Set diagonal to 0 since distance to itself is not needed
    distance_matrix = 1 / (
        overlap_matrix + 1e-6
    )  # Add small value to avoid division by zero
    return distance_matrix


def build_graph(samples):
    """
    Build a graph where nodes are samples and edges are weighted by inverse peptide overlap.
    """
    # Create peptide index to map peptide names to indices in the binary matrix
    all_peptides = set(p for sample in samples for p in sample)
    peptide_index = {peptide: idx for idx, peptide in enumerate(all_peptides)}

    # Create binary matrix representation of samples
    binary_matrix = build_binary_matrix(samples, peptide_index)

    # Compute the pairwise distance matrix (using binary matrix operations)
    distance_matrix = compute_distance_matrix(binary_matrix)

    # Create a graph
    G = nx.Graph()
    n_samples = len(samples)
    for i in range(n_samples):
        G.add_node(i, peptides=samples[i])

    # Add edges with computed distances
    for i in range(n_samples):
        for j in range(i + 1, n_samples):
            dist = distance_matrix[i, j]
            G.add_edge(i, j, weight=dist)

    return G


def prune_graph(G, min_neighbors=3, avg_neighbors=6):
    """
    Prune edges while ensuring:
    - Each node has at least `min_neighbors` edges.
    - Average degree across the graph is around `avg_neighbors`.
    - Graph remains connected.
    """
    # Start with full graph and sort edges by weight (smallest distance first)
    edges_sorted = sorted(G.edges(data=True), key=lambda x: x[2]["weight"])

    # Build minimal required edges (KNN)
    pruned_G = nx.Graph()
    pruned_G.add_nodes_from(G.nodes(data=True))

    for node in G.nodes:
        neighbors = sorted(G[node].items(), key=lambda x: x[1]["weight"])
        for neighbor, attr in neighbors[:min_neighbors]:
            pruned_G.add_edge(node, neighbor, weight=attr["weight"])

    # Add more edges until avg degree is around avg_neighbors
    for u, v, attr in edges_sorted:
        if pruned_G.has_edge(u, v):
            continue
        pruned_G.add_edge(u, v, weight=attr["weight"])
        avg_deg = sum(dict(pruned_G.degree()).values()) / pruned_G.number_of_nodes()
        if avg_deg >= avg_neighbors:
            break

    # Ensure graph is connected
    if not nx.is_connected(pruned_G):
        for u, v, attr in edges_sorted:
            if pruned_G.has_edge(u, v):
                continue
            pruned_G.add_edge(u, v, weight=attr["weight"])
            if nx.is_connected(pruned_G):
                break

    return pruned_G


def export_to_vis_js(G):
    # Copy the nodes and edges into this fiddle: https://jsfiddle.net/apn1rfsc/
    # Generate vis.js nodes
    nodes = ",\n    ".join([f"{{id: {n}, label: 'Node {n}'}}" for n in G.nodes()])

    # Generate vis.js edges
    edges = ",\n    ".join([f"{{from: {u}, to: {v}}}" for u, v in G.edges()])

    # Final output
    vis_js_code = f"""
    // create an array with nodes
    var nodes = new vis.DataSet([
        {nodes}
    ]);

    // create an array with edges
    var edges = new vis.DataSet([
        {edges}
    ]);
    """

    print(vis_js_code)


if __name__ == "__main__":
    # Example usage:
    import random

    # Generate a pool of 100 unique peptides
    all_peptides = [f"peptide{i}" for i in range(1, 1001)]

    # Generate 20 samples, each with 10 to 30 peptides
    samples = [
        set(random.sample(all_peptides, random.randint(10, 300)))
        for _ in range(2000)
    ]

    print("Creating and pruning graph")

    G_full = build_graph(samples)
    G_pruned = prune_graph(G_full, min_neighbors=3, avg_neighbors=6)

    export_to_vis_js(G_pruned)
