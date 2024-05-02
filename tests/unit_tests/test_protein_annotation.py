import picked_group_fdr.protein_annotation as protein_annotation


def test_has_gene_names_true():
    proteinAnnotations = dict()
    for p in ["A", "B", "C", "D", "E"]:
        proteinAnnotations[f"protein{p}"] = protein_annotation.ProteinAnnotation(
            id=f"protein{p}", gene_name=f"gene{p}", fasta_header=f"fasta{p}"
        )
    assert (
        protein_annotation.has_gene_names(proteinAnnotations, min_ratio_with_genes=0.5)
        == True
    )


def test_has_gene_names_false():
    proteinAnnotations = dict()
    for p in ["A", "B", "C", "D", "E"]:
        proteinAnnotations[f"protein{p}"] = protein_annotation.ProteinAnnotation(
            id=f"protein{p}", gene_name="", fasta_header=f"fasta{p}"
        )
    assert (
        protein_annotation.has_gene_names(proteinAnnotations, min_ratio_with_genes=0.5)
        == False
    )


def test_has_gene_names_false_below_half():
    proteinAnnotations = dict()
    for p in ["A", "B", "C", "D"]:
        proteinAnnotations[f"protein{p}"] = protein_annotation.ProteinAnnotation(
            id=f"protein{p}", gene_name=f"gene{p}", fasta_header=f"fasta{p}"
        )
    for p in ["A", "B", "C", "D", "E"]:
        proteinAnnotations[f"protein_no_gene{p}"] = (
            protein_annotation.ProteinAnnotation(
                id=f"protein_no_gene{p}",
                gene_name="",
                fasta_header=f"fasta_no_gene{p}",
            )
        )

    assert (
        protein_annotation.has_gene_names(proteinAnnotations, min_ratio_with_genes=0.5)
        == False
    )
