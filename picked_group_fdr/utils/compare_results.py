import sys

from ..parsers import maxquant


def main(argv):
    protein_groups_file1 = argv[0]
    protein_groups_file2 = argv[1]
    
    protein_groups1 = maxquant.parse_mq_protein_groups_file(protein_groups_file1)
    protein_groups2 = maxquant.parse_mq_protein_groups_file(protein_groups_file2)
    
    leading_protein_groups1 = sorted(list(get_leading_proteins(protein_groups1)), key = lambda x: x[1])
    leading_protein_groups2 = {pg: (q, score) for pg, q, score in get_leading_proteins(protein_groups2)}
    #leading_protein_groups2 = {p: (q, score) for pg, q, score in get_leading_proteins(protein_groups2) for p in pg.split(";")} # for razor approaches
    
    not_found = 0
    for pg, qval1, score1 in leading_protein_groups1:
        qval2, score2 = leading_protein_groups2.get(pg, (-1, -1))
        print(qval1, score1, qval2, score2, pg, sep='\t')
        if qval2 == -1:
            not_found += 1
    print(f"Protein groups not found: {not_found}")


def get_leading_proteins(protein_group_results):
    for pgr in protein_group_results:
        if pgr.qValue > 0.01:
            continue
        proteins = pgr.proteinIds.split(";")
        num_peptides = list(map(int, pgr.peptideCountsUnique.split(";")))
        
        leading_proteins = sorted(protein for protein, n in zip(proteins, num_peptides) if n == max(num_peptides))
        leading_proteins = ";".join(leading_proteins)
        yield leading_proteins, pgr.qValue, pgr.score

if __name__ == "__main__":
    main(sys.argv[1:])
