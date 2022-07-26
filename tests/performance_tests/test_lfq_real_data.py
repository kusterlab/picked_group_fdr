from picked_group_fdr.quantification import main
import sys

print(sys.version)

def run_example_776():
    data_dir = "/media/processing_results/bierdimpfl/workDir/776"
    #data_dir = "/root/data"
    #evidence_file = "evidence_Q8WZ42.txt"
    evidence_file = "evidence_top10.txt"
    #evidence_file = "evidence.txt"
    argv = f"--mq_evidence {data_dir}/{evidence_file} --mq_protein_groups {data_dir}/percolator/proteinGroups.txt --protein_groups_out {data_dir}/percolator/proteinGroups_test.txt --fasta {data_dir}/uniprot_HomoSapiens_Swissprot_Canonical_20210528.fasta".split()
    main(argv)


def run_example():
    data_dir = "/media/processing_results/bierdimpfl/workDir/964"
    #evidence_file = "evidence_test_with_extra_random_precursors.txt"
    #protein_groups_file = "proteinGroups_test.txt"
    evidence_file = "evidence_svd_convergence_fail_filtered.txt"
    protein_groups_file = "proteinGroups_svd_convergence_fail.txt"
    argv = f"--mq_evidence {data_dir}/percolator/{evidence_file} --mq_protein_groups {data_dir}/percolator/{protein_groups_file} --protein_groups_out {data_dir}/percolator/proteinGroups_test_out.txt --fasta {data_dir}/10090_uniprot-all_canonical_Mouse_20200224.fasta --fasta_use_uniprot_id".split()
    main(argv)


if __name__ == "__main__":
    run_example()
