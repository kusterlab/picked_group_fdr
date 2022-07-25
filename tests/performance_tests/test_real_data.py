from picked_group_fdr.picked_group_fdr import main


def run_example():
    data_dir = "/media/processing_results/bierdimpfl/workDir/964"
    #evidence_file = "evidence_test_frozen.txt"
    #evidence_file = "evidence_test_with_extra_random_precursors.txt"
    evidence_file = "evidence_svd_convergence_fail_filtered_with_extra_random_precursors.txt"
    argv = f"--mq_evidence {data_dir}/percolator/{evidence_file} --protein_groups_out {data_dir}/percolator/proteinGroups_test_out.txt --fasta {data_dir}/10090_uniprot-all_canonical_Mouse_20200224.fasta --do_quant --enzyme trypsinp --min-length 7 --max-length 60 --special-aas KR --cleavages 2".split()
    main(argv)


if __name__ == "__main__":
    run_example()
