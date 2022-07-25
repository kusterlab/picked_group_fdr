from picked_group_fdr.picked_group_fdr import main


def run_example():
    argv = f"--mq_evidence /home/matthewt/git/picked_group_fdr/data/itraq_example//out/percolator/evidence.txt --protein_groups_out /home/matthewt/git/picked_group_fdr/data/itraq_example//out/percolator/proteinGroups.txt --do_quant --fasta /home/matthewt/git/picked_group_fdr/data/itraq_example//human.fasta --enzyme trypsinp --min-length 7 --max-length 60 --special-aas KR --cleavages 2".split()
    main(argv)


if __name__ == "__main__":
    run_example()
