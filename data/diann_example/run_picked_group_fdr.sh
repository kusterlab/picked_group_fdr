python3 -u -m picked_group_fdr \
      --fasta ./diann/human_yeast_ecoli.fasta \
      --fasta_use_uniprot_id \
      --diann_reports ./diann/report.parquet \
      --protein_groups_out ./results/report.pg_matrix.tsv \
      --do_quant \
      --lfq_min_peptide_ratios 1 \
      --methods diann


python3 -u -m picked_group_fdr \
      --fasta ./diann/human_yeast_ecoli.fasta \
      --fasta_use_uniprot_id \
      --diann_reports ./diann/report.parquet \
      --protein_groups_out ./results/combined_protein.tsv \
      --output_format fragpipe \
      --do_quant \
      --lfq_min_peptide_ratios 1 \
      --methods diann

python3 -u -m picked_group_fdr \
      --fasta ./diann/human_yeast_ecoli.fasta \
      --fasta_use_uniprot_id \
      --diann_reports ./diann/report.parquet \
      --protein_groups_out ./results/proteinGroups.txt \
      --output_format maxquant \
      --do_quant \
      --lfq_min_peptide_ratios 1 \
      --methods diann
