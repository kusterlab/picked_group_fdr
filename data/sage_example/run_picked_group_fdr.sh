RESULT_DIR=./results

mkdir -p ${RESULT_DIR}

fasta_files="./iprg2016_with_labels.fasta"

sage_results_file="./results.sage.tsv"
sage_lfq_file="./lfq.tsv"

# Option 1: directly create combined_protein.tsv (FragPipe format) using lfq.tsv
python3 -u -m picked_group_fdr \
    --fasta ${fasta_files} \
    --sage_results ${sage_results_file} \
    --sage_lfq_tsv ${sage_lfq_file} \
    --protein_groups_out ${RESULT_DIR}/combined_protein.tsv \
    --output_format fragpipe \
    --do_quant \
    --lfq_min_peptide_ratios 1 \
    --methods sage

# Option 2: directly create proteinGroups.txt (MQ format) using lfq.tsv
python3 -u -m picked_group_fdr \
    --fasta ${fasta_files} \
    --sage_results ${sage_results_file} \
    --sage_lfq_tsv ${sage_lfq_file} \
    --protein_groups_out ${RESULT_DIR}/proteinGroups_with_quant.tsv \
    --do_quant \
    --methods sage

# Option 3: create proteinGroups.txt (MQ format) without quantification
python3 -u -m picked_group_fdr \
    --sage_results ${sage_results_file} \
    --protein_groups_out ${RESULT_DIR}/proteinGroups.txt \
    --methods sage

# Option 4: use intermediate proteinGroups.txt file (from Option 3) to create combined_protein.tsv (FragPipe format) using lfq.tsv
python3 -u -m picked_group_fdr.pipeline.sage_quantification \
    --fasta ${fasta_files} \
    --sage_results ${sage_results_file} \
    --protein_groups ${RESULT_DIR}/proteinGroups.txt \
    --sage_lfq_tsv ${sage_lfq_file} \
    --protein_groups_out ${RESULT_DIR}/combined_protein.tsv \
    --output_format fragpipe

# Option 5: use intermediate proteinGroups.txt (from Option 3) file to create proteinGroups.txt (MQ format) using lfq.tsv
python3 -u -m picked_group_fdr.pipeline.sage_quantification \
    --fasta ${fasta_files} \
    --sage_results ${sage_results_file} \
    --protein_groups ${RESULT_DIR}/proteinGroups.txt \
    --sage_lfq_tsv ${sage_lfq_file} \
    --protein_groups_out ${RESULT_DIR}/proteinGroups_with_quant.tsv