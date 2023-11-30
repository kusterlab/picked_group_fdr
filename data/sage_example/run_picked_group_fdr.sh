RESULT_DIR=./results

mkdir -p ${RESULT_DIR}

fasta_files="./iprg2016_with_labels.fasta"

sage_results_file="./results.sage.tsv"

python3 -u -m picked_group_fdr \
    --sage_results ${sage_results_file} \
    --protein_groups_out ${RESULT_DIR}/proteinGroups.txt \
    --methods sage

# directly create combined_protein.tsv using lfq.tsv
python3 -u -m picked_group_fdr.pipeline.sage_quantification \
    --fasta ${fasta_files} \
    --sage_results ${sage_results_file} \
    --protein_groups ${RESULT_DIR}/proteinGroups.txt \
    --sage_lfq_tsv ./lfq.tsv \
    --output_folder ${RESULT_DIR}
