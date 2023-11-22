RESULT_DIR=./results

mkdir -p ${RESULT_DIR}

fasta_files="./iprg2016_with_labels.fasta"

fragpipe_psm_files="./A1/psm.tsv ./A2/psm.tsv ./A3/psm.tsv \
    ./B1/psm.tsv ./B2/psm.tsv ./B3/psm.tsv \
    ./C1/psm.tsv ./C2/psm.tsv ./C3/psm.tsv"

python3 -u -m picked_group_fdr \
    --fasta ${fasta_files} \
    --fragpipe_psm ${fragpipe_psm_files} \
    --protein_groups_out ${RESULT_DIR}/proteinGroups.txt \
    --methods fragpipe

python3 -u -m picked_group_fdr.pipeline.update_fragpipe_results \
    --fasta ${fasta_files} \
    --fragpipe_psm ${fragpipe_psm_files} \
    --protein_groups ${RESULT_DIR}/proteinGroups.txt \
    --output_folder ${RESULT_DIR}
