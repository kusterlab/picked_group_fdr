DATA_DIR="data/lfq_example"
RESULT_DIR="tests/system_tests/test_pipeline"

mkdir -p ${RESULT_DIR}

python -u -m picked_group_fdr.pipeline.andromeda2pin \
    ${DATA_DIR}/evidence.txt \
    --outputTab ${RESULT_DIR}/andromeda.tab \
    --databases ${DATA_DIR}/db.fasta \
    --enzyme trypsinp \
    --min-length 6

OMP_NUM_THREADS=1 \
    python -u -m picked_group_fdr.pipeline.run_mokapot \
        0.01 0.01 ${RESULT_DIR} 1

python -u -m picked_group_fdr.pipeline.update_evidence_from_pout \
    --mq_evidence ${DATA_DIR}/evidence.txt \
    --perc_results ${RESULT_DIR}/andromeda.mokapot.psms.txt \
        ${RESULT_DIR}/andromeda.mokapot.decoys.psms.txt \
    --mq_evidence_out ${RESULT_DIR}/evidence.txt

OMP_DYNAMIC=FALSE OMP_NUM_THREADS=1 \
    python3 -u -m picked_group_fdr \
        --mq_evidence ${RESULT_DIR}/evidence.txt \
        --protein_groups_out ${RESULT_DIR}/proteinGroups.txt \
        --do_quant \
        --fasta ${DATA_DIR}/db.fasta \
        --methods picked_protein_group_mq_input \
        --num_threads 1 \
        --enzyme trypsinp \
        --min-length 6 \
        --special-aas KR

OMP_DYNAMIC=FALSE OMP_NUM_THREADS=1 \
    python3 -u -m picked_group_fdr \
        --mq_evidence ${RESULT_DIR}/evidence.txt \
        --protein_groups_out ${RESULT_DIR}/geneGroups.txt \
        --do_quant \
        --gene_level \
        --fasta ${DATA_DIR}/db.fasta \
        --methods picked_protein_group_mq_input \
        --num_threads 1 \
        --enzyme trypsinp \
        --min-length 6 \
        --special-aas KR \
        --suppress_missing_peptide_warning

python3 -u -m picked_group_fdr.pipeline.filter_fdr_maxquant \
    --mq_msms ${RESULT_DIR}/evidence.txt \
    --mq_msms_out ${RESULT_DIR}/evidence_fdr0.01.txt \
    --mq_protein_groups ${RESULT_DIR}/proteinGroups.txt \
    --mq_protein_groups_out ${RESULT_DIR}/proteinGroups_fdr0.01.txt \
    --fdr_cutoff 0.01 \
    --psm_level_fdr