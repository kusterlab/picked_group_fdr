DATA_DIR="data/lfq_example"
RESULT_DIR="tests/system_tests/test_pipeline"

# exit on first error
set -e

mkdir -p ${RESULT_DIR}

rm -f ${RESULT_DIR}/andromeda.{tab,mokapot.psms.txt,mokapot.decoy.psms.txt}
rm -f ${RESULT_DIR}/evidence.txt ${RESULT_DIR}/proteinGroups.txt

python -u -m picked_group_fdr.pipeline.andromeda2pin \
    ${DATA_DIR}/evidence.txt \
    --outputTab ${RESULT_DIR}/andromeda.tab \
    --databases ${DATA_DIR}/db.fasta \
    --enzyme trypsinp \
    --min-length 6
diff -q ${RESULT_DIR}/andromeda.tab ${RESULT_DIR}/andromeda.reference.tab

OMP_NUM_THREADS=1 \
    python -u -m picked_group_fdr.pipeline.run_mokapot \
        0.01 0.01 ${RESULT_DIR} 1
diff -q ${RESULT_DIR}/andromeda.mokapot.psms.txt ${RESULT_DIR}/andromeda.mokapot.psms.reference.txt
diff -q ${RESULT_DIR}/andromeda.mokapot.decoy.psms.txt ${RESULT_DIR}/andromeda.mokapot.decoy.psms.reference.txt

python -u -m picked_group_fdr.pipeline.update_evidence_from_pout \
    --mq_evidence ${DATA_DIR}/evidence.txt \
    --perc_results ${RESULT_DIR}/andromeda.mokapot.psms.txt \
        ${RESULT_DIR}/andromeda.mokapot.decoy.psms.txt \
    --mq_evidence_out ${RESULT_DIR}/evidence.txt
diff -q ${RESULT_DIR}/evidence.txt ${RESULT_DIR}/evidence.reference.txt

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
diff -q ${RESULT_DIR}/proteinGroups.txt ${RESULT_DIR}/proteinGroups.reference.txt

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
# TODO: reproducible order of proteins in gene groups
# diff -q ${RESULT_DIR}/geneGroups.txt ${RESULT_DIR}/geneGroups.reference.txt

python3 -u -m picked_group_fdr.pipeline.filter_fdr_maxquant \
    --mq_msms ${RESULT_DIR}/evidence.txt \
    --mq_msms_out ${RESULT_DIR}/evidence_fdr0.01.txt \
    --mq_protein_groups ${RESULT_DIR}/proteinGroups.txt \
    --mq_protein_groups_out ${RESULT_DIR}/proteinGroups_fdr0.01.txt \
    --fdr_cutoff 0.01 \
    --psm_level_fdr
diff -q ${RESULT_DIR}/evidence_fdr0.01.txt ${RESULT_DIR}/evidence_fdr0.01.reference.txt
diff -q ${RESULT_DIR}/proteinGroups_fdr0.01.txt ${RESULT_DIR}/proteinGroups_fdr0.01.reference.txt
