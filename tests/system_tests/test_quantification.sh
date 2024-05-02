DATA_DIR="data/lfq_example"
RESULT_DIR="tests/system_tests/test_quantification"

# exit on first error
set -e

mkdir -p ${RESULT_DIR}

# MaxQuant input
python3 -um picked_group_fdr --mq_evidence ${DATA_DIR}/evidence.txt \
    --fasta ${DATA_DIR}/db.fasta \
    --enzyme trypsinp \
    --min-length 6 \
    --protein_groups_out ${RESULT_DIR}/proteinGroups.txt \
    --method picked_protein_group_mq_input
diff ${RESULT_DIR}/proteinGroups.txt ${RESULT_DIR}/proteinGroups.reference.txt

python3 -um picked_group_fdr.quantification --mq_evidence ${DATA_DIR}/evidence.txt \
    --mq_protein_groups ${RESULT_DIR}/proteinGroups.txt \
    --fasta ${DATA_DIR}/db.fasta \
    --enzyme trypsinp \
    --min-length 6 \
    --protein_groups_out ${RESULT_DIR}/proteinGroups_with_quant.txt
diff ${RESULT_DIR}/proteinGroups_with_quant.txt ${RESULT_DIR}/proteinGroups_with_quant.reference.txt


# peptide to protein map input
python3 -um picked_group_fdr --mq_evidence ${DATA_DIR}/evidence.txt \
    --peptide_protein_map ${RESULT_DIR}/../test_digest/test_digest_peptide_protein_map.txt \
    --protein_groups_out ${RESULT_DIR}/proteinGroups_peptide_protein_map.txt \
    --method picked_protein_group_mq_input
diff ${RESULT_DIR}/proteinGroups_peptide_protein_map.txt ${RESULT_DIR}/proteinGroups_peptide_protein_map.reference.txt

python3 -um picked_group_fdr.quantification --mq_evidence ${DATA_DIR}/evidence.txt \
    --mq_protein_groups ${RESULT_DIR}/proteinGroups_peptide_protein_map.txt \
    --peptide_protein_map ${RESULT_DIR}/../test_digest/test_digest_peptide_protein_map.txt \
    --protein_groups_out ${RESULT_DIR}/proteinGroups_peptide_protein_map_with_quant.txt
diff ${RESULT_DIR}/proteinGroups_peptide_protein_map_with_quant.txt ${RESULT_DIR}/proteinGroups_peptide_protein_map_with_quant.reference.txt