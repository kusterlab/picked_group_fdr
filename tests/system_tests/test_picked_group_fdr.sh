DATA_DIR="data/lfq_example"
RESULT_DIR="tests/system_tests/test_picked_group_fdr"

mkdir -p ${RESULT_DIR}

# MaxQuant input
python3 -um picked_group_fdr --mq_evidence ${DATA_DIR}/evidence.txt \
    --fasta ${DATA_DIR}/db.fasta \
    --enzyme trypsinp \
    --min-length 6 \
    --protein_groups_out ${RESULT_DIR}/proteinGroups.txt \
    --method picked_protein_group_mq_input \
    --do_quant
diff ${RESULT_DIR}/proteinGroups.txt ${RESULT_DIR}/proteinGroups.reference.txt

# Percolator input
python3 -um picked_group_fdr \
    --perc_evidence ${DATA_DIR}/percolator_target.psms ${DATA_DIR}/percolator_decoy.psms \
    --method picked_protein_group_no_remap \
    --protein_groups_out ${RESULT_DIR}/proteinGroups_percolator.txt
diff ${RESULT_DIR}/proteinGroups_percolator.txt ${RESULT_DIR}/proteinGroups_percolator.reference.txt

# peptide to protein map input
python3 -um picked_group_fdr --mq_evidence ${DATA_DIR}/evidence.txt \
    --peptide_protein_map ${RESULT_DIR}/../test_digest/test_digest_peptide_protein_map.txt \
    --protein_groups_out ${RESULT_DIR}/proteinGroups_peptide_protein_map.txt \
    --method picked_protein_group_mq_input \
    --do_quant
diff ${RESULT_DIR}/proteinGroups_peptide_protein_map.txt ${RESULT_DIR}/proteinGroups_peptide_protein_map.reference.txt

# multiple fasta

# multiple proteases