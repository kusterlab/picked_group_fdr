DATA_DIR="data/lfq_example"
RESULT_DIR="tests/system_tests/test_digest"

# exit on first error
set -e

mkdir -p ${RESULT_DIR}

# generate prosit input file
python -m picked_group_fdr.digest \
    --prosit_input ${RESULT_DIR}/test_digest_prosit_file.txt \
    --fasta ${DATA_DIR}/db.fasta
diff ${RESULT_DIR}/test_digest_prosit_file.txt ${RESULT_DIR}/test_digest_prosit_file.reference.txt

# multiple digestion parameters
python -m picked_group_fdr.digest \
    --prosit_input ${RESULT_DIR}/test_digest_multiple_proteases_prosit_file.txt \
    --enzyme trypsin trypsinp \
    --cleavages 2 3 \
    --fasta ${DATA_DIR}/db.fasta
diff ${RESULT_DIR}/test_digest_multiple_proteases_prosit_file.txt ${RESULT_DIR}/test_digest_multiple_proteases_prosit_file.reference.txt

# peptide to protein map
python -m picked_group_fdr.digest \
    --peptide_protein_map ${RESULT_DIR}/test_digest_peptide_protein_map.txt \
    --enzyme trypsinp \
    --cleavages 2 \
    --min-length 6 \
    --fasta ${DATA_DIR}/db.fasta
diff ${RESULT_DIR}/test_digest_peptide_protein_map.txt ${RESULT_DIR}/test_digest_peptide_protein_map.reference.txt

# multiple digestion parameters
python -m picked_group_fdr.digest \
    --peptide_protein_map ${RESULT_DIR}/test_digest_multiple_proteases_peptide_protein_map.txt \
    --enzyme trypsin trypsinp \
    --cleavages 2 3 \
    --min-length 6 \
    --fasta ${DATA_DIR}/db.fasta
diff ${RESULT_DIR}/test_digest_multiple_proteases_peptide_protein_map.txt ${RESULT_DIR}/test_digest_multiple_proteases_peptide_protein_map.reference.txt

# ibaq peptide count
python -m picked_group_fdr.digest \
    --ibaq_map ${RESULT_DIR}/test_digest_multiple_proteases_ibaq_map.txt \
    --enzyme trypsin trypsinp \
    --cleavages 2 3 \
    --fasta ${DATA_DIR}/db.fasta
diff ${RESULT_DIR}/test_digest_multiple_proteases_ibaq_map.txt ${RESULT_DIR}/test_digest_multiple_proteases_ibaq_map.reference.txt
