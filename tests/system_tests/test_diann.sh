DATA_DIR="data/diann_example"
RESULT_DIR="tests/system_tests/test_diann"

# exit on first error
set -ex

mkdir -p ${RESULT_DIR}

fasta_files="${DATA_DIR}/diann/human_yeast_ecoli.fasta"

diann_report_file="${DATA_DIR}/diann/report.parquet"


# Option 1: directly create creport.pg_matrix.tsv (DIA-NN format)
python3 -u -m picked_group_fdr \
    --fasta ${fasta_files} \
    --fasta_use_uniprot_id \
    --diann_reports ${diann_report_file} \
    --protein_groups_out ${RESULT_DIR}/report.pg_matrix.tsv \
    --do_quant \
    --lfq_min_peptide_ratios 1 \
    --methods diann
diff -q ${RESULT_DIR}/report.pg_matrix.tsv ${RESULT_DIR}/report.pg_matrix.reference.tsv

# Option 2: create a combined_protein.tsv (FragPipe format)
python3 -u -m picked_group_fdr \
    --fasta ${fasta_files} \
    --fasta_use_uniprot_id \
    --diann_reports ${diann_report_file} \
    --protein_groups_out ${RESULT_DIR}/combined_protein.tsv \
    --do_quant \
    --lfq_min_peptide_ratios 1 \
    --output_format fragpipe \
    --methods diann
diff -q ${RESULT_DIR}/combined_protein.tsv ${RESULT_DIR}/combined_protein.reference.tsv

# Option 3: create a proteinGroups.txt (MaxQuant format)
python3 -u -m picked_group_fdr \
    --fasta ${fasta_files} \
    --fasta_use_uniprot_id \
    --diann_reports ${diann_report_file} \
    --protein_groups_out ${RESULT_DIR}/proteinGroups.txt \
    --do_quant \
    --lfq_min_peptide_ratios 1 \
    --output_format maxquant \
    --methods diann
diff -q ${RESULT_DIR}/proteinGroups.txt ${RESULT_DIR}/proteinGroups.reference.txt
