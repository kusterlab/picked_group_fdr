DATA_DIR="data/fragpipe_example"
RESULT_DIR="tests/system_tests/test_fragpipe"

# exit on first error
set -e

mkdir -p ${RESULT_DIR}

fasta_files="${DATA_DIR}/iprg2016_with_labels.fasta"

fragpipe_psm_files="${DATA_DIR}/fragpipe/A1/psm.tsv \
    ${DATA_DIR}/fragpipe/A2/psm.tsv \
    ${DATA_DIR}/fragpipe/A3/psm.tsv \
    ${DATA_DIR}/fragpipe/B1/psm.tsv \
    ${DATA_DIR}/fragpipe/B2/psm.tsv \
    ${DATA_DIR}/fragpipe/B3/psm.tsv \
    ${DATA_DIR}/fragpipe/C1/psm.tsv \
    ${DATA_DIR}/fragpipe/C2/psm.tsv \
    ${DATA_DIR}/fragpipe/C3/psm.tsv"

# Option 1: directly create combined_protein.tsv using combined_ion.tsv from a previous IonQuant run
python3 -u -m picked_group_fdr \
    --fasta ${fasta_files} \
    --fragpipe_psm ${fragpipe_psm_files} \
    --combined_ion ${DATA_DIR}/fragpipe/combined_ion.tsv \
    --protein_groups_out ${RESULT_DIR}/combined_protein.tsv \
    --do_quant \
    --lfq_min_peptide_ratios 1 \
    --methods fragpipe
diff ${RESULT_DIR}/combined_protein.tsv ${RESULT_DIR}/combined_protein.reference.tsv

# Option 2: create protein groups without quant
python3 -u -m picked_group_fdr \
    --fasta ${fasta_files} \
    --fragpipe_psm ${fragpipe_psm_files} \
    --protein_groups_out ${RESULT_DIR}/proteinGroups.txt \
    --methods fragpipe
diff ${RESULT_DIR}/proteinGroups.txt ${RESULT_DIR}/proteinGroups.reference.txt

# Option 3: use intermediate proteinGroups.txt (from Option 2) to create combined_protein.tsv using combined_ion.tsv from a previous IonQuant run
python3 -u -m picked_group_fdr.pipeline.update_fragpipe_results \
    --fasta ${fasta_files} \
    --fragpipe_psm ${fragpipe_psm_files} \
    --protein_groups ${RESULT_DIR}/proteinGroups.txt \
    --combined_ion ${DATA_DIR}/fragpipe/combined_ion.tsv \
    --output_folder ${RESULT_DIR}_from_combined_ion
diff ${RESULT_DIR}_from_combined_ion/combined_protein.tsv ${RESULT_DIR}_from_combined_ion/combined_protein.reference.tsv

# Option 4: update psm.tsv and protein.tsv using the intermediate proteinGroups.txt (from Option 2) to be requantified with IonQuant (see ./run_ionquant.sh)
python3 -u -m picked_group_fdr.pipeline.update_fragpipe_results \
    --fasta ${fasta_files} \
    --fragpipe_psm ${fragpipe_psm_files} \
    --protein_groups ${RESULT_DIR}/proteinGroups.txt \
    --output_folder ${RESULT_DIR}
diff ${RESULT_DIR}/A1/protein.tsv ${RESULT_DIR}/A1/protein.reference.tsv
diff ${RESULT_DIR}/A1/psm.tsv ${RESULT_DIR}/A1/psm.reference.tsv