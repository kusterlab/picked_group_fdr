WORK_DIR=./fragpipe
RESULT_DIR=./results

cd ${WORK_DIR}

mkdir -p ${RESULT_DIR}

fasta_files="../iprg2016_with_labels.fasta"

fragpipe_psm_files="./A1/psm.tsv ./A2/psm.tsv ./A3/psm.tsv \
    ./B1/psm.tsv ./B2/psm.tsv ./B3/psm.tsv \
    ./C1/psm.tsv ./C2/psm.tsv ./C3/psm.tsv"

# Option 1: directly create combined_protein.tsv using combined_ion.tsv from a previous IonQuant run
python3 -u -m picked_group_fdr \
    --fasta ${fasta_files} \
    --fragpipe_psm ${fragpipe_psm_files} \
    --combined_ion ./combined_ion.tsv \
    --protein_groups_out ${RESULT_DIR}/combined_protein.tsv \
    --do_quant \
    --lfq_min_peptide_ratios 1 \
    --methods fragpipe

# Option 2: create protein groups without quant
python3 -u -m picked_group_fdr \
    --fasta ${fasta_files} \
    --fragpipe_psm ${fragpipe_psm_files} \
    --protein_groups_out ${RESULT_DIR}/proteinGroups.txt \
    --methods fragpipe

# Option 3: use intermediate proteinGroups.txt (from Option 2) to create combined_protein.tsv using combined_ion.tsv from a previous IonQuant run
python3 -u -m picked_group_fdr.pipeline.update_fragpipe_results \
    --fasta ${fasta_files} \
    --fragpipe_psm ${fragpipe_psm_files} \
    --protein_groups ${RESULT_DIR}/proteinGroups.txt \
    --combined_ion ./combined_ion.tsv \
    --output_folder ${RESULT_DIR}_from_combined_ion

# Option 4: update psm.tsv and protein.tsv using the intermediate proteinGroups.txt (from Option 2) to be requantified with IonQuant (see ./run_ionquant.sh)
python3 -u -m picked_group_fdr.pipeline.update_fragpipe_results \
    --fasta ${fasta_files} \
    --fragpipe_psm ${fragpipe_psm_files} \
    --protein_groups ${RESULT_DIR}/proteinGroups.txt \
    --output_folder ${RESULT_DIR}

