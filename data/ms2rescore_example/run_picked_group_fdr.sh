RESULT_DIR=./results

mkdir -p ${RESULT_DIR}

fasta_files="./sage/iprg2016_with_labels.fasta"

sage_ms2rescore_results_file="./sage/results.sage.ms2rescore.mokapot.psms.txt ./sage/results.sage.ms2rescore.mokapot.decoy.psms.txt"

# Option 1: directly create combined_protein.tsv (FragPipe format) using lfq.tsv
python3 -u -m picked_group_fdr \
    --fasta ${fasta_files} \
    --perc_evidence ${sage_ms2rescore_results_file} \
    --protein_groups_out ${RESULT_DIR}/combined_protein.tsv \
    --output_format fragpipe \
    --methods picked_protein_group_no_remap