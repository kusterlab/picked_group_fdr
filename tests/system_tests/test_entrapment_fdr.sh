DATA_DIR="data/lfq_example"
RESULT_DIR="tests/system_tests/test_pipeline"

# exit on first error
set -e

python -m picked_group_fdr.pipeline.entrapment_fdr \
    --protein_groups_files ${RESULT_DIR}/proteinGroups.txt \
    --figure_base_fn ${RESULT_DIR}/protein_entrapment_fdr

python -m picked_group_fdr.pipeline.entrapment_fdr \
    --protein_col "Leading proteins" \
    --peptides_files ${RESULT_DIR}/evidence.txt \
    --figure_base_fn ${RESULT_DIR}/rescored_peptide_entrapment_fdr

python -m picked_group_fdr.pipeline.entrapment_fdr \
    --protein_col "Leading proteins" \
    --peptides_files ${DATA_DIR}/evidence.txt \
    --figure_base_fn ${RESULT_DIR}/andromeda_peptide_entrapment_fdr