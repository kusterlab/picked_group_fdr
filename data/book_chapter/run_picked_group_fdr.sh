RESULT_DIR=./example_data/command_line_results

mkdir -p ${RESULT_DIR}

input_folders=(
    "./example_data/trypsin" 
    "./example_data/chymotrypsin" 
    "./example_data/LysC" 
    "./example_data/LysN"
)

enzymes="trypsinp chymotrypsin+ lys-cp lys-n"
special_aas="KR none K K"
miscleavages="2 4 2 2"

# enzymes="trypsinp lys-cp"
# special_aas="KR K"
# miscleavages="2 2"

# enzymes="trypsinp"
# special_aas="KR"
# miscleavages=2

# enzymes="chymotrypsin+"
# special_aas="none"
# miscleavages=4

# enzymes="lys-cp"
# special_aas="K"
# miscleavages=2

# enzymes="lys-n"
# special_aas="K"
# miscleavages=2

fasta_files="./example_data/fasta/Homo_sapiens.GRCh38.pep.abinitio.fa \
    ./example_data/fasta/Homo_sapiens.GRCh38.pep.all.fa \
    ./example_data/fasta/UP000005640_9606.fasta \
    ./example_data/fasta/UP000005640_9606_additional.fasta"

mq_evidence_files=""
rescored_mq_evidence_files=""
for input_folder in "${input_folders[@]}"
do
    mq_evidence_files+="${input_folder}/maxquant/combined/txt/evidence.txt "
    rescored_mq_evidence_files+="${input_folder}/picked_group_fdr/evidence.txt "
done

# without_rescoring
python3 -u -m picked_group_fdr \
    --fasta ${fasta_files} \
    --mq_evidence ${mq_evidence_files} \
    --enzyme ${enzymes} \
    --special-aas ${special_aas} \
    --cleavages ${miscleavages} \
    --protein_groups_out ${RESULT_DIR}/original.proteinGroups.txt \
    --methods picked_protein_group_mq_input

# with rescoring
python3 -u -m picked_group_fdr \
    --fasta ${fasta_files} \
    --mq_evidence ${rescored_mq_evidence_files} \
    --enzyme ${enzymes} \
    --special-aas ${special_aas} \
    --cleavages ${miscleavages} \
    --protein_groups_out ${RESULT_DIR}/rescore.proteinGroups.txt \
    --methods picked_protein_group_mq_input

# with rescoring and quant
python3 -u -m picked_group_fdr \
    --fasta ${fasta_files} \
    --mq_evidence ${rescored_mq_evidence_files} \
    --enzyme ${enzymes} \
    --special-aas ${special_aas} \
    --cleavages ${miscleavages} \
    --protein_groups_out ${RESULT_DIR}/rescore.proteinGroups.with_quant.txt \
    --methods picked_protein_group_mq_input \
    --lfq_min_peptide_ratios 1 \
    --do_quant

python3 -u -m picked_group_fdr.pipeline.filter_fdr_maxquant \
    --mq_protein_groups ${RESULT_DIR}/original.proteinGroups.txt \
    --mq_protein_groups_out ${RESULT_DIR}/original.proteinGroups.fdr1.txt \
    --fdr_cutoff 0.01

python3 -u -m picked_group_fdr.pipeline.filter_fdr_maxquant \
    --mq_protein_groups ${RESULT_DIR}/rescore.proteinGroups.txt \
    --mq_protein_groups_out ${RESULT_DIR}/rescore.proteinGroups.fdr1.txt \
    --fdr_cutoff 0.01

python3 -u -m picked_group_fdr.pipeline.filter_fdr_maxquant \
    --mq_protein_groups ${RESULT_DIR}/rescore.proteinGroups.with_quant.txt \
    --mq_protein_groups_out ${RESULT_DIR}/rescore.proteinGroups.with_quant.fdr1.txt \
    --fdr_cutoff 0.01