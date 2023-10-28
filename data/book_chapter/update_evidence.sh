input_folders=(
    "./example_data/trypsin" 
    "./example_data/chymotrypsin" 
    "./example_data/LysC" 
    "./example_data/LysN"
)

for input_folder in "${input_folders[@]}"
do
    mkdir -p ${input_folder}/picked_group_fdr

    python3 -u -m picked_group_fdr.pipeline.update_evidence_from_pout \
        --mq_evidence ${input_folder}/maxquant/combined/txt/evidence.txt \
        --perc_results ${input_folder}/oktoberfest/results/percolator/rescore.percolator.psms.txt \
            ${input_folder}/oktoberfest/results/percolator/rescore.percolator.decoy.psms.txt \
        --mq_evidence_out ${input_folder}/picked_group_fdr/evidence.txt \
        --pout_input_type prosit
done