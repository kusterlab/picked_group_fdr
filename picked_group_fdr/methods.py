from .competition import PickedStrategy, PickedGroupStrategy, ClassicStrategy
from .scoring import ProteinScoringStrategy
from .grouping import RescuedSubsetGrouping, NoGrouping, MQNativeGrouping, SubsetGrouping, RescuedMQNativeGrouping


def get_methods(args):
    """
    pickedStrategy: PickedStrategy() = picked FDR
                    ClassicStrategy() = classic FDR
    
    grouping:       MQNativeGrouping() = MaxQuant grouping from proteinGroups.txt
                    SubsetGrouping() = Emulate protein grouping of MaxQuant based on evidence.txt
                                       (currently does not work with simulated datasets since peptideToProteinMap does not contain entrapment labels)
                    NoGrouping() = No protein grouping, each protein is in its own group
                    +Rescued = Rescue protein groups by only considering peptides below 1% protein FDR threshold
    """
    configs = list()
    
    ### WELL-CALIBRATED METHODS
    
    # final method
    configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("bestPEP"), 'grouping' : RescuedSubsetGrouping(), 'label' : 'Picked Protein Group FDR'})
    
    # picked FDR
    # N.B. Add "remap" to the scoreType if the fasta database used for protein grouping is different from the one used by Percolator
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("Perc remap bestPEP"), 'grouping' : RescuedSubsetGrouping()})
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("bestPEP"), 'grouping' : RescuedSubsetGrouping()})
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("Perc bestPEP"), 'grouping' : RescuedSubsetGrouping()})
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("Perc bestPEP"), 'grouping' : RescuedMQNativeGrouping()})
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("Perc multPEP"), 'grouping' : MQNativeGrouping()})
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("Perc multPEP"), 'grouping' : SubsetGrouping()})
    #configs.append({'pickedStrategy' : PickedStrategy(), 'scoreType' : ProteinScoringStrategy("bestPEP"), 'grouping' : NoGrouping()})
    
    # classic FDR
    #configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("Perc bestPEP"), 'grouping' : RescuedSubsetGrouping()})
    #configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("Perc bestPEP"), 'grouping' : RescuedMQNativeGrouping()})
    #configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("Andromeda"), 'grouping' : MQNativeGrouping()})
    #configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("bestPEP"), 'grouping' : SubsetGrouping()})
    
    # proteotypicity score experiments
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("bestPEP proteotypicity"), 'grouping' : RescuedSubsetGrouping()})
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("Perc bestPEP proteotypicity"), 'grouping' : SubsetGrouping()})
    
    
    ### EXISTING APPROACHES ###
    
    # MaxQuant scores
    #configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("MQ_protein with_shared", mq_protein_groups_file = args.mq_protein_groups), 'grouping' : MQNativeGrouping()})
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("MQ_protein with_shared", mq_protein_groups_file = args.mq_protein_groups), 'grouping' : MQNativeGrouping()})
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("Andromeda razor with_shared"), 'grouping' : MQNativeGrouping()})
    
    # MaxQuant emulator
    #configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("multPEP razor"), 'grouping' : SubsetGrouping()})
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("bestPEP"), 'grouping' : RescuedSubsetGrouping()})
    #configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("bestPEP razor"), 'grouping' : SubsetGrouping()})
    #configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("bestPEP razor"), 'grouping' : SubsetGrouping()})
    
    # MaxQuant emulator with Percolator PEP
    #configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("Perc bestPEP razor"), 'grouping' : SubsetGrouping()})
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("Perc bestPEP razor"), 'grouping' : SubsetGrouping()})
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("Perc multPEP razor"), 'grouping' : SubsetGrouping()})
    
    # Percolator 3 / PrDB approach: picked, best score, no shared peptides, no grouping
    #configs.append({'pickedStrategy' : PickedStrategy(), 'scoreType' : ProteinScoringStrategy("Perc bestPEP"), 'grouping' : NoGrouping()})
    #configs.append({'pickedStrategy' : PickedStrategy(), 'scoreType' : ProteinScoringStrategy("bestPEP"), 'grouping' : NoGrouping()})
    #configs.append({'pickedStrategy' : PickedStrategy(), 'scoreType' : ProteinScoringStrategy("Andromeda"), 'grouping' : NoGrouping()})
    
    
    ### MANUSCRIPT FIGURES ###
    
    #################
    ### Figure 3A ###
    #################
    
    # performance on Wang base with Percolator vs MaxQuant
    # for database in {prdb_uniprot,prdb_swissprot+iso,prdb_swissprot}; do python3 -um picked_group_fdr --mq_evidence ${DATA_DIR}/Dongxue_original_PrDB_runs/evidence_0.01peptideFDR_per_rawfile.txt --perc_evidence ${DATA_DIR}/Dongxue_original_PrDB_runs/peptides.all.txt --fasta ${DATA_DIR}/fasta/${database}.fasta --suppress_missing_peptide_warning --protein_groups_out ${DATA_DIR}/Dongxue_original_PrDB_runs/proteinGroups_${database}.txt | tee ${DATA_DIR}/Dongxue_original_PrDB_runs/proteinGroups_${database}.log; done
    
    #configs.append({'pickedStrategy' : PickedStrategy(), 'scoreType' : ProteinScoringStrategy("Perc remap bestPEP"), 'grouping' : NoGrouping(), 'label' : 'Savitski' })
    #configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("multPEP razor"), 'grouping' : SubsetGrouping(), 'label' : 'MaxQuant' })


    #################
    ### Figure 4B ###
    #################
    
    # calibration (s0.5)
    # python3 -m picked_group_fdr --mq_evidence ${DATA_DIR}/Dongxue_mimic_final_runs/mimic_s0.5_m4_peptFDR0.01_MQ/evidence_all_cleaned.txt --perc_evidence ${DATA_DIR}/Dongxue_mimic_final_runs/mimic_s0.5_m4_peptFDR0.01_MQ/peptides.all.txt --fasta ${DATA_DIR}/fasta/uniprot-proteome_human_190327.with_mimic_s0.5_m4.fasta --enzyme trypsinp --cleavages 0 --protein_groups_out ${DATA_DIR}/Dongxue_mimic_final_runs/mimic_s0.5_m4_peptFDR0.01_MQ/proteinGroups.txt --figure_base_fn results/221004_manuscript_figures/razor-vs-discard
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("Perc bestPEP"), 'grouping' : SubsetGrouping(), 'label' : 'Discard + Picked group' })
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("Perc bestPEP razor"), 'grouping' : SubsetGrouping(), 'label' : 'Razor + Picked group' })
    
    
    #################
    ### Figure 5B ###
    #################
    
    # number of shared peptides
    # for database in {prdb_uniprot,prdb_swissprot+iso,prdb_swissprot}; do python3 -m picked_group_fdr --perc_evidence ${DATA_DIR}/Dongxue_original_PrDB_runs/peptides_0.01peptideFDR_MQ_per_rawfile.all.txt --fasta ${DATA_DIR}/fasta/${database}.fasta --suppress_missing_peptide_warning; done
#    configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("Perc bestPEP remap"), 'grouping' : NoGrouping()})
#    configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("Perc bestPEP remap"), 'grouping' : SubsetGrouping()})
#    configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("Perc bestPEP remap"), 'grouping' : RescuedSubsetGrouping()})
    
    
    #################
    ### Figure 5C ###
    #################
    
    # calibration final
    
    # python3 -um picked_group_fdr --mq_evidence ${DATA_DIR}/Dongxue_mimic_final_runs/mimic_s0.5_m4_peptFDR0.01_MQ/evidence_all_cleaned.txt --perc_evidence ${DATA_DIR}/Dongxue_mimic_final_runs/mimic_s0.5_m4_peptFDR0.01_MQ/peptides.all.txt --fasta ${DATA_DIR}/fasta/uniprot-proteome_human_190327.with_mimic_s0.5_m4.fasta --enzyme trypsinp --cleavages 0 --figure_base_fn results/221004_manuscript_figures/protein --protein_groups_out ${DATA_DIR}/Dongxue_mimic_final_runs/mimic_s0.5_m4_peptFDR0.01_MQ/proteinGroups.txt | tee ${DATA_DIR}/Dongxue_mimic_final_runs/mimic_s0.5_m4_peptFDR0.01_MQ/proteinGroups.log
    
#    configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("multPEP razor"), 'grouping' : SubsetGrouping(), 'label' : 'MaxQuant' })
#    configs.append({'pickedStrategy' : PickedStrategy(), 'scoreType' : ProteinScoringStrategy("Perc bestPEP"), 'grouping' : NoGrouping(), 'label' : 'Savitski' })
#    configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("Perc bestPEP"), 'grouping' : RescuedSubsetGrouping(), 'label' : 'Classic Protein Group FDR' })
#    configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("Perc bestPEP"), 'grouping' : RescuedSubsetGrouping(), 'label' : 'Picked Protein Group FDR' })
    
    
#    configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("multPEP razor"), 'grouping' : SubsetGrouping(), 'label' : 'MaxQuant' })
#    configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("Perc remap multPEP razor"), 'grouping' : SubsetGrouping(), 'label' : 'MaxQuant Perc PEP' })
#    configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("multPEP razor"), 'grouping' : SubsetGrouping(), 'label' : 'MaxQuant picked' })
#    configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("Perc remap multPEP razor"), 'grouping' : SubsetGrouping(), 'label' : 'MaxQuant Perc PEP picked' })
#    configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("Perc remap bestPEP razor"), 'grouping' : SubsetGrouping(), 'label' : 'MaxQuant Best Perc PEP picked' })
#    configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("Perc remap bestPEP"), 'grouping' : SubsetGrouping(), 'label' : 'Classic Protein Group FDR no rescue' })
#    configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("Perc remap bestPEP"), 'grouping' : RescuedSubsetGrouping(), 'label' : 'Classic Protein Group FDR'})
#    configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("Perc remap bestPEP"), 'grouping' : SubsetGrouping(), 'label' : 'Picked Protein Group FDR no rescue' })
#    configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("Perc remap bestPEP"), 'grouping' : RescuedSubsetGrouping(), 'label' : 'Picked Protein Group FDR'})
    
    
    #################
    ### Figure 5D ###
    #################
    
    # performance on Wang base with Percolator vs MaxQuant
    # for database in {prdb_uniprot,prdb_swissprot+iso,prdb_swissprot}; do python3 -m picked_group_fdr --mq_evidence ${DATA_DIR}/Dongxue_original_PrDB_runs/evidence_0.01peptideFDR_MQ_per_rawfile.txt --perc_evidence ${DATA_DIR}/Dongxue_original_PrDB_runs/peptides_0.01peptideFDR_MQ_per_rawfile.all.txt --fasta ${DATA_DIR}/fasta/${database}.fasta --suppress_missing_peptide_warning --protein_groups_out ${DATA_DIR}/Dongxue_original_PrDB_runs/proteinGroups_${database}.txt | tee ${DATA_DIR}/Dongxue_original_PrDB_runs/proteinGroups_${database}.log; done
    
#    configs.append({'pickedStrategy' : PickedStrategy(), 'scoreType' : ProteinScoringStrategy("Perc remap bestPEP"), 'grouping' : NoGrouping(), 'label' : 'Savitski' })
#    configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("multPEP razor"), 'grouping' : SubsetGrouping(), 'label' : 'MaxQuant' })
#    configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("Perc remap bestPEP"), 'grouping' : RescuedSubsetGrouping(), 'label' : 'Classic Protein Group FDR' })
#    configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("Perc remap bestPEP"), 'grouping' : RescuedSubsetGrouping(), 'label' : 'Picked Protein Group FDR' })
    
    
    ################
    ### Figure 6 ###
    ################
    
    # entrapment simulation
    # for NUMEXP in {10,100,200,300,400,500}; do python3 -um picked_group_fdr --mq_evidence ${DATA_DIR}/simulations/evidence.${NUMEXP}_experiments_peptFDR0.01_no_recalibration.txt --peptide_protein_map ${DATA_DIR}/fasta/uniprot-proteome_human_190327.peptide_to_protein_map.long_protein_ids.txt --mq_protein_groups ${DATA_DIR}/simulations/proteinGroups.${NUMEXP}_experiments_peptFDR0.01_no_recalibration.txt --figure_base_fn results/221004_manuscript_figures/simulated-${NUMEXP}-experiments-peptFDR0.01-no-recalibration | tee results/221004_manuscript_figures/simulated-${NUMEXP}-experiments-peptFDR0.01-no-recalibration.log; done
    #configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("multPEP razor"), 'grouping' : MQNativeGrouping(), 'label' : 'MaxQuant'})
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("bestPEP razor"), 'grouping' : MQNativeGrouping(), 'label' : 'Razor + Picked group'})
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("bestPEP"), 'grouping' : MQNativeGrouping(), 'label' : 'Discard + Picked group'})
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("bestPEP"), 'grouping' : RescuedMQNativeGrouping(), 'label' : 'Picked Protein Group FDR'})
    
    #configs.append({'pickedStrategy' : PickedGroupStrategy("all"), 'scoreType' : ProteinScoringStrategy("Perc remap bestPEP"), 'grouping' : RescuedSubsetGrouping(), 'label' : 'Picked (all) Protein Group FDR'})
    #configs.append({'pickedStrategy' : PickedGroupStrategy("majority"), 'scoreType' : ProteinScoringStrategy("Perc remap bestPEP"), 'grouping' : RescuedSubsetGrouping(), 'label' : 'Picked (majority) Protein Group FDR'})
    #configs.append({'pickedStrategy' : PickedGroupStrategy("leading"), 'scoreType' : ProteinScoringStrategy("Perc remap bestPEP"), 'grouping' : RescuedSubsetGrouping(), 'label' : 'Picked (leading) Protein Group FDR'})
    
    
    #################
    ### Figure 7A ###
    #################
    
    # ProteomicsDB protein groups
    # for database in {swissprot,swissprot+iso,uniprot}; do python3 -um picked_group_fdr --perc_evidence ${DATA_DIR}/ProteomicsDB/percolator.${database}.all.peptides.txt --protein_groups_out ${DATA_DIR}/ProteomicsDB/proteinGroups_${database}.txt --fasta ${DATA_DIR}/fasta/prdb_${database}.fasta --suppress_missing_peptide_warning 2>&1 | tee ${DATA_DIR}/ProteomicsDB/proteinGroups_${database}.log; done
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("Perc bestPEP"), 'grouping' : RescuedSubsetGrouping(), 'label' : 'Picked Protein Group FDR'})
    #configs.append({'pickedStrategy' : PickedStrategy(), 'scoreType' : ProteinScoringStrategy("Perc bestPEP"), 'grouping' : NoGrouping(), 'label' : 'Savitski'})
    
    
    ###############################
    ### Supplementary Figure 1C ###
    ###############################
    
    # MQ PEP is not well calibrated for Savitski
    # python3 -um picked_group_fdr --mq_evidence ${DATA_DIR}/Dongxue_mimic_final_runs/mimic_s0.5_m4_peptFDR0.01_MQ/evidence_all_cleaned.txt --perc_evidence ${DATA_DIR}/Dongxue_mimic_final_runs/mimic_s0.5_m4_peptFDR0.01_MQ/peptides.all.txt --fasta ${DATA_DIR}/fasta/uniprot-proteome_human_190327.with_mimic_s0.5_m4.fasta --enzyme trypsinp --cleavages 0 --protein_groups_out results/221004_manuscript_figures/proteinGroups_protein-mq-pep --figure_base_fn results/221004_manuscript_figures/Supplementary_Figure1c | tee results/221004_manuscript_figures/protein-mq-pep.log
    #configs.append({'pickedStrategy' : PickedStrategy(), 'scoreType' : ProteinScoringStrategy("multPEP"), 'grouping' : NoGrouping(), 'label' : 'Savitski + mult MQ PEP' })
    #configs.append({'pickedStrategy' : PickedStrategy(), 'scoreType' : ProteinScoringStrategy("bestPEP"), 'grouping' : NoGrouping(), 'label' : 'Savitski + best MQ PEP' })
    #configs.append({'pickedStrategy' : PickedStrategy(), 'scoreType' : ProteinScoringStrategy("Perc bestPEP"), 'grouping' : NoGrouping(), 'label' : 'Savitski' })
    #configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("multPEP razor"), 'grouping' : SubsetGrouping(), 'label' : 'MaxQuant' })
    
    ##############################
    ### Supplementary Figure 2 ###
    ##############################
    
    # comparing MQ, Savitski MQ PEP and Savitski Perc PEP
    # 1% peptide-level FDR
    # for database in {prdb_uniprot,prdb_swissprot+iso,prdb_swissprot}; do python3 -um picked_group_fdr --mq_evidence ${DATA_DIR}/Dongxue_original_PrDB_runs/evidence_0.01peptideFDR_MQ_per_rawfile.txt --perc_evidence ${DATA_DIR}/Dongxue_original_PrDB_runs/peptides_0.01peptideFDR_MQ_per_rawfile.all.txt --fasta ${DATA_DIR}/fasta/${database}.fasta --suppress_missing_peptide_warning --protein_groups_out results/221004_manuscript_figures/proteinGroups_method_comparison_${database}.txt | tee results/221004_manuscript_figures/proteinGroups_method_comparison_${database}.log; done
    # 100% peptide-level FDR
    # for database in {prdb_uniprot,prdb_swissprot+iso,prdb_swissprot}; do python3 -um picked_group_fdr --mq_evidence ${DATA_DIR}/Dongxue_original_PrDB_runs/evidence.txt --perc_evidence ${DATA_DIR}/Dongxue_original_PrDB_runs/peptides.all.txt --fasta ${DATA_DIR}/fasta/${database}.fasta --suppress_missing_peptide_warning --protein_groups_out results/221004_manuscript_figures/proteinGroups_method_comparison_100FDR_${database}.txt | tee results/221004_manuscript_figures/proteinGroups_method_comparison_100FDR_${database}.log; done
    
    #configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("multPEP razor"), 'grouping' : SubsetGrouping(), 'label' : 'MaxQuant' })
    #configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("Perc remap bestPEP"), 'grouping' : NoGrouping(), 'label' : 'Savitski + Classic FDR' })
    #configs.append({'pickedStrategy' : PickedStrategy(), 'scoreType' : ProteinScoringStrategy("Perc remap bestPEP"), 'grouping' : NoGrouping(), 'label' : 'Savitski' })
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("Perc remap bestPEP razor"), 'grouping' : SubsetGrouping(), 'label' : 'Razor + Picked group'})
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("Perc remap bestPEP"), 'grouping' : SubsetGrouping(), 'label' : 'Discard + Picked group'})
    #configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("Perc remap bestPEP"), 'grouping' : RescuedSubsetGrouping(), 'label' : 'Classic Protein Group FDR'})
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("Perc remap bestPEP"), 'grouping' : RescuedSubsetGrouping(), 'label' : 'Picked Protein Group FDR'})
    
    ##############################
    ### Supplementary Figure 3 ###
    ##############################
    
    # comparing MQ, Savitski MQ PEP and Savitski Perc PEP
    # 1% peptide-level FDR
    # for database in {prdb_uniprot,prdb_swissprot+iso,prdb_swissprot}; do python3 -um picked_group_fdr --mq_evidence ${DATA_DIR}/Dongxue_original_PrDB_runs/evidence_0.01peptideFDR_MQ_per_rawfile.txt --perc_evidence ${DATA_DIR}/Dongxue_original_PrDB_runs/peptides_0.01peptideFDR_MQ_per_rawfile.all.txt --fasta ${DATA_DIR}/fasta/${database}.fasta --suppress_missing_peptide_warning; done
    
    #configs.append({'pickedStrategy' : PickedStrategy(), 'scoreType' : ProteinScoringStrategy("Perc remap bestPEP"), 'grouping' : NoGrouping(), 'label' : 'Savitski' })
    
    ##############################
    ### Supplementary Figure 4 ###
    ##############################
    
    # supplement: razor never works
    # python3 -um picked_group_fdr --mq_evidence ${DATA_DIR}/Dongxue_mimic_final_runs/mimic_s0.5_m4_peptFDR0.01_MQ/evidence_all_cleaned.txt --perc_evidence ${DATA_DIR}/Dongxue_mimic_final_runs/mimic_s0.5_m4_peptFDR0.01_MQ/peptides.all.txt --fasta ${DATA_DIR}/fasta/uniprot-proteome_human_190327.with_mimic_s0.5_m4.fasta --enzyme trypsinp --cleavages 0 --figure_base_fn results/221004_manuscript_figures/protein-razor --protein_groups_out results/221004_manuscript_figures/proteinGroups_protein-razor.txt | tee results/221004_manuscript_figures/proteinGroups_protein-razor.log
    
    #configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("multPEP razor"), 'grouping' : SubsetGrouping(), 'label' : 'MaxQuant'})
    #configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("bestPEP razor"), 'grouping' : SubsetGrouping(), 'label' : 'MaxQuant + best MQ PEP'})
    #configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("Perc bestPEP razor"), 'grouping' : SubsetGrouping(), 'label' : 'MaxQuant + best Perc PEP'})
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("multPEP razor"), 'grouping' : SubsetGrouping(), 'label' : 'MaxQuant + Picked group'})
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("bestPEP razor"), 'grouping' : SubsetGrouping(), 'label' : 'MaxQuant + best MQ PEP + Picked group'})
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("Perc bestPEP razor"), 'grouping' : SubsetGrouping(), 'label' : 'MaxQuant + best Perc PEP + Picked group'})

    ##############################
    ### Supplementary Figure 5 ###
    ##############################
    
    # supplement: calibration s0.5
    # python3 -um picked_group_fdr --mq_evidence ${DATA_DIR}/Dongxue_mimic_final_runs/mimic_s0.5_m4_peptFDR0.01_MQ/evidence_all_cleaned.txt --perc_evidence ${DATA_DIR}/Dongxue_mimic_final_runs/mimic_s0.5_m4_peptFDR0.01_MQ/peptides.all.txt --fasta ${DATA_DIR}/fasta/uniprot-proteome_human_190327.with_mimic_s0.5_m4.fasta --enzyme trypsinp --cleavages 0 --figure_base_fn results/221004_manuscript_figures/protein-s0-5-with-alts --protein_groups_out results/221004_manuscript_figures/proteinGroups_protein-s0-5-with-alts.txt | tee results/221004_manuscript_figures/proteinGroups_protein-s0-5-with-alts.log
    
    # supplement: calibration s0.04
    # python3 -um picked_group_fdr --mq_evidence ${DATA_DIR}/Dongxue_mimic_final_runs/mimic_s0.04_m4_peptFDR0.01_MQ/evidence_all_cleaned.txt --perc_evidence ${DATA_DIR}/Dongxue_mimic_final_runs/mimic_s0.04_m4_peptFDR0.01_MQ/peptides.all.txt --fasta ${DATA_DIR}/fasta/uniprot-proteome_human_190327.with_mimic_s0.04_m4.fasta --enzyme trypsinp --cleavages 0 --figure_base_fn results/221004_manuscript_figures/protein-s0-04-with-alts --protein_groups_out results/221004_manuscript_figures/proteinGroups_protein-s0-04-with-alts.txt | tee results/221004_manuscript_figures/proteinGroups_protein-s0-04-with-alts.log
    
    #configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("multPEP razor"), 'grouping' : SubsetGrouping(), 'label' : 'MaxQuant' })
    #configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("Perc remap bestPEP"), 'grouping' : NoGrouping(), 'label' : 'Savitski + Classic FDR' })
    #configs.append({'pickedStrategy' : PickedStrategy(), 'scoreType' : ProteinScoringStrategy("Perc remap bestPEP"), 'grouping' : NoGrouping(), 'label' : 'Savitski' })
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("Perc remap bestPEP razor"), 'grouping' : SubsetGrouping(), 'label' : 'Razor + Picked group'})
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("Perc remap bestPEP"), 'grouping' : SubsetGrouping(), 'label' : 'Discard + Picked group'})
    #configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("Perc remap bestPEP"), 'grouping' : RescuedSubsetGrouping(), 'label' : 'Classic Protein Group FDR'})
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("Perc remap bestPEP"), 'grouping' : RescuedSubsetGrouping(), 'label' : 'Picked Protein Group FDR'})

    
    #################################
    ### Supplementary Figure 6A-D ###
    #################################
    
    # supplement: entrapment simulation verification - simulation
    # for fdr in {1,01}; do for ratio in {5,04}; do python3 -um picked_group_fdr --mq_evidence ${DATA_DIR}/simulations/evidence.mimic_s0.${ratio}_m4_peptFDR0.${fdr}.30_experiments.txt --peptide_protein_map ${DATA_DIR}/fasta/uniprot-proteome_human_190327.with_mimic_s0.${ratio}_m4.missing_last_aa.peptide_to_protein_map.txt --mq_protein_groups ${DATA_DIR}/simulations/proteinGroups.mimic_s0.${ratio}_m4_peptFDR0.${fdr}.30_experiments.txt --figure_base_fn results/221004_manuscript_figures/simulated-verification-s0-${ratio}-fdr0-${fdr} --protein_groups_out results/221004_manuscript_figures/proteinGroups_simulated-verification-s0-${ratio}-fdr0-${fdr}.txt | tee results/221004_manuscript_figures/proteinGroups_simulated-verification-s0-${ratio}-fdr0-${fdr}.log; done; done
    
    #configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("multPEP razor"), 'grouping' : SubsetGrouping(), 'label' : 'MaxQuant' })
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("bestPEP razor"), 'grouping' : SubsetGrouping(), 'label' : 'Razor + Picked group'})
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("bestPEP"), 'grouping' : SubsetGrouping(), 'label' : 'Discard + Picked group'})
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("bestPEP"), 'grouping' : RescuedSubsetGrouping(), 'label' : 'Picked Protein Group FDR' })


    #################################
    ### Supplementary Figure 6E-H ###
    #################################
    
    # supplement: entrapment simulation verification - entrapment
    # for fdr in {1,01_MQ}; do for ratio in {5,04}; do python3 -um picked_group_fdr --mq_evidence ${DATA_DIR}/Dongxue_mimic_final_runs/mimic_s0.${ratio}_m4_peptFDR0.${fdr}/evidence_all_cleaned.txt --perc_evidence ${DATA_DIR}/Dongxue_mimic_final_runs/mimic_s0.${ratio}_m4_peptFDR0.${fdr}/peptides.all.txt --fasta ${DATA_DIR}/fasta/uniprot-proteome_human_190327.with_mimic_s0.${ratio}_m4.fasta --enzyme trypsinp --cleavages 0 --figure_base_fn results/221004_manuscript_figures/simulation-entrapment-verification-s0-${ratio}-fdr0-${fdr} --protein_groups_out results/221004_manuscript_figures/proteinGroups_simulation-entrapment-verification-s0-${ratio}-fdr0-${fdr}.txt | tee results/221004_manuscript_figures/proteinGroups_simulation-entrapment-verification-s0-${ratio}-fdr0-${fdr}.log; done; done
    
    #configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("multPEP razor"), 'grouping' : SubsetGrouping(), 'label' : 'MaxQuant' })
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("Perc bestPEP razor"), 'grouping' : SubsetGrouping(), 'label' : 'Razor + Picked group'})
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("Perc bestPEP"), 'grouping' : SubsetGrouping(), 'label' : 'Discard + Picked group'})
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("Perc bestPEP"), 'grouping' : RescuedSubsetGrouping(), 'label' : 'Picked Protein Group FDR' })
    

    ###################################
    ### Supplementary Figures 7,8,9 ###
    ###################################
    
    # supplementary figure 7: entrapment simulation without 1% per experiment peptide-FDR
    # for NUMEXP in {100,200,500}; do python3 -um picked_group_fdr --mq_evidence ${DATA_DIR}/simulations/evidence.${NUMEXP}_experiments_peptFDR0.01_no_recalibration.txt --peptide_protein_map ${DATA_DIR}/fasta/uniprot-proteome_human_190327.peptide_to_protein_map.long_protein_ids.txt --mq_protein_groups ${DATA_DIR}/simulations/proteinGroups.${NUMEXP}_experiments_peptFDR0.01_no_recalibration.txt --figure_base_fn results/221004_manuscript_figures/simulated-${NUMEXP}-experiments-peptFDR0.01-no-recalibration --protein_groups_out results/221004_manuscript_figures/proteinGroups_simulated-${NUMEXP}-experiments-peptFDR0.01-no-recalibration.txt | tee results/221004_manuscript_figures/proteinGroups_simulated-${NUMEXP}-experiments-peptFDR0.01-no-recalibration.log; done    
    
    # supplementary figure 8: entrapment simulation without 10% per experiment peptide-FDR
    # for NUMEXP in {100,200,500}; do python3 -um picked_group_fdr --mq_evidence ${DATA_DIR}/simulations/100-200-500-experiments/evidence.${NUMEXP}_experiments_peptFDR0.1_no_recalibration.txt --peptide_protein_map ${DATA_DIR}/fasta/uniprot-proteome_human_190327.peptide_to_protein_map.txt --mq_protein_groups ${DATA_DIR}/simulations/100-200-500-experiments/proteinGroups.${NUMEXP}_experiments_peptFDR0.1_no_recalibration.txt --figure_base_fn results/221004_manuscript_figures/simulated-${NUMEXP}-experiments-peptFDR0.1-no-recalibration --protein_groups_out results/221004_manuscript_figures/proteinGroups_simulated-${NUMEXP}-experiments-peptFDR0.1-no-recalibration.txt | tee results/221004_manuscript_figures/proteinGroups_simulated-${NUMEXP}-experiments-peptFDR0.1-no-recalibration.log; done
    
    # supplementary figure 9: entrapment simulation with 1% global peptide-FDR
    # for NUMEXP in {100,200,500}; do python3 -um picked_group_fdr --mq_evidence ${DATA_DIR}/simulations/100-200-500-experiments/evidence.${NUMEXP}_experiments_peptFDR0.01_with_recalibration.txt --peptide_protein_map ${DATA_DIR}/fasta/uniprot-proteome_human_190327.peptide_to_protein_map.txt --mq_protein_groups ${DATA_DIR}/simulations/100-200-500-experiments/proteinGroups.${NUMEXP}_experiments_peptFDR0.01_with_recalibration.txt --figure_base_fn results/221004_manuscript_figures/simulated-${NUMEXP}-experiments-peptFDR0.01-with-recalibration --protein_groups_out results/221004_manuscript_figures/proteinGroups_simulated-${NUMEXP}-experiments-peptFDR0.01-with-recalibration.txt | tee results/221004_manuscript_figures/proteinGroups_simulated-${NUMEXP}-experiments-peptFDR0.01-with-recalibration.log; done
    
    #configs.append({'pickedStrategy' : ClassicStrategy(), 'scoreType' : ProteinScoringStrategy("multPEP razor"), 'grouping' : MQNativeGrouping(), 'label' : 'MaxQuant'})
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("bestPEP razor"), 'grouping' : MQNativeGrouping(), 'label' : 'Razor + Picked group'})
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("bestPEP"), 'grouping' : MQNativeGrouping(), 'label' : 'Discard + Picked group'})
    #configs.append({'pickedStrategy' : PickedGroupStrategy(), 'scoreType' : ProteinScoringStrategy("bestPEP"), 'grouping' : RescuedMQNativeGrouping(), 'label' : 'Picked Protein Group FDR' })

    return configs


def short_description(scoreType, groupingStrategy, pickedStrategy, rescue_step, sep='_'):
    scoreLabel = scoreType.short_description()
    groupingLabel = groupingStrategy.short_description(rescue_step = rescue_step)
    razorLabel = scoreType.short_description_razor()
    fdrLabel = pickedStrategy.short_description()
    return f"{scoreLabel}{sep}{groupingLabel}{sep}{razorLabel}{sep}{fdrLabel}"


def long_description(scoreType, groupingStrategy, pickedStrategy, rescue_step, sep=', '):
    scoreLabel = scoreType.long_description()
    groupingLabel = groupingStrategy.long_description(rescue_step = rescue_step)
    razorLabel = scoreType.long_description_razor()
    fdrLabel = pickedStrategy.long_description()
    return f"{scoreLabel}{sep}{groupingLabel}{sep}{razorLabel}{sep}{fdrLabel}"

