import os
import toml

from .competition import ProteinCompetitionStrategyFactory
from .scoring import ProteinScoringStrategy
from .grouping import ProteinGroupingStrategyFactory
    

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
    
    if args.methods:
        for method in args.methods.split(","):
            configs.append(parse_method_toml(method))
        return configs
    
    # final method
    configs.append(parse_method_toml('picked_protein_group'))    
    
    ### MANUSCRIPT FIGURES ###
    
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


def parse_method_toml(method: str):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    method_toml_file = os.path.join(dir_path, 'methods', f'{method}.toml')
    method = toml.load(method_toml_file)
    
    picked_strategy = ProteinCompetitionStrategyFactory(method['pickedStrategy'])
    score_type = method['scoreType']
    if method['sharedPeptides'] == "razor":
        score_type += " razor"
    scoring_strategy = ProteinScoringStrategy(score_type)
    grouping_strategy = ProteinGroupingStrategyFactory(method['grouping'])
        
    return {'pickedStrategy': picked_strategy, 'scoreType': scoring_strategy, 'grouping': grouping_strategy, 'label': method['label'] }


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

