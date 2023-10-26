import logging
import os

import mokapot
import numpy as np
from joblib import parallel_backend

import picked_group_fdr.pipeline.andromeda2pin as andromeda2pin
import picked_group_fdr.pipeline.update_evidence_from_pout as update_evidence
import picked_group_fdr.pipeline.merge_pout as merge_pout
import picked_group_fdr.picked_group_fdr as picked_group_fdr


logger = logging.getLogger()

def run_picked_group_fdr_all(evidence_files, pout_files, fasta_file, output_dir, digest_params, input_type):
    try:
        if len(output_dir) == 0:
            raise RuntimeError("Please specify an output folder")
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        
        if input_type == "rescoring":
            run_update_evidence_rescoring(evidence_files, pout_files, output_dir)
            run_picked_group_fdr(output_dir, fasta_file, digest_params)
        elif input_type == "percolator":
            run_merge_pout(pout_files, fasta_file, output_dir, digest_params)
            run_picked_group_fdr_percolator_input(output_dir, fasta_file)
        else:
            run_andromeda_to_pin(evidence_files, fasta_file, output_dir, digest_params)
            run_mokapot(output_dir)
            run_update_evidence(evidence_files, output_dir)
            run_picked_group_fdr(output_dir, fasta_file, digest_params)
        
    except SystemExit as e:
        logger.error(f"Error while running Picked Group FDR, exited with error code: {e}.")
    except Exception as e:
        logger.error(f"Error while running Picked Group FDR: {e}")


def run_andromeda_to_pin(evidence_files, fasta_file, output_dir, digest_params):
    andromeda2pin.main(evidence_files + ['--outputTab', f'{output_dir}/pin.tab', '--databases', fasta_file] + digest_params)
    

def run_mokapot(output_dir):
    np.random.seed(0) # TODO: Make seed configurable
    psms = mokapot.read_pin(f'{output_dir}/pin.tab')
    with parallel_backend("threading"):
        results, models = mokapot.brew(psms)
    results.to_txt(dest_dir=output_dir, decoys = True)
    

def run_update_evidence(evidence_files, output_dir):
    update_evidence.main(
        ['--mq_evidence'] + evidence_files + 
        ['--perc_results', f'{output_dir}/mokapot.psms.txt', f'{output_dir}/mokapot.decoys.psms.txt', 
         '--mq_evidence_out', f'{output_dir}/evidence_percolator.txt'])


def run_update_evidence_rescoring(evidence_files, pout_files, output_dir):
    update_evidence.main(
        ['--mq_evidence'] + evidence_files + 
        ['--perc_results'] + pout_files +
        ['--mq_evidence_out', f'{output_dir}/evidence_percolator.txt',
         '--pout_input_type', 'prosit'])


def run_picked_group_fdr(output_dir, fasta_file, digest_params):
    picked_group_fdr.main(
        ['--mq_evidence', f'{output_dir}/evidence_percolator.txt',
         '--methods', 'picked_protein_group_mq_input', 
         '--do_quant',
         '--protein_groups_out', f'{output_dir}/proteinGroups_percolator.txt',
         '--fasta', fasta_file] + digest_params)


def run_merge_pout(pout_files, fasta_file, output_dir, digest_params):
    merge_pout.main(
        ['--perc_results'] + pout_files +
        ['--perc_merged', f'{output_dir}/pout_merged.txt', 
         '--fasta', fasta_file] + digest_params)


def run_picked_group_fdr_percolator_input(output_dir, fasta_file):
    picked_group_fdr.main(
        ['--perc_evidence', f'{output_dir}/pout_merged.txt',
         '--methods', 'picked_protein_group_no_remap', 
         '--protein_groups_out', f'{output_dir}/proteinGroups_percolator.txt',
         '--fasta', fasta_file])