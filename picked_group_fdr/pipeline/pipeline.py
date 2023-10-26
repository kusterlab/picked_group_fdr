import logging
import os
from typing import List

import mokapot
import numpy as np
from joblib import parallel_backend

from ..digestion_params import DigestionParams, digestion_params_list_to_arg_list

from .. import picked_group_fdr
from . import andromeda2pin, merge_pout
from . import update_evidence_from_pout as update_evidence


logger = logging.getLogger()


def run_picked_group_fdr_all(
    evidence_files: List[str],
    pout_files: List[str],
    fasta_files: List[str],
    output_dir: str,
    digest_params_list: List[DigestionParams],
    input_type: str,
    do_quant: bool,
    lfq_min_peptide_ratios: int,
):
    digest_params = digestion_params_list_to_arg_list(digest_params_list)
    try:
        if len(output_dir) == 0:
            raise RuntimeError("Please specify an output folder")
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

        if input_type in ["rescoring", "mq"]:
            pout_input_type = "prosit"
            if input_type == "mq":
                pout_input_type = "andromeda"
                pin_files = run_andromeda_to_pin(
                    evidence_files, fasta_files, output_dir, digest_params_list
                )
                pout_files = run_mokapot(pin_files, output_dir)

            evidence_files_rescored = run_update_evidence(
                evidence_files, pout_files, output_dir, pout_input_type
            )
            run_picked_group_fdr(
                evidence_files_rescored,
                output_dir,
                fasta_files,
                digest_params,
                do_quant,
                lfq_min_peptide_ratios,
            )
        elif input_type == "percolator_remap":  # currently not accessible by the GUI
            run_merge_pout_remap(pout_files, fasta_files, output_dir, digest_params)
            run_picked_group_fdr_percolator_input_remap(output_dir, fasta_files)
        elif input_type == "percolator":
            run_picked_group_fdr_percolator_input(pout_files, output_dir)
        else:
            logger.error(
                f"Error while running Picked Group FDR, unknown input type: {input_type}."
            )

    except SystemExit as e:
        logger.error(
            f"Error while running Picked Group FDR, exited with error code: {e}."
        )
    except Exception as e:
        logger.error(f"Error while running Picked Group FDR: {e}")


def run_andromeda_to_pin(
    evidence_files: List[str],
    fasta_files: List[str],
    output_dir: str,
    digest_params_list: List[DigestionParams],
):
    pin_files = list()
    for idx, (evidence_file, digest_params) in enumerate(
        zip(evidence_files, digest_params_list)
    ):
        digest_params_str = digestion_params_list_to_arg_list([digest_params])
        pin_file = f"{output_dir}/pin_{idx}.tab"
        andromeda2pin.main(
            [evidence_file, "--outputTab", pin_file, "--databases"]
            + fasta_files
            + digest_params_str
        )
        pin_files.append(pin_file)
    return pin_files


def run_mokapot(pin_files: List[str], output_dir: str):
    pout_files = []
    for idx, pin_file in enumerate(pin_files):
        np.random.seed(0)  # TODO: Make seed configurable
        psms = mokapot.read_pin(pin_file)
        with parallel_backend("threading"):
            results, _ = mokapot.brew(psms)
        results.to_txt(file_root=f"{idx}", dest_dir=output_dir, decoys=True)
        pout_files.extend(
            [
                f"{output_dir}/{idx}.mokapot.psms.txt",
                f"{output_dir}/{idx}.mokapot.decoys.psms.txt",
            ]
        )
    return pout_files


def run_update_evidence(
    evidence_files: List[str],
    pout_files: List[str],
    output_dir: str,
    pout_input_type: str,
):
    evidence_files_rescored = list()
    for idx, evidence_file in enumerate(evidence_files):
        evidence_file_rescored = f"{output_dir}/evidence_{idx}.txt"
        update_evidence.main(
            ["--mq_evidence", evidence_file, "--perc_results"]
            + pout_files
            + [
                "--mq_evidence_out",
                evidence_file_rescored,
                "--pout_input_type",
                pout_input_type,
            ]
        )
        evidence_files_rescored.append(evidence_file_rescored)
    return evidence_files_rescored


def run_picked_group_fdr(
    evidence_files: List[str],
    output_dir: str,
    fasta_files: List[str],
    digest_params: List[str],
    do_quant: bool,
    lfq_min_peptide_ratios: int,
):
    quant_flags = []
    if do_quant:
        quant_flags = [
            "--do_quant",
            "--lfq_min_peptide_ratios",
            str(lfq_min_peptide_ratios),
        ]
    picked_group_fdr.main(
        ["--mq_evidence"]
        + evidence_files
        + [
            "--methods",
            "picked_protein_group_mq_input",
            "--do_quant",
            "--protein_groups_out",
            f"{output_dir}/proteinGroups_percolator.txt",
            "--fasta",
        ]
        + fasta_files
        + digest_params
        + quant_flags
    )


def run_merge_pout(
    pout_files: List[str],
    output_dir: str,
):
    merge_pout.main(
        ["--perc_results"]
        + pout_files
        + ["--perc_merged", f"{output_dir}/pout_merged.txt"]
    )


def run_merge_pout_remap(
    pout_files: List[str],
    fasta_files: List[str],
    output_dir: str,
    digest_params: List[str],
):
    merge_pout.main(
        ["--perc_results"]
        + pout_files
        + ["--perc_merged", f"{output_dir}/pout_merged.txt", "--fasta"]
        + fasta_files
        + digest_params
    )


def run_picked_group_fdr_percolator_input(pout_files: List[str], output_dir: str):
    picked_group_fdr.main(
        [
            "--perc_evidence",
        ]
        + pout_files
        + [
            "--methods",
            "picked_protein_group_no_remap",
            "--protein_groups_out",
            f"{output_dir}/proteinGroups_percolator.txt",
        ]
    )


def run_picked_group_fdr_percolator_input_remap(
    output_dir: str, fasta_files: List[str]
):
    picked_group_fdr.main(
        [
            "--perc_evidence",
            f"{output_dir}/pout_merged.txt",
            "--methods",
            "picked_protein_group_no_remap",
            "--protein_groups_out",
            f"{output_dir}/proteinGroups_percolator.txt",
            "--fasta",
        ]
        + fasta_files
    )
