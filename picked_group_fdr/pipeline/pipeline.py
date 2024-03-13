import logging
import os
from typing import List

import mokapot
import numpy as np
from joblib import parallel_backend

from ..digestion_params import DigestionParams, digestion_params_list_to_arg_list

from .. import picked_group_fdr
from . import andromeda2pin, merge_pout, filter_fdr_maxquant
from . import update_evidence_from_pout as update_evidence


logger = logging.getLogger(__name__)


def run_picked_group_fdr_all(
    evidence_files: List[str],
    pout_files: List[str],
    fasta_files: List[str],
    output_dir: str,
    digest_params_list: List[DigestionParams],
    input_type: str,
    do_quant: bool,
    lfq_min_peptide_ratios: int,
    fdr_cutoff: float,
    suppress_missing_peptide_warning: bool = False,
):
    try:
        if len(output_dir) == 0:
            raise RuntimeError("Please specify an output folder")
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

        protein_groups_out = f"{output_dir}/proteinGroups.txt"
        protein_groups_filtered_out = (
            f"{output_dir}/proteinGroups.fdr{fdr_cutoff*100:g}.txt"
        )
        if input_type in ["rescoring", "mq"]:
            pout_input_type = "prosit"
            if input_type == "mq":
                pout_input_type = "andromeda"
                pin_files = run_andromeda_to_pin(
                    evidence_files,
                    fasta_files,
                    output_dir,
                    digest_params_list,
                    suppress_missing_peptide_warning,
                )
                pout_files = run_mokapot(pin_files, output_dir)

            evidence_files_rescored = [
                f"{output_dir}/evidence_{idx}.txt" for idx in range(len(evidence_files))
            ]
            run_update_evidence(
                evidence_files,
                pout_files,
                evidence_files_rescored,
                pout_input_type,
                suppress_missing_peptide_warning,
            )
            run_picked_group_fdr(
                evidence_files_rescored,
                protein_groups_out,
                fasta_files,
                digest_params_list,
                do_quant,
                lfq_min_peptide_ratios,
                suppress_missing_peptide_warning,
            )
        elif input_type == "percolator_remap":  # currently not accessible by the GUI
            pout_merged = f"{output_dir}/pout_merged.txt"
            run_merge_pout_remap(
                pout_files,
                fasta_files,
                output_dir,
                digest_params_list,
                suppress_missing_peptide_warning,
            )
            run_picked_group_fdr_percolator_input_remap(
                pout_merged,
                fasta_files,
                protein_groups_out,
                suppress_missing_peptide_warning,
            )
        elif input_type == "percolator":
            run_picked_group_fdr_percolator_input(
                pout_files, protein_groups_out, suppress_missing_peptide_warning
            )
        else:
            logger.error(
                f"Error while running Picked Group FDR, unknown input type: {input_type}."
            )

        run_filter_fdr_maxquant(
            [protein_groups_out], protein_groups_filtered_out, fdr_cutoff
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
    suppress_missing_peptide_warning: bool,
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
            + suppress_missing_peptide_warning_flag(
                suppress_missing_peptide_warning,
            )
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
                f"{output_dir}/{idx}.mokapot.decoy.psms.txt",
            ]
        )
    return pout_files


def run_update_evidence(
    evidence_files: List[str],
    pout_files: List[str],
    evidence_files_rescored: List[str],
    pout_input_type: str,
    suppress_missing_peptide_warning: bool,
):
    if len(evidence_files) != len(evidence_files_rescored):
        logger.error("Unequal number of input and output evidence files.")

    for evidence_file, evidence_file_rescored in zip(
        evidence_files, evidence_files_rescored
    ):
        update_evidence.main(
            ["--mq_evidence", evidence_file, "--perc_results"]
            + pout_files
            + [
                "--mq_evidence_out",
                evidence_file_rescored,
                "--pout_input_type",
                pout_input_type,
            ]
            + suppress_missing_peptide_warning_flag(
                suppress_missing_peptide_warning,
            )
        )


def run_picked_group_fdr(
    evidence_files: List[str],
    protein_groups_out: str,
    fasta_files: List[str],
    digest_params_list: List[DigestionParams],
    do_quant: bool,
    lfq_min_peptide_ratios: int,
    suppress_missing_peptide_warning: bool,
):
    digest_params_str = digestion_params_list_to_arg_list(digest_params_list)
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
            protein_groups_out,
            "--fasta",
        ]
        + fasta_files
        + digest_params_str
        + quant_flags
        + suppress_missing_peptide_warning_flag(
            suppress_missing_peptide_warning,
        )
    )


def run_merge_pout(
    pout_files: List[str],
    pout_merged: str,
    suppress_missing_peptide_warning: bool,
):
    merge_pout.main(
        ["--perc_results"]
        + pout_files
        + ["--perc_merged", pout_merged]
        + suppress_missing_peptide_warning_flag(
            suppress_missing_peptide_warning,
        )
    )


def run_merge_pout_remap(
    pout_files: List[str],
    fasta_files: List[str],
    pout_merged: str,
    digest_params_list: List[DigestionParams],
    suppress_missing_peptide_warning: bool,
):
    digest_params_str = digestion_params_list_to_arg_list(digest_params_list)
    merge_pout.main(
        ["--perc_results"]
        + pout_files
        + ["--perc_merged", pout_merged, "--fasta"]
        + fasta_files
        + digest_params_str
        + suppress_missing_peptide_warning_flag(
            suppress_missing_peptide_warning,
        )
    )


def run_picked_group_fdr_percolator_input(
    pout_files: List[str],
    protein_groups_out: str,
    suppress_missing_peptide_warning: bool,
):
    picked_group_fdr.main(
        [
            "--perc_evidence",
        ]
        + pout_files
        + [
            "--methods",
            "picked_protein_group_no_remap",
            "--protein_groups_out",
            protein_groups_out,
        ]
        + suppress_missing_peptide_warning_flag(
            suppress_missing_peptide_warning,
        )
    )


def run_picked_group_fdr_percolator_input_remap(
    pout_merged: str,
    fasta_files: List[str],
    protein_groups_out: str,
    suppress_missing_peptide_warning: bool,
):
    picked_group_fdr.main(
        [
            "--perc_evidence",
            pout_merged,
            "--methods",
            "picked_protein_group_no_remap",
            "--protein_groups_out",
            protein_groups_out,
            "--fasta",
        ]
        + fasta_files
        + suppress_missing_peptide_warning_flag(
            suppress_missing_peptide_warning,
        )
    )


def run_filter_fdr_maxquant(
    protein_groups_files: List[str], protein_groups_out: str, fdr_cutoff: float = 0.01
):
    filter_fdr_maxquant.main(
        ["--mq_protein_groups"]
        + protein_groups_files
        + [
            "--mq_protein_groups_out",
            protein_groups_out,
            "--fdr_cutoff",
            str(fdr_cutoff),
        ]
    )


def suppress_missing_peptide_warning_flag(
    suppress_missing_peptide_warning,
) -> List[str]:
    if suppress_missing_peptide_warning:
        return ["--suppress_missing_peptide_warning"]
    return []
