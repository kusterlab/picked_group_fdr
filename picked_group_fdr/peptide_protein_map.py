import logging
import argparse
from typing import List

from . import digest
from . import entrapment
from . import protein_annotation
from .digestion_params import DigestionParams, get_digestion_params_list

logger = logging.getLogger(__name__)


def get_peptide_to_protein_maps(
    fasta_file: str,
    peptide_protein_map_files: List[str],
    digestion_params_list: List[DigestionParams],
    mq_protein_groups_file: str,
    **kwargs,
):
    peptide_to_protein_maps = list()
    if fasta_file:
        for digestion_params in digestion_params_list:
            peptide_to_protein_maps.append(
                digest.get_peptide_to_protein_map_from_params(
                    fasta_file, [digestion_params], **kwargs
                )
            )
            entrapment.mark_entrapment_proteins(
                peptide_to_protein_maps[-1], mq_protein_groups_file
            )
    elif peptide_protein_map_files:
        logger.info("Loading peptide to protein map")
        for peptide_protein_map_file in peptide_protein_map_files:
            peptide_to_protein_maps.append(
                digest.get_peptide_to_protein_map_from_file(
                    peptide_protein_map_file, use_hash_key=False
                )
            )
    else:
        raise ValueError(
            (
                "No fasta or peptide to protein mapping file detected, please"
                "specify either the --fasta or --peptide_protein_map flags."
            )
        )
    return peptide_to_protein_maps


def get_peptide_to_protein_maps_from_args(
    args: argparse.Namespace,
    use_pseudo_genes: bool,
) -> List[digest.PeptideToProteinMap]:
    parse_id = digest.parse_until_first_space
    if args.gene_level and not use_pseudo_genes:
        parse_id = protein_annotation.parse_gene_name_func
    elif args.fasta_use_uniprot_id:
        parse_id = protein_annotation.parse_uniprot_id

    digestion_params_list = get_digestion_params_list(args)
    return get_peptide_to_protein_maps(
        args.fasta,
        args.peptide_protein_map,
        digestion_params_list,
        args.mq_protein_groups,
        parse_id=parse_id,
    )
