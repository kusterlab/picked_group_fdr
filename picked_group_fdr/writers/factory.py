from __future__ import annotations

import argparse
from typing import Callable, Dict, List

from .. import digest
from .. import protein_annotation
from .. import writers
from ..columns.triqler import init_triqler_params
from ..protein_groups import ProteinGroups

# for type hints only
from .. import results


def get_protein_groups_output_writer(
    protein_group_results: results.ProteinGroupResults,
    output_format: str,
    args: argparse.Namespace,
    protein_annotations: Dict[str, protein_annotation.ProteinAnnotation],
    parse_id: Callable,
    peptide_to_protein_maps: List[digest.PeptideToProteinMap],
):
    if args.output_format == "fragpipe":
        protein_groups = ProteinGroups.from_protein_group_results(protein_group_results)
        protein_groups.create_index()
        return writers.FragPipeCombinedProteinWriter(
            protein_groups,
            protein_annotations,
            args.lfq_min_peptide_ratios,
            args.lfq_stabilize_large_ratios,
            args.num_threads,
        )
    elif args.output_format == "maxquant":
        db = "target" if args.fasta_contains_decoys else "concat"
        protein_sequences = digest.get_protein_sequences(
            args.fasta, db=db, parse_id=parse_id
        )
        num_ibaq_peptides_per_protein = (
            digest.get_num_ibaq_peptides_per_protein_from_args(
                args, peptide_to_protein_maps
            )
        )
        triqler_params = init_triqler_params(args.file_list_file)
        return writers.MaxQuantProteinGroupsWriter(
            num_ibaq_peptides_per_protein,
            protein_annotations,
            protein_sequences,
            args.lfq_min_peptide_ratios,
            args.lfq_stabilize_large_ratios,
            args.num_threads,
            triqler_params,
        )
    
    raise ValueError(f"Unknown output format: {args.output_format}.")