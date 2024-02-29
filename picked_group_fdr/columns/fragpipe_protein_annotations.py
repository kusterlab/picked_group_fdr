from __future__ import annotations

from typing import Dict
import logging

from picked_group_fdr import protein_groups as pg

from ..protein_annotation import ProteinAnnotation  # TODO: get rid of this import
from .base import ProteinGroupColumns

# for type hints only
from .. import results


logger = logging.getLogger(__name__)

FRAGPIPE_PROTEIN_ANNOTATION_HEADERS = [
    "Protein",
    "Protein ID",
    "Entry Name",
    "Gene",
    "Length",
    "Organism",
    "Protein Description",
    "Protein Existence",
]


class FragpipeProteinAnnotationsColumns(ProteinGroupColumns):
    protein_groups: pg.ProteinGroups
    protein_annotations: Dict[str, ProteinAnnotation]

    def __init__(
        self,
        protein_groups: pg.ProteinGroups,
        protein_annotations: Dict[str, ProteinAnnotation],
    ):
        self.protein_groups = protein_groups
        self.protein_annotations = protein_annotations

    def append_headers(
        self,
        protein_group_results: results.ProteinGroupResults,
    ):
        protein_group_results.append_headers(FRAGPIPE_PROTEIN_ANNOTATION_HEADERS)

    def append_columns(
        self,
        protein_group_results: results.ProteinGroupResults,
        post_err_prob_cutoff: float,
    ):
        logger.info("Adding protein annotations")
        for pgr in protein_group_results:
            row_protein_groups = self.protein_groups.get_protein_groups(
                pgr.proteinIds.split(";")
            )
            leading_protein = row_protein_groups[0][0]

            protein_annotation = self.protein_annotations.get(
                leading_protein,
                ProteinAnnotation(id=leading_protein, fasta_header=leading_protein),
            )
            protein_annotation_list = [
                protein_annotation.id,
                protein_annotation.uniprot_id,
                protein_annotation.entry_name,
                protein_annotation.gene_name,
                protein_annotation.length,
                protein_annotation.organism,
                protein_annotation.description,
                protein_annotation.existence,
            ]
            protein_annotation_list = [
                a if a is not None else "" for a in protein_annotation_list
            ]
            pgr.extend(protein_annotation_list)
