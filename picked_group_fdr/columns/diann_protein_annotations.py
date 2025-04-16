from __future__ import annotations

from typing import Dict
import logging

from picked_group_fdr import protein_groups as pg

from ..protein_annotation import ProteinAnnotation  # TODO: get rid of this import
from .base import ProteinGroupColumns

# for type hints only
from .. import results


logger = logging.getLogger(__name__)

DIANN_PROTEIN_ANNOTATION_HEADERS = [
    "Protein.Group",
    "Protein.Names",
    "Genes",
    "First.Protein.Description",
]


class DiannProteinAnnotationsColumns(ProteinGroupColumns):
    protein_annotations: Dict[str, ProteinAnnotation]

    def __init__(self, protein_annotations):
        self.protein_annotations = protein_annotations

    def append_headers(self, protein_group_results: results.ProteinGroupResults):
        protein_group_results.append_headers(DIANN_PROTEIN_ANNOTATION_HEADERS)

    def append_columns(
        self,
        protein_group_results: results.ProteinGroupResults,
        post_err_prob_cutoff: float,
    ):
        logger.info("Adding protein annotations")
        for pgr in protein_group_results:
            (
                protein_uniprot_ids,
                protein_names,
                gene_names,
                first_protein_description,
            ) = (list(), list(), list(), list())
            for p in pgr.proteinIds.split(";"):
                if p not in self.protein_annotations:
                    continue

                protein_annotation = self.protein_annotations[p]
                if protein_annotation.uniprot_id not in protein_uniprot_ids:
                    protein_uniprot_ids.append(protein_annotation.uniprot_id)

                if protein_annotation.entry_name not in protein_names:
                    protein_names.append(protein_annotation.entry_name)

                if (
                    protein_annotation.gene_name is not None
                    and protein_annotation.gene_name not in gene_names
                ):
                    gene_names.append(protein_annotation.gene_name)

                if len(first_protein_description) == 0:
                    first_protein_description.append(protein_annotation.description)

            protein_uniprot_ids = ";".join(protein_uniprot_ids)
            protein_names = ";".join(protein_names)
            gene_names = ";".join(gene_names)
            first_protein_description = ";".join(first_protein_description)

            pgr.extend(
                [
                    protein_uniprot_ids,
                    protein_names,
                    gene_names,
                    first_protein_description,
                ]
            )
