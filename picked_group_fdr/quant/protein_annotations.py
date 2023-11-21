from typing import Dict
import logging

from ..parsers.maxquant import MQ_PROTEIN_ANNOTATION_HEADERS
from ..protein_annotation import ProteinAnnotation
from ..results import ProteinGroupResults
from .base import ProteinGroupColumns


logger = logging.getLogger(__name__)


class ProteinAnnotationsColumns(ProteinGroupColumns):
    protein_annotations: Dict[str, ProteinAnnotation]

    def __init__(self, protein_annotations):
        self.protein_annotations = protein_annotations

    def append_headers(self, protein_group_results: ProteinGroupResults, experiments):
        protein_group_results.append_headers(
            MQ_PROTEIN_ANNOTATION_HEADERS
        )

    def append_columns(
        self,
        protein_group_results: ProteinGroupResults,
        experimentToIdxMap,
        postErrProbCutoff,
    ):
        logger.info("Adding protein annotations")
        for pgr in protein_group_results:
            proteinNames, geneNames, fastaHeaders = list(), list(), list()
            for p in pgr.proteinIds.split(";"):
                if p not in self.protein_annotations:
                    continue

                protein_annotation = self.protein_annotations[p]
                if protein_annotation.id not in proteinNames:
                    proteinNames.append(protein_annotation.id)
                if (
                    protein_annotation.gene_name is not None
                    and protein_annotation.gene_name not in geneNames
                ):
                    geneNames.append(protein_annotation.gene_name)
                if protein_annotation.fasta_header not in fastaHeaders:
                    fastaHeaders.append(protein_annotation.fasta_header)

            proteinNames = ";".join(proteinNames)
            geneNames = ";".join(geneNames)
            fastaHeaders = ";".join(fastaHeaders)

            pgr.extend([proteinNames, geneNames, fastaHeaders])
