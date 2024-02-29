from __future__ import annotations

import logging

from .base import ProteinGroupColumns

# for type hints only
from .. import results

logger = logging.getLogger(__name__)


class IndistinguishableProteinsColumns(ProteinGroupColumns):
    def append_headers(
        self,
        protein_group_results: results.ProteinGroupResults,
    ) -> None:
        protein_group_results.append_header("Indistinguishable Proteins")

    def append_columns(
        self,
        protein_group_results: results.ProteinGroupResults,
        post_err_prob_cutoff: float,
    ) -> None:
        logger.info("Doing quantification: Indistinguishable Proteins")
        for pgr in protein_group_results:
            highest_peptide_count = 0
            indistinguishable_proteins = []
            for protein_id, peptide_counts in zip(
                pgr.proteinIds.split(";"), pgr.peptideCountsUnique.split(";")
            ):
                if int(peptide_counts) >= highest_peptide_count:
                    if highest_peptide_count > 0:
                        indistinguishable_proteins.append(protein_id)
                    highest_peptide_count = int(peptide_counts)
            pgr.append(", ".join(indistinguishable_proteins))
