from typing import List, Dict, Set
import logging

from .parsers import protein_groups as pgp


logger = logging.getLogger(__name__)


def mark_entrapment_proteins(
    peptide_to_protein_map: Dict, mq_protein_groups_file: str
) -> None:
    """
    necessary for simulated datasets with entrapments
    since the peptide-protein map is made without knowledge of which proteins
    are entrapment proteins we update the peptide_to_protein_map by adding a
    "_entrapment" suffix to the simulated entrapment proteins
    """
    if not mq_protein_groups_file:
        return

    entrapment_proteins = get_entrapment_proteins(mq_protein_groups_file)
    if len(entrapment_proteins) > 0:
        logger.info(
            "Updating peptide to protein map with entrapment protein identifiers"
        )
        for peptide, proteins in peptide_to_protein_map.items():
            peptide_to_protein_map[peptide] = mark_entrapment(
                proteins, entrapment_proteins
            )


def get_entrapment_proteins(mq_protein_groups_file: str) -> Set[str]:
    entrapment_proteins = set()
    for protein_group, _ in pgp.parse_protein_groups_file_single(
        mq_protein_groups_file
    ):
        for protein in protein_group:
            if "_entrapment" in protein:
                entrapment_proteins.add(protein)
    return entrapment_proteins


def mark_entrapment(proteins: List[str], entrapment_proteins: Set[str]) -> List[str]:
    return list(
        map(
            lambda protein: protein + "_entrapment"
            if protein + "_entrapment" in entrapment_proteins
            else protein,
            proteins,
        )
    )


def is_entrapment(protein_group: List[str]) -> bool:
    return len(
        [
            1
            for x in protein_group
            if "_entrapment" in x or "Random_" in x or "REV__" in x or "mimic" in x
        ]
    ) == len(protein_group)
