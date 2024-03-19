from dataclasses import dataclass
from typing import Callable, Dict, Iterator, Optional, Tuple
import os
import logging

from picked_group_fdr import digest


logger = logging.getLogger(__name__)


"""
https://www.uniprot.org/help/protein_existence

1. Experimental evidence at protein level
2. Experimental evidence at transcript level
3. Protein inferred from homology
4. Protein predicted
5. Protein uncertain
"""


@dataclass
class ProteinAnnotation:
    """
    - id (sp|P00167|CYB5_HUMAN)
    - fasta_header (sp|P00167|CYB5_HUMAN Cytochrome b5 OS=Homo sapiens OX=9606 GN=CYB5A PE=1 SV=2)
    - uniprot_id (P00167)
    - entry_name (CYB5_HUMAN)
    - gene_name (CYB5A)
    - length
    - organism (Homo sapiens OX=9606)
    - description (Cytochrome b5)
    - existence (1:Experimental evidence at protein level)
    """

    id: str
    fasta_header: str
    uniprot_id: Optional[str] = None
    entry_name: Optional[str] = None
    gene_name: Optional[str] = None
    length: int = 0
    organism: Optional[str] = None
    description: Optional[str] = None
    existence: Optional[int] = None


def parse_protein_name_func(fasta_header: str) -> str:
    return " ".join(fasta_header.split(" OS=")[0].split(" ")[1:])


def parse_organism(fasta_header: str) -> Optional[str]:
    if " OS=" in fasta_header:
        return fasta_header.split(" OS=")[1].split(" GN=")[0]
    return None


def parse_protein_existence_level(fasta_header: str) -> Optional[str]:
    if " PE=" in fasta_header:
        return int(fasta_header.split(" PE=")[1].split(" ")[0])
    return None


def parse_gene_name_func(fasta_header: str) -> Optional[str]:
    if " GN=" in fasta_header:
        return fasta_header.split(" GN=")[1].split(" ")[0]
    return None


def parse_uniprot_id(fasta_id: str) -> str:
    protein_id = digest.parse_until_first_space(fasta_id)
    if "|" in protein_id:
        return protein_id.split("|")[1]
    else:
        return protein_id


def parse_entry_name(fasta_header: str) -> str:
    protein_id = digest.parse_until_first_space(fasta_header)
    if "|" in protein_id and protein_id.count("|") >= 2:
        return protein_id.split("|")[2]
    else:
        return protein_id


def parse_fasta_header(fasta_header: str) -> str:
    return fasta_header


def read_fasta_proteins(
    file_path: str,
    db: str = "concat",
    parse_id: Callable[[str], str] = digest.parse_until_first_space,
    parse_protein_name: Callable[[str], str] = parse_protein_name_func,
    parse_gene_name: Callable[[str], str] = parse_gene_name_func,
) -> Iterator[ProteinAnnotation]:
    if not os.path.isfile(file_path):
        raise FileNotFoundError(
            f"Could not find fasta file {file_path}. Please make sure you provided a valid fasta file."
        )

    for fasta_header, protein_sequence in digest.read_fasta(
        file_path, db=db, parse_id=parse_fasta_header
    ):
        yield ProteinAnnotation(
            id=parse_id(fasta_header),
            fasta_header=fasta_header,
            uniprot_id=parse_uniprot_id(fasta_header),
            entry_name=parse_entry_name(fasta_header),
            gene_name=parse_gene_name(fasta_header),
            length=len(protein_sequence),
            organism=parse_organism(fasta_header),
            description=parse_protein_name(fasta_header),
            existence=parse_protein_existence_level(fasta_header),
        )


def get_protein_annotations_single(
    fasta_file: str, **kwargs
) -> Dict[str, ProteinAnnotation]:
    protein_annotations = dict()

    if not fasta_file:
        return protein_annotations

    for protein_annotation in read_fasta_proteins(fasta_file, **kwargs):
        if protein_annotation.id not in protein_annotations:
            protein_annotations[protein_annotation.id] = protein_annotation

    return protein_annotations


def get_protein_annotations_multiple(
    fasta_files: Iterator[str], **kwargs
) -> Dict[str, ProteinAnnotation]:
    protein_annotations = dict()
    for fasta_file in fasta_files:
        protein_annotations = {
            **protein_annotations,
            **get_protein_annotations_single(fasta_file, **kwargs),
        }

    return protein_annotations


def has_gene_names(
    protein_annotations: Dict[str, ProteinAnnotation], min_ratio_with_genes: float
) -> bool:
    counts = sum(
        1
        for pa in protein_annotations.values()
        if pa.gene_name is not None and len(pa.gene_name) > 0
    )
    return counts / len(protein_annotations) > min_ratio_with_genes


def get_protein_annotations(
    fasta: str, fasta_contains_decoys: bool, use_gene_level: bool
) -> Tuple[Dict[str, ProteinAnnotation], bool]:
    protein_annotations = dict()
    use_pseudo_genes = False
    if fasta is None:
        return protein_annotations, use_pseudo_genes
    
    db = "target" if fasta_contains_decoys else "concat"
    protein_annotations = get_protein_annotations_multiple(
        fasta, db=db, parse_id=digest.parse_until_first_space
    )
    if use_gene_level:
        if has_gene_names(protein_annotations, min_ratio_with_genes=0.5):
            protein_annotations = get_protein_annotations_multiple(
                fasta, db=db, parse_id=parse_gene_name_func
            )
        else:
            logger.warning(
                (
                    "Found >50% of proteins without gene names in the "
                    "fasta file, will infer pseudo-genes based on "
                    "shared peptides instead."
                )
            )
            use_pseudo_genes = True
    return protein_annotations, use_pseudo_genes
