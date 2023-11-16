from __future__ import annotations

import collections
import logging
from dataclasses import dataclass, field
from typing import List

import numpy as np

from . import helpers
from .parsers import tsv

# for type hints only
from .peptide_info import ProteinGroupPeptideInfos
from .protein_groups import ProteinGroups


logger = logging.getLogger(__name__)


PROTEIN_GROUP_HEADERS = [
    "Protein IDs",
    "Majority protein IDs",
    "Peptide counts (unique)",
    "Protein names",
    "Gene names",
    "Fasta headers",
    "Best peptide",
    "Number of proteins",
    "Q-value",
    "Score",
    "Reverse",
    "Potential contaminant",
]


@dataclass
class ProteinGroupResult:
    proteinIds: str = ""
    majorityProteinIds: str = ""
    peptideCountsUnique: str = ""
    proteinNames: str = ""
    geneNames: str = ""
    fastaHeaders: str = ""
    bestPeptide: str = ""
    numberOfProteins: int = 0
    qValue: float = np.nan
    score: float = np.nan
    reverse: str = ""
    potentialContaminant: str = ""
    precursorQuants: List[float] = field(default_factory=list)
    extraColumns: List[float] = field(default_factory=list)
    
    def extend(self, e):
        self.extraColumns.extend(e)
    
    def append(self, a):
        self.extraColumns.append(a)
        
    @classmethod
    def from_protein_group(cls, proteinGroup, peptideScores, reportedFdr, proteinScore, scoreCutoff, proteinAnnotations, keep_all_proteins):
        numUniquePeptidesPerProtein = cls._get_peptide_counts(peptideScores, scoreCutoff)
        peptideCountsUnique = [numUniquePeptidesPerProtein[p] for p in proteinGroup]
        if sum(peptideCountsUnique) == 0 and not keep_all_proteins:
            return None

        proteinGroup, peptideCountsUnique = zip(*[(p, num_peptides) for p, num_peptides in zip(proteinGroup, peptideCountsUnique) if num_peptides > 0 or keep_all_proteins])
        
        bestPeptide = sorted([(p[0], p[1]) for p in peptideScores])[0][1]
        majorityProteinIds = ";".join([p for p, num_peptides in zip(proteinGroup, peptideCountsUnique) if num_peptides >= max(peptideCountsUnique) / 2])
        numberOfProteins = len(proteinGroup)

        peptideCountsUnique = ";".join(map(str, peptideCountsUnique))
        proteinIds = ";".join(proteinGroup)
                
        proteinNames, geneNames, fastaHeaders = list(), list(), list()
        for p in proteinGroup:
            if p in proteinAnnotations:
                proteinName, geneName, fastaHeader = proteinAnnotations[p]
                if proteinName not in proteinNames:
                    proteinNames.append(proteinName)
                if geneName not in geneNames:
                    geneNames.append(geneName)
                if fastaHeader not in fastaHeaders:
                    fastaHeaders.append(fastaHeader)
        
        proteinNames = ";".join(proteinNames)
        geneNames = ";".join(geneNames)
        fastaHeaders = ";".join(fastaHeaders)
        
        qValue = reportedFdr
        score = proteinScore
        reverse = '+' if helpers.isDecoy(proteinGroup) else ''
        potentialContaminant = '+' if helpers.is_contaminant(proteinGroup) else ''
        return cls(proteinIds, majorityProteinIds, peptideCountsUnique, 
                             proteinNames, geneNames, fastaHeaders, bestPeptide, 
                             numberOfProteins, qValue, score, reverse, potentialContaminant, [], [])
    
    @staticmethod
    def _get_peptide_counts(scorePeptidePairs, scoreCutoff):
        proteinPeptideCount = collections.defaultdict(int)
        seenPeptides = set()
        for PEP, peptide, proteins in sorted(scorePeptidePairs):
            if PEP > scoreCutoff:
                break
            if peptide not in seenPeptides:
                seenPeptides.add(peptide)
                for protein in proteins:
                    proteinPeptideCount[protein] += 1
        return proteinPeptideCount
        
    def to_list(self):
        return [self.proteinIds, 
                self.majorityProteinIds, 
                self.peptideCountsUnique,
                self.proteinNames,
                self.geneNames,
                self.fastaHeaders,
                self.bestPeptide,
                self.numberOfProteins,
                self.qValue,
                self.score,
                self.reverse,
                self.potentialContaminant] + ["%.0f" % (x) if not type(x) == str else x for x in self.extraColumns]


class ProteinGroupResults:
    headers = List[str]
    protein_group_results: List[ProteinGroupResult] = field(default_factory=list)
    
    def __init__(self, protein_group_results: ProteinGroupResults=None):
        """
        NOTE: cannot use empty list as default argument: https://docs.python-guide.org/writing/gotchas/
        """
        self.protein_group_results = []
        if protein_group_results is not None: 
            self.protein_group_results = protein_group_results
        self.headers = PROTEIN_GROUP_HEADERS
        
    def __iter__(self):
        return iter(self.protein_group_results)
    
    def __next__(self):
        return next(self.protein_group_results)
    
    def __len__(self):
        return len(self.protein_group_results)
    
    def __getitem__(self, indices):
        return self.protein_group_results[indices]
    
    def append_header(self, header: str) -> None:
        self.headers.append(header)
    
    def write(self, output_file: str) -> None:
        writer = tsv.get_tsv_writer(output_file)
        writer.writerow(self.headers)
        for proteinRow in self.protein_group_results:
            writer.writerow(proteinRow.to_list())

    @classmethod
    def from_protein_groups(cls, 
            proteinGroups: ProteinGroups, 
            proteinGroupPeptideInfos: ProteinGroupPeptideInfos, 
            proteinScores: List[float], 
            reportedQvals: List[float], 
            scoreCutoff: float, 
            proteinAnnotations,
            keep_all_proteins: bool) -> ProteinGroupResults:
        protein_group_results = list()
        for proteinGroup, peptideScores, proteinScore, reportedFdr in zip(proteinGroups, proteinGroupPeptideInfos, proteinScores, reportedQvals):
            if helpers.isObsolete(proteinGroup):
                continue
            pgr = ProteinGroupResult.from_protein_group(proteinGroup, peptideScores, reportedFdr, proteinScore, scoreCutoff, proteinAnnotations, keep_all_proteins)
            if pgr is not None: # protein groups can get filtered out when they do not have any PSMs below the PSM FDR cutoff
                protein_group_results.append(pgr)
        
        return cls(protein_group_results)
