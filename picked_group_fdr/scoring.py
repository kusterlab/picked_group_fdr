from abc import ABC, abstractmethod
from typing import Callable, Dict, List
import logging
import hashlib

import numpy as np

from . import helpers
from . import fdr
from . import digest
from .parsers import parsers
from .parsers import psm
from .parsers import percolator
from .parsers import maxquant
from .parsers import fragpipe
from .parsers import sage
from .observed_peptides import ObservedPeptides
from .protein_groups import ProteinGroups
from .peptide_info import ProteinGroupPeptideInfos, PeptideInfoList


logger = logging.getLogger(__name__)


"""
N.B.: Higher protein score indicates more confident identification

scoreType:  MQ_protein = MaxQuant protein score from proteinGroups.txt (see below)
            multPEP = emulate MaxQuant protein score with multPEP with dividing of PEPs by fixed value
            bestPEP = best scoring peptide based on PEP
            Andromeda = raw Andromeda score from evidence.txt

razor: True = use Occam's razor to decide which protein to assign to a shared peptide

useSharedPeptides: True = use peptides shared by multiple proteins / protein groups
                   False = discard peptides shared by multiple proteins / protein groups
"""


class ProteinScore(ABC):
    @abstractmethod
    def calculate_score(self, scorePeptidePairs):
        pass
    
    @abstractmethod
    def can_do_protein_group_rescue(self):
        pass
    
    @abstractmethod
    def get_score_column(self, percolator_input):
        pass
        
    @abstractmethod
    def short_description(self):
        pass
    
    def optimize_hyperparameters(self, protein_groups, protein_group_peptide_infos):
        pass


class MQProteinScore(ProteinScore):
    mq_protein_groups_file: str
    
    def __init__(self, mq_protein_groups_file):
        self.mq_protein_groups_file = mq_protein_groups_file
    
    def calculate_score(self, scorePeptidePairs):
        raise NotImplementedError
    
    def can_do_protein_group_rescue(self):
        return False
    
    def get_score_column(self, percolator_input):
        return None
    
    def short_description(self):
        return 'm'
    
    def long_description(self):
        return 'multiplication of'
        
    def get_protein_scores_from_file(self):
        protein_group_peptide_infos = list()
        for protein_group, protein_score in parsers.parse_protein_groups_file_single(self.mq_protein_groups_file):
            protein_group_peptide_infos.append([(protein_score, "NA", protein_group)])
        return protein_group_peptide_infos


class BestAndromedaScore(ProteinScore):
    def calculate_score(self, score_peptide_pairs):
        return max([y[0] for y in score_peptide_pairs]) if len(score_peptide_pairs) > 0 else -100.0
    
    def can_do_protein_group_rescue(self):
        return False
    
    def get_score_column(self, percolator_input):
        return 'score'
    
    def short_description(self):
        return 'b'

    def long_description(self):
        return 'best'


class BestPEPScore(ProteinScore):
    def calculate_score(self, score_peptide_pairs):
        #return max([-1*y[0] for y in scorePeptidePairs]) if len(scorePeptidePairs) > 0 else -100.0
        return max([-1*np.log10(y[0] + np.nextafter(0,1)) for y in score_peptide_pairs]) if len(score_peptide_pairs) > 0 else -100.0
    
    def can_do_protein_group_rescue(self):
        return True
    
    def get_score_column(self, percolator_input):
        if percolator_input:
            return 'posterior_error_prob'
        else:
            return 'pep'
    
    def short_description(self):
        return 'b'

    def long_description(self):
        return 'best'


class MultPEPScore(ProteinScore):
    div: float
    
    def __init__(self):
        self.div = 1.0
    
    def optimize_hyperparameters(self, protein_groups, protein_group_peptide_infos):
        protein_score_tuples = list()
        for protein_group, protein_group_score_list in zip(protein_groups, protein_group_peptide_infos):
            mult_pep, num_peptides = self._get_score_and_num_peptides(protein_group_score_list)
            if num_peptides > 0 and not np.isnan(mult_pep):
                protein_score_tuples.append([mult_pep, num_peptides, helpers.is_decoy(protein_group)])
        
        protein_score_tuples = np.array(protein_score_tuples)
        min_range, max_range = 0.1, 1.0
        for step_size in [1e-1, 1e-2, 1e-3, 1e-4]:
            num_identified_targets, div = self._get_optimal_div(protein_score_tuples, np.arange(min_range, max_range, step_size))
            min_range, max_range = div - step_size*0.9, div + step_size*0.9
        
        self.div = div
        logger.info(f"Optimal division factor for multPEP score: {self.div:.4f} (#targets = {num_identified_targets})")
    
    def _get_optimal_div(self, protein_score_tuples, div_range: np.array):    
        most_identified_targets = (-np.inf, 1.0)
        for div in div_range:
            protein_scores = protein_score_tuples[:,0] + protein_score_tuples[:,1] * np.log10(div)
            
            sorted_idxs = np.argsort(protein_scores)[::-1]
            num_decoys = protein_score_tuples[:,2][sorted_idxs].cumsum()
            num_targets = (-protein_score_tuples[:,2] + 1)[sorted_idxs].cumsum()
            
            fdrs = np.divide(num_decoys+1, num_targets+1)
            qvals = fdr.fdrs_to_qvals(fdrs)
            num_identified_targets = fdr.count_below_threshold(qvals, 0.01, protein_score_tuples[:,2][sorted_idxs])
            if num_identified_targets > most_identified_targets[0]:
                most_identified_targets = (num_identified_targets, div)
        return most_identified_targets

    def calculate_score(self, score_peptide_pairs):
        mult_pep, num_peptides = self._get_score_and_num_peptides(score_peptide_pairs)
        mult_pep += np.log10(self.div) * num_peptides
        if num_peptides == 0 or np.isnan(mult_pep):
            mult_pep = -100.0
        return mult_pep
    
    def _get_score_and_num_peptides(self, score_peptide_pairs):
        mult_pep = 0.0
        seen_peptides = set()
        for PEP, peptide, _ in sorted(score_peptide_pairs):
            if peptide not in seen_peptides:
                seen_peptides.add(peptide)
                mult_pep -= np.log10(PEP + np.nextafter(0,1))
        return mult_pep, len(seen_peptides)
        
    def can_do_protein_group_rescue(self):
        return True
    
    def get_score_column(self, percolator_input):
        if percolator_input:
            return 'posterior_error_prob'
        else:
            return 'pep'
    
    def short_description(self):
        return 'm'
    
    def long_description(self):
        return 'multiplication of'


class ScoreOrigin(ABC):
    @abstractmethod
    def get_evidence_file(self, args):
        pass

    @abstractmethod
    def get_evidence_parser(self):
        pass
    
    @abstractmethod
    def remaps_peptides_to_proteins(self):
        pass
    
    @abstractmethod
    def can_do_quantification(self):
        pass

    @abstractmethod
    def short_description(self):
        pass
    
    @abstractmethod
    def long_description(self):
        pass


class PercolatorInput(ScoreOrigin):
    def get_evidence_file(self, args):
        return args.perc_evidence
    
    def get_evidence_parser(self):
        return percolator.parse_percolator_out_file
    
    def remaps_peptides_to_proteins(self):
        return False
    
    def can_do_quantification(self):
        return False
    
    def short_description(self):
        return 'p'
    
    def long_description(self):
        return 'Percolator'


class PercolatorInputRemapped(PercolatorInput):
    def remaps_peptides_to_proteins(self):
        return True
        

class MaxQuantInput(ScoreOrigin):
    def get_evidence_file(self, args):
        return args.mq_evidence
    
    def get_evidence_parser(self):
        return maxquant.parse_mq_evidence_file
    
    def remaps_peptides_to_proteins(self):
        return True
    
    def can_do_quantification(self):
        return True
    
    def short_description(self):
        return 'm'
    
    def long_description(self):
        return 'MaxQuant'


class FragPipeInput(ScoreOrigin):
    def get_evidence_file(self, args):
        return args.fragpipe_psm
    
    def get_evidence_parser(self):
        return fragpipe.parse_fragpipe_psm_file
    
    def remaps_peptides_to_proteins(self):
        return False
    
    def can_do_quantification(self):
        return True
    
    def short_description(self):
        return 'f'
    
    def long_description(self):
        return 'FragPipe'


class SageInput(ScoreOrigin):
    def get_evidence_file(self, args):
        return args.sage_results
    
    def get_evidence_parser(self):
        return sage.parse_sage_results_file
    
    def remaps_peptides_to_proteins(self):
        return False
    
    def can_do_quantification(self):
        return True
    
    def short_description(self):
        return 's'
    
    def long_description(self):
        return 'Sage'


class ProteinScoringStrategy:
    use_proteotypicity: bool
    use_razor: bool
    use_shared_peptides: bool
    protein_score: ProteinScore
    peptide_score_cutoff: float
    score_origin: ScoreOrigin
    peptide_counts_per_protein: Dict[str, int]
    best_peptide_score_per_protein: Dict[str, float]
    
    def __init__(self, score_description, mq_protein_groups_file = "", peptide_qval_cutoff = 0.01):
        if "multPEP" in score_description:
            self.protein_score = MultPEPScore()
        elif "bestPEP" in score_description:
            self.protein_score = BestPEPScore()
        elif "Andromeda" in score_description:
            self.protein_score = BestAndromedaScore()
        elif "MQ_protein" in score_description:
            self.protein_score = MQProteinScore(mq_protein_groups_file)
        else:
            raise NotImplementedError
        
        self.use_razor = "razor" in score_description
        self.use_proteotypicity = "proteotypicity" in score_description
        self.use_shared_peptides = "with_shared" in score_description
        self.peptide_qval_cutoff = peptide_qval_cutoff
        
        if "Perc" in score_description:
            # Add "remap" to the score_type if the fasta database used for protein grouping is different from the one used by Percolator
            if "remap" in score_description:
                self.score_origin = PercolatorInputRemapped()
            else:
                self.score_origin = PercolatorInput()
        elif "FragPipe" in score_description:
            self.score_origin = FragPipeInput()
        elif "Sage" in score_description:
            self.score_origin = SageInput()
        else:
            self.score_origin = MaxQuantInput()
        
    def get_evidence_file(self, args) -> str:
        return self.score_origin.get_evidence_file(args)

    def get_evidence_parser(self) -> Callable:
        return self.score_origin.get_evidence_parser()
    
    def remaps_peptides_to_proteins(self) -> bool:
        return self.score_origin.remaps_peptides_to_proteins()
    
    def can_do_quantification(self) -> bool:
        return self.score_origin.can_do_quantification()
    
    def optimize_hyperparameters(self, protein_groups, protein_group_peptide_infos) -> float:
        return self.protein_score.optimize_hyperparameters(protein_groups, protein_group_peptide_infos)
    
    def calculate_score(self, score_peptide_pairs) -> float:
        return self.protein_score.calculate_score(score_peptide_pairs)
    
    def get_score_column(self) -> str:
        return self.protein_score.get_score_column(self.score_origin.short_description() == 'p')
    
    def can_do_protein_group_rescue(self) -> bool:
        return self.protein_score.can_do_protein_group_rescue()
    
    def short_description(self) -> str:
        return self.protein_score.short_description() + self.score_origin.short_description() + 'P'
    
    def short_description_razor(self) -> str:
        return "rS" if self.use_razor else "dS"
    
    def long_description(self) -> str:
        return f"{self.protein_score.long_description()} {self.score_origin.long_description()} PEP"
    
    def long_description_razor(self) -> str:
        return "razor peptides" if self.use_razor else "discard shared peptides"
        
    def filter_proteins(self, proteins) -> List[str]:
        if self.use_razor:
            return self._retain_protein_with_most_observed_peptides(proteins)
        else:
            return proteins
    
    def set_peptide_counts_per_protein(self, 
            peptide_info_list: PeptideInfoList) -> None:
        if self.use_razor:
            observed_peptides = ObservedPeptides()
            observed_peptides.create(peptide_info_list)
            self.peptide_counts_per_protein = observed_peptides.get_peptide_counts_per_protein()
            self.best_peptide_score_per_protein = observed_peptides.get_best_peptide_score_per_protein()
            
    def _retain_protein_with_most_observed_peptides(self, proteins: List[str]) -> List[str]:
        """Retains only the protein with the most observed peptides. 
        Ties are first broken on best scoring (potentially shared) peptide and 
        otherwise by the md5 hash of the protein identifier. The latter
        ensures that we randomly select a protein, but that this happens
        consistently across different peptides."""
        num_peptides_per_protein_pairs = [(self.peptide_counts_per_protein.get(protein, 0), 
                                       -1*self.best_peptide_score_per_protein.get(protein, 1.0),
                                       hashlib.md5(protein.encode('utf-8')).hexdigest(),
                                       protein) for protein in proteins]
        pair_with_most_observed_peptides = sorted(num_peptides_per_protein_pairs, reverse = True)[0]
        protein_with_most_observed_peptides = pair_with_most_observed_peptides[-1]
        return [protein_with_most_observed_peptides]
    
    def collect_peptide_scores_per_protein(self,
            protein_groups: ProteinGroups, 
            peptide_info_list: PeptideInfoList, 
            suppress_missing_protein_warning: bool = False) -> ProteinGroupPeptideInfos:
        """Groups peptides with associated scores by protein
        
        :param protein_groups: ProteinGroups object
        :param peptide_info_list: Dict of peptide -> (score, proteins)
        :param suppress_missing_protein_warning: suppresses the warning for missing proteins
            in the proteiprotein_groupsnGroups object. This is set during the rescuing grouping procedure
            since some protein groups will have been filtered out in the rescuing step.
        :returns: lists of (score, peptide, proteins) tuples per protein group
        """
        if not self.get_score_column():
            return self.protein_score.get_protein_scores_from_file()
        
        logger.info("Assigning peptides to protein groups")
        shared_peptides, unique_peptides = 0, 0
        protein_group_peptide_infos = [list() for _ in range(len(protein_groups))]
        post_err_probs = list()
        for peptide, (score, proteins) in peptide_info_list.items():
            proteins = self.filter_proteins(proteins) # filtering for razor peptide approach
            
            protein_group_idxs = protein_groups.get_protein_group_idxs(proteins)
            if len(protein_group_idxs) == 0 and not suppress_missing_protein_warning:
                raise Exception(f"Could not find any of the proteins {proteins} in the ProteinGroups object, check if the identifier format is the same. \
                                  1st protein group in ProteinGroups object: {protein_groups.protein_groups[0]}")
            
            if not self.use_shared_peptides and helpers.is_shared_peptide(protein_group_idxs): # ignore shared peptides
                shared_peptides += 1
                continue
            
            unique_peptides += 1
            for protein_group_idx in protein_group_idxs:
                protein_group_peptide_infos[protein_group_idx].append((score, peptide, proteins))
            
            if not helpers.is_decoy(proteins) and not helpers.is_mbr(score):
                post_err_probs.append(score)
        
        self.peptide_score_cutoff = fdr.calc_post_err_prob_cutoff(post_err_probs, self.peptide_qval_cutoff)
        logger.info(f"#Precursors: Shared peptides = {shared_peptides}; Unique peptides = {unique_peptides}")
        return protein_group_peptide_infos
    

def compare_razor_peptides(mq_evidence_file, peptide_to_protein_map, protein_groups, score_type):
    """Compares the chosen protein by MaxQuant according to the razor peptide rule with our implementation of the razor peptide rule"""
    score_type = ProteinScoringStrategy("multPEP razor")
    for peptide_row in psm.parse_evidence_file_single(mq_evidence_file, score_type = score_type):
        peptide, tmp_proteins, _, _ = peptide_row
        
        proteins = digest.get_proteins(peptide_to_protein_map, helpers.clean_peptide(peptide))
        
        leading_proteins = protein_groups.get_leading_proteins(proteins)
        if len(leading_proteins) > 1:
            predicted_razor = score_type.filter_proteins(proteins) # filtering for razor peptide approach
            logger.debug(f"{tmp_proteins[0] == predicted_razor[0]} {tmp_proteins[0]} {predicted_razor[0]}")
