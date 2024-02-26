from abc import ABC, abstractmethod

from .parsers import fragpipe
from .parsers import maxquant
from .parsers import percolator
from .parsers import sage
from .quant import fragpipe as fragpipe_quant
from .quant import maxquant as mq_quant
from .quant import sage as sage_quant


class ScoreOrigin(ABC):
    @abstractmethod
    def get_evidence_file(self, args):
        pass

    @abstractmethod
    def get_evidence_parser(self):
        pass

    @abstractmethod
    def get_quantification_file(self, args):
        pass

    @abstractmethod
    def get_quantification_parser(self):
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

    def get_quantification_file(self, args):
        raise NotImplementedError("Cannot perform quantification using Percolator input.")

    def get_quantification_parser(self):
        raise NotImplementedError("Cannot perform quantification using Percolator input.")

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

    def get_quantification_file(self, args):
        return args.mq_evidence

    def get_quantification_parser(self):
        return mq_quant.add_precursor_quants

    def remaps_peptides_to_proteins(self):
        return True

    def can_do_quantification(self):
        return True

    def short_description(self):
        return 'm'

    def long_description(self):
        return 'MaxQuant'


class MaxQuantInputNoRemap(MaxQuantInput):
    def remaps_peptides_to_proteins(self):
        return False


class FragPipeInput(ScoreOrigin):
    def get_evidence_file(self, args):
        return args.fragpipe_psm

    def get_evidence_parser(self):
        return fragpipe.parse_fragpipe_psm_file

    def get_quantification_file(self, args):
        return args.combined_ion

    def get_quantification_parser(self):
        return fragpipe_quant.add_precursor_quants_multiple

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

    def get_quantification_file(self, args):
        return args.sage_lfq_tsv

    def get_quantification_parser(self):
        return sage_quant.add_precursor_quants_multiple

    def remaps_peptides_to_proteins(self):
        return False

    def can_do_quantification(self):
        return True

    def short_description(self):
        return 's'

    def long_description(self):
        return 'Sage'