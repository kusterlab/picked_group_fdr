from typing import Protocol, Callable
import argparse

from .parsers import fragpipe
from .parsers import maxquant
from .parsers import percolator
from .parsers import sage
from .parsers import diann

from .quant import fragpipe as fragpipe_quant
from .quant import maxquant as mq_quant
from .quant import sage as sage_quant


class ScoreOrigin(Protocol):
    def get_evidence_file(self, args: argparse.Namespace) -> str: ...
    def get_evidence_parser(self) -> Callable: ...
    def get_quantification_file(self, args: argparse.Namespace) -> str: ...
    def get_quantification_parser(self) -> Callable: ...
    def remaps_peptides_to_proteins(self) -> bool: ...
    def can_do_quantification(self) -> bool: ...
    def short_description(self) -> str: ...
    def long_description(self) -> str: ...


class PercolatorInput:
    def get_evidence_file(self, args):
        return args.perc_evidence

    def get_evidence_parser(self):
        return percolator.parse_percolator_out_file

    def get_quantification_file(self, args):
        raise NotImplementedError(
            "Cannot perform quantification using Percolator input."
        )

    def get_quantification_parser(self):
        raise NotImplementedError(
            "Cannot perform quantification using Percolator input."
        )

    def remaps_peptides_to_proteins(self):
        return False

    def can_do_quantification(self):
        return False

    def short_description(self):
        return "p"

    def long_description(self):
        return "Percolator"


class PercolatorInputRemapped(PercolatorInput):
    def remaps_peptides_to_proteins(self):
        return True


class MaxQuantInput:
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
        return "m"

    def long_description(self):
        return "MaxQuant"


class MaxQuantInputNoRemap(MaxQuantInput):
    def remaps_peptides_to_proteins(self):
        return False


class FragPipeInput:
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
        return "f"

    def long_description(self):
        return "FragPipe"


class SageInput:
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
        return "s"

    def long_description(self):
        return "Sage"


class DiannInput:
    def get_evidence_file(self, args):
        return args.diann_reports

    def get_evidence_parser(self):
        return diann.parse_diann_report_file

    def get_quantification_file(self, args):
        return args.diann_reports

    def get_quantification_parser(self):
        return mq_quant.add_precursor_quants

    def remaps_peptides_to_proteins(self):
        return False

    def can_do_quantification(self):
        return True

    def short_description(self):
        return "d"

    def long_description(self):
        return "DIA-NN"
