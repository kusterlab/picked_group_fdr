import pytest
import argparse
from unittest.mock import patch

from picked_group_fdr import picked_group_fdr
from picked_group_fdr.results import ProteinGroupResult


def test_run_picked_group_fdr():
    # Mocking the required objects and functions
    args = argparse.Namespace(
        fasta="example.fasta",
        fasta_contains_decoys=False,
        gene_level=True,
        methods="method1,method2",
        figure_base_fn="figures",
        plot_figures=True,
    )
    method_configs = [{"name": "method1"}, {"name": "method2"}]
    protein_annotations = {"annotation": "example"}
    use_pseudo_genes = False
    peptide_to_protein_maps = [None]

    with patch(
        "picked_group_fdr.protein_annotation.get_protein_annotations",
        return_value=(protein_annotations, use_pseudo_genes),
    ) as mock_get_annotations, patch(
        "picked_group_fdr.methods.get_methods", return_value=method_configs
    ) as mock_get_methods, patch(
        "picked_group_fdr.methods.requires_peptide_to_protein_map", return_value=True
    ) as mock_requires_ptp_map, patch(
        "picked_group_fdr.peptide_protein_map.get_peptide_to_protein_maps_from_args",
        return_value=peptide_to_protein_maps,
    ) as mock_get_ptp_map, patch(
        "picked_group_fdr.picked_group_fdr.run_method"
    ) as mock_run_method, patch(
        "picked_group_fdr.plotter.PlotterFactory.get_plotter"
    ) as mock_get_plotter, patch(
        "picked_group_fdr.picked_group_fdr.logger"
    ) as mock_logger:

        # Call the function
        picked_group_fdr.run_picked_group_fdr(args)

        # Assert that the necessary functions were called with the correct arguments
        mock_get_annotations.assert_called_once_with("example.fasta", False, True)
        mock_get_methods.assert_called_once_with("method1,method2", use_pseudo_genes)
        mock_requires_ptp_map.assert_called_once_with(method_configs)
        mock_get_ptp_map.assert_called_once_with(args, use_pseudo_genes)
        mock_get_plotter.assert_called_once_with("figures", True)

        # Assert that run_method was called for each method configuration
        assert mock_run_method.call_count == 2
        mock_run_method.assert_any_call(
            args,
            method_configs[0],
            peptide_to_protein_maps,
            protein_annotations,
            mock_get_plotter.return_value,
            use_pseudo_genes,
            apply_filename_suffix=True,
        )
        mock_run_method.assert_any_call(
            args,
            method_configs[1],
            peptide_to_protein_maps,
            protein_annotations,
            mock_get_plotter.return_value,
            use_pseudo_genes,
            apply_filename_suffix=True,
        )

        # Assert that logger.info was called with the correct message
        mock_logger.info.assert_called_once_with(
            "PickedGroupFDR execution took 0.0 seconds wall clock time"
        )


class TestGetProteinGroupResults:
    def test_get_protein_group_results_with_defaults(self):
        results = picked_group_fdr.get_protein_group_results(
            {"PEPA": (1e-5, ["PROTA", "PROTB"])}
        )

        assert results.protein_group_results == [
            ProteinGroupResult(
                proteinIds="PROTB;PROTA",
                majorityProteinIds="PROTB;PROTA",
                peptideCountsUnique="1;1",
                bestPeptide="PEPA",
                numberOfProteins=2,
                qValue=0.5,
                score=5.0,
                reverse="",
                potentialContaminant="",
                precursorQuants=[],
                extraColumns=[],
            )
        ]


@pytest.fixture
def arg_parser():
    """Fixture to create ArgumentParser instance."""
    from argparse import ArgumentParser

    return ArgumentParser()


class TestParseArgs:
    def test_parse_args_with_mq_evidence(self, arg_parser):
        """Test parse_args with mq_evidence argument."""
        argv = ["--mq_evidence", "file1.txt"]
        args = picked_group_fdr.parse_args(argv)
        assert args.mq_evidence == ["file1.txt"]

    def test_parse_args_with_perc_evidence(self, arg_parser):
        """Test parse_args with perc_evidence argument."""
        argv = ["--perc_evidence", "file2.txt"]
        args = picked_group_fdr.parse_args(argv)
        assert args.perc_evidence == ["file2.txt"]

    def test_parse_args_with_fragpipe_psm(self, arg_parser):
        """Test parse_args with fragpipe_psm argument."""
        argv = ["--fragpipe_psm", "file3.txt"]
        args = picked_group_fdr.parse_args(argv)
        assert args.fragpipe_psm == ["file3.txt"]

    def test_parse_args_with_combined_ion(self, arg_parser):
        """Test parse_args with combined_ion argument."""
        argv = ["--combined_ion", "file4.txt"]
        args = picked_group_fdr.parse_args(argv)
        assert args.combined_ion == ["file4.txt"]

    def test_parse_args_with_sage_results(self, arg_parser):
        """Test parse_args with sage_results argument."""
        argv = ["--sage_results", "file5.txt"]
        args = picked_group_fdr.parse_args(argv)
        assert args.sage_results == ["file5.txt"]

    def test_parse_args_with_sage_lfq_tsv(self, arg_parser):
        """Test parse_args with sage_lfq_tsv argument."""
        argv = ["--sage_lfq_tsv", "file6.txt"]
        args = picked_group_fdr.parse_args(argv)
        assert args.sage_lfq_tsv == ["file6.txt"]

    def test_parse_args_with_protein_groups_out(self, arg_parser):
        """Test parse_args with protein_groups_out argument."""
        argv = ["--protein_groups_out", "output.txt"]
        args = picked_group_fdr.parse_args(argv)
        assert args.protein_groups_out == "output.txt"

    def test_parse_args_with_output_format(self, arg_parser):
        """Test parse_args with output_format argument."""
        argv = ["--output_format", "fragpipe"]
        args = picked_group_fdr.parse_args(argv)
        assert args.output_format == "fragpipe"

    def test_parse_args_with_fasta(self, arg_parser):
        """Test parse_args with fasta argument."""
        argv = ["--fasta", "fasta_file.txt"]
        args = picked_group_fdr.parse_args(argv)
        assert args.fasta == ["fasta_file.txt"]

    def test_parse_args_with_mq_protein_groups(self, arg_parser):
        """Test parse_args with mq_protein_groups argument."""
        argv = ["--mq_protein_groups", "mq_protein_groups.txt"]
        args = picked_group_fdr.parse_args(argv)
        assert args.mq_protein_groups == "mq_protein_groups.txt"

    def test_parse_args_with_methods(self, arg_parser):
        """Test parse_args with methods argument."""
        argv = ["--methods", "method1"]
        args = picked_group_fdr.parse_args(argv)
        assert args.methods == "method1"

    def test_parse_args_with_peptide_protein_map(self, arg_parser):
        """Test parse_args with peptide_protein_map argument."""
        argv = ["--peptide_protein_map", "map_file.txt"]
        args = picked_group_fdr.parse_args(argv)
        assert args.peptide_protein_map == ["map_file.txt"]

    def test_parse_args_with_keep_all_proteins(self, arg_parser):
        """Test parse_args with keep_all_proteins argument."""
        argv = ["--keep_all_proteins"]
        args = picked_group_fdr.parse_args(argv)
        assert args.keep_all_proteins

    def test_parse_args_with_gene_level(self, arg_parser):
        """Test parse_args with gene_level argument."""
        argv = ["--gene_level"]
        args = picked_group_fdr.parse_args(argv)
        assert args.gene_level

    def test_parse_args_with_do_quant(self, arg_parser):
        """Test parse_args with do_quant argument."""
        argv = ["--do_quant"]
        args = picked_group_fdr.parse_args(argv)
        assert args.do_quant

    def test_parse_args_with_suppress_missing_peptide_warning(self, arg_parser):
        """Test parse_args with suppress_missing_peptide_warning argument."""
        argv = ["--suppress_missing_peptide_warning"]
        args = picked_group_fdr.parse_args(argv)
        assert args.suppress_missing_peptide_warning

    def test_parse_args_with_figure_base_fn(self, arg_parser):
        """Test parse_args with figure_base_fn argument."""
        argv = ["--figure_base_fn", "figure_base"]
        args = picked_group_fdr.parse_args(argv)
        assert args.figure_base_fn == "figure_base"

    def test_parse_args_with_plot_figures(self, arg_parser):
        """Test parse_args with plot_figures argument."""
        argv = ["--plot_figures"]
        args = picked_group_fdr.parse_args(argv)
        assert args.plot_figures
