import sys
import logging
import os

import numpy as np

from .. import __version__, __copyright__
from ..parsers import parsers
from ..parsers import protein_groups as pgp
from .. import helpers
from .. import fdr
from ..scoring import BestAndromedaScore
from ..plotter import PlotterFactory
from ..picked_group_fdr import ArgumentParserWithLogger

# hacky way to get package logger when running as module
logger = logging.getLogger(__package__ + "." + __file__)


# python -m picked_group_fdr.pipeline.entrapment_fdr \
# --protein_groups_files proteinGroups.txt --plot_figures
def parse_args(argv):
    import argparse

    apars = ArgumentParserWithLogger(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    apars.add_argument(
        "--peptides_files",
        default=None,
        metavar="M",
        nargs="+",
        required=False,
        help="MaxQuant msms.txt or evidence.txt file(s).",
    )

    apars.add_argument(
        "--protein_groups_files",
        default=None,
        metavar="PG",
        nargs="+",
        required=False,
        help="proteinGroups.txt filtered at 100% FDR file(s)",
    )

    apars.add_argument(
        "--protein_col",
        default="Protein IDs",
        metavar="O",
        required=False,
        help="Column name with protein identifiers.",
    )

    apars.add_argument(
        "--peptide_col",
        default="Modified sequence",
        metavar="O",
        required=False,
        help="Column name for peptide sequences.",
    )

    apars.add_argument(
        "--score_col",
        default="Score",
        metavar="O",
        required=False,
        help="Column name with scores.",
    )

    apars.add_argument(
        "--plot_labels",
        default=None,
        metavar="PG",
        nargs="+",
        required=False,
        help="""Label(s) for the plots. Can be used to indicate which input files 
                belong to the same experiment. Must have the same length as the number
                of input files.""",
    )

    apars.add_argument(
        "--is_decoy_file",
        default=None,
        metavar="PG",
        nargs="+",
        required=False,
        help="""Boolean for each input file indicating if it is a decoy input or not.
                Must have the same length as the number of input files.""",
    )

    apars.add_argument(
        "--fdr_cutoff", default=0.01, metavar="C", type=float, help="FDR threshold."
    )

    apars.add_argument(
        "--absent_ratio",
        default=1.0,
        metavar="C",
        type=float,
        help="""Ratio of entrapment sequences to total target sequences. E.g. if 9
                entrapment databases were used, this ratio is 0.9. This ratio is used
                to correct the peptide level entrapment FDR. Do not apply this for 
                protein-level entrapment FDRs!""",
    )

    apars.add_argument(
        "--figure_base_fn",
        default=None,
        metavar="F",
        help="""Base file name for calibration and performance figures.""",
    )

    apars.add_argument(
        "--plot_figures", help="Plot figures with matplotlib.", action="store_true"
    )

    # ------------------------------------------------
    args = apars.parse_args(argv)

    return args


def main(argv):
    logger.info(f"EntrapmentFDR version {__version__}\n{__copyright__}")
    logger.info(
        f'Issued command: {os.path.basename(__file__)} {" ".join(map(str, argv))}'
    )

    args = parse_args(argv)

    if args.protein_groups_files:
        num_input_files = len(args.protein_groups_files)
    elif args.peptides_files:
        num_input_files = len(args.peptides_files)
    else:
        raise ValueError(
            "Either --protein_groups_files or --peptides_files needs to be set."
        )

    if args.plot_labels:
        if num_input_files != len(args.plot_labels):
            raise ValueError(
                "--protein_groups_files and --plot_labels need to have the same number of arguments"
            )
    else:
        args.plot_labels = list(range(num_input_files))

    if args.is_decoy_file:
        if num_input_files != len(args.is_decoy_file):
            raise ValueError(
                "--protein_groups_files and --is_decoy_file need to have the same number of arguments"
            )
        args.is_decoy_file = list(map(helpers.string_to_bool, args.is_decoy_file))
    else:
        args.is_decoy_file = [False] * num_input_files

    plotter = PlotterFactory.get_plotter(args.figure_base_fn, args.plot_figures)

    unique_labels = sorted(set(args.plot_labels))
    for label in unique_labels:
        logger.info(f"Processing label: {label}")
        if args.protein_groups_files:
            protein_group_files, is_decoy_file = zip(
                *[
                    (f, d)
                    for f, l, d in zip(
                        args.protein_groups_files, args.plot_labels, args.is_decoy_file
                    )
                    if l == label
                ]
            )
            protein_groups, protein_scores = zip(
                *list(
                    pgp.parse_protein_groups_file_multiple(
                        protein_group_files,
                        are_decoy_file=is_decoy_file,
                        protein_column=args.protein_col,
                        score_column=args.score_col,
                    )
                )
            )

            score_group_tuples = list(zip(protein_groups, protein_scores))
            np.random.shuffle(score_group_tuples)
            score_group_tuples = sorted(
                score_group_tuples, key=lambda x: x[1], reverse=True
            )
            protein_groups, protein_scores = zip(*score_group_tuples)

            reported_qvals, observed_qvals = fdr.calculate_protein_fdrs(
                protein_groups, protein_scores, args.fdr_cutoff
            )
        elif args.peptides_files:
            peptide_files, is_decoy_file = zip(
                *[
                    (f, d)
                    for f, l, d in zip(
                        args.peptides_files, args.plot_labels, args.is_decoy_file
                    )
                    if l == label
                ]
            )
            peptides, proteins, _, peptide_scores = zip(
                *list(
                    parsers.parse_peptides_files_multiple(
                        peptide_files,
                        are_decoy_file=is_decoy_file,
                        peptide_column=args.peptide_col,
                        protein_column=args.protein_col,
                        score_column=args.score_col,
                    )
                )
            )

            score_group_tuples = list(zip(peptide_scores, peptides, proteins))
            np.random.shuffle(score_group_tuples)
            score_group_tuples = sorted(
                score_group_tuples, key=lambda x: x[0], reverse=True
            )

            reported_qvals, observed_qvals = fdr.calculate_peptide_fdrs(
                score_group_tuples,
                score_type=BestAndromedaScore(),
                peptide_fdr_threshold=args.fdr_cutoff,
            )

        plotter.label = label
        plotter.plot_qval_calibration_and_performance(
            reported_qvals, observed_qvals, absent_ratio=args.absent_ratio
        )

    plotter.decorate_plots()
    plotter.save_plots()
    plotter.show()


if __name__ == "__main__":
    main(sys.argv[1:])
