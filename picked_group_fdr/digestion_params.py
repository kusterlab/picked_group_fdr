from argparse import Namespace
from typing import List


ENZYME_DEFAULT = "trypsin"
CLEAVAGES_DEFAULT = 2
MIN_PEPLEN_DEFAULT = 7
MAX_PEPLEN_DEFAULT = 60
SPECIAL_AAS_DEFAULT = "KR"
DIGESTION_DEFAULT = "full"


class DigestionParams:
    enzyme: str
    digestion: str
    min_length: int
    max_length: int
    cleavages: int
    special_aas: str
    methionine_cleavage: bool
    db: str
    use_hash_key: bool

    def __init__(
        self,
        enzyme = ENZYME_DEFAULT,
        digestion = DIGESTION_DEFAULT,
        min_length = MIN_PEPLEN_DEFAULT,
        max_length = MAX_PEPLEN_DEFAULT,
        cleavages = CLEAVAGES_DEFAULT,
        special_aas = SPECIAL_AAS_DEFAULT,
        fasta_contains_decoys = False,
    ):
        self.enzyme = enzyme

        self.digestion = digestion
        if self.enzyme == "no_enzyme":
            self.digestion = "none"

        self.min_length = min_length
        self.max_length = max_length
        self.cleavages = cleavages

        self.special_aas = list()
        if special_aas != "none":
            self.special_aas = list(special_aas)

        self.methionine_cleavage = True
        self.db = "target" if fasta_contains_decoys else "concat"

        self.use_hash_key = self.digestion == "none"


def digestion_params_list_to_arg_list(
    digestion_params_list: List[DigestionParams],
) -> List[str]:
    return (
        ["--min-length"]
        + [str(p.min_length) for p in digestion_params_list]
        + ["--max-length"]
        + [str(p.max_length) for p in digestion_params_list]
        + ["--cleavages"]
        + [str(p.cleavages) for p in digestion_params_list]
        + ["--enzyme"]
        + [p.enzyme for p in digestion_params_list]
        + ["--digestion"]
        + [p.digestion for p in digestion_params_list]
        + ["--special-aas"]
        + ["".join(p.special_aas) for p in digestion_params_list]
    )


def get_digestion_params_list(args: Namespace) -> List[DigestionParams]:
    """Takes the parsed arguments from argparse and returns a list of DigestionParams.

    Args:
        args (Namespace): arguments from argparse.ArgumentParser.

    Raises:
        ValueError: if digestion parameters of length > 1 are of unequal length.

    Returns:
        List[DigestionParams]: list of DigestionParams with length of longest digestion parameter argument.
    """
    params_list = [
        args.enzyme,
        args.digestion,
        args.min_length,
        args.max_length,
        args.cleavages,
        args.special_aas,
        [args.fasta_contains_decoys],
    ]
    param_lengths = [len(p) for p in params_list if len(p) != 1]
    if len(set(param_lengths)) > 1:
        raise ValueError("Received digestion parameters of unequal length.")

    max_params = max(param_lengths) if len(param_lengths) > 0 else 1
    params_list_updated = []
    for param in params_list:
        if len(param) == 1:
            param = param * max_params
        params_list_updated.append(param)

    return [DigestionParams(*p) for p in zip(*params_list_updated)]


def add_digestion_arguments(apars):
    apars.add_argument(
        "-e",
        "--enzyme",
        default=[ENZYME_DEFAULT],
        metavar="E",
        nargs="+",
        help="""Enzyme used for digestion. Available enzymes are 
                "trypsin","trypsinp","no_enzyme","elastase","pepsin", 
                "proteinasek","thermolysin","chymotrypsin","chymotrypsin+",
                "lys-n","lys-c","lys-cp","arg-c","asp-n","glu-c".""",
    )

    apars.add_argument(
        "-c",
        "--cleavages",
        default=[CLEAVAGES_DEFAULT],
        metavar="C",
        type=int,
        nargs="+",
        help="""Number of allowed miss cleavages used in the search engine.""",
    )

    apars.add_argument(
        "-l",
        "--min-length",
        default=[MIN_PEPLEN_DEFAULT],
        metavar="L",
        type=int,
        nargs="+",
        help="""Minimum peptide length allowed used in the search engine.""",
    )

    apars.add_argument(
        "-t",
        "--max-length",
        default=[MAX_PEPLEN_DEFAULT],
        metavar="L",
        type=int,
        nargs="+",
        help="""Maximum peptide length allowed used in the search engine.""",
    )

    apars.add_argument(
        "--special-aas",
        default=[SPECIAL_AAS_DEFAULT],
        metavar="S",
        nargs="+",
        help="""Special AAs that MaxQuant uses for decoy generation. 
                Amino acids are written as a single string with all 
                amino acids, e.g. "RK". To specify no amino acids, 
                supply the string "none".""",
    )

    apars.add_argument(
        "--digestion",
        default=[DIGESTION_DEFAULT],
        metavar="D",
        nargs="+",
        help="""Digestion mode ('full', 'semi' or 'none').
                                                    """,
    )

    apars.add_argument(
        "--fasta_contains_decoys",
        help="Set this flag if your fasta file already contains decoy protein sequences.",
        action="store_true",
    )

    apars.add_argument(
        "--fasta_use_uniprot_id",
        help="""Parse protein identifiers in the fasta file as UniProt IDs, 
                i.e. Q9UM47 for the protein identifier sp|Q9UM47|NOTC3_HUMAN""",
        action="store_true",
    )
