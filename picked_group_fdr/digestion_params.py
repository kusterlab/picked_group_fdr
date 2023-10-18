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
        enzyme,
        digestion,
        min_length,
        max_length,
        cleavages,
        special_aas,
        fasta_contains_decoys,
    ):
        self.enzyme = enzyme
        self.digestion = digestion
        self.min_length = min_length
        self.max_length = max_length
        self.cleavages = cleavages
        self.special_aas = list(special_aas)

        self.methionine_cleavage = True
        self.db = "target" if fasta_contains_decoys else "concat"
        if self.enzyme == "no_enzyme":
            self.digestion = "none"

        self.use_hash_key = self.digestion == "none"


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
        [args.fasta_contains_decoys]
    ]
    param_lengths = [len(p) for p in params_list if len(p) != 1]
    if len(set(param_lengths)) > 1:
        raise ValueError("Received digestion parameters of unequal length.")

    max_params = max(param_lengths) if len(param_lengths) > 0 else 1
    params_list_updated = []
    for param in params_list:
        if len(param) == 1:
            param = param*max_params
        params_list_updated.append(param)

    return [DigestionParams(*p) for p in zip(*params_list_updated)]


def add_digestion_arguments(apars):
    apars.add_argument('-e', '--enzyme', default = [ENZYME_DEFAULT], metavar='E', nargs="+",
                                         help='''Type of enzyme "no_enzyme","elastase","pepsin",
                                                         "proteinasek","thermolysin","chymotrypsin",
                                                         "lys-n","lys-c","arg-c","asp-n","glu-c","trypsin",
                                                         "trypsinp".
                                                    ''')

    apars.add_argument('-c', '--cleavages', default = [CLEAVAGES_DEFAULT], metavar='C', type=int, nargs="+",
                                         help='''Number of allowed miss cleavages used in the search 
                                                         engine (Only valid when using option -F).
                                                    ''')

    apars.add_argument('-l', '--min-length', default = [MIN_PEPLEN_DEFAULT], metavar='L', type=int, nargs="+",
                                         help='''Minimum peptide length allowed used in the search 
                                                         engine (Only valid when using option -F).
                                                    ''')

    apars.add_argument('-t', '--max-length', default = [MAX_PEPLEN_DEFAULT], metavar='L', type=int, nargs="+",
                                         help='''Maximum peptide length allowed used in the search 
                                                         engine (Only valid when using option -F).
                                                    ''')

    apars.add_argument('--special-aas', default = [SPECIAL_AAS_DEFAULT], metavar='S', nargs="+",
                                         help='''Special AAs that MaxQuant uses for decoy generation.
                                                    ''')

    apars.add_argument('--digestion', default = [DIGESTION_DEFAULT], metavar='D', nargs="+",
                                         help='''Digestion mode ('full', 'semi' or 'none').
                                                    ''')

    apars.add_argument('--fasta_contains_decoys',
                         help='Set this flag if your fasta file already contains decoy protein sequences.',
                         action='store_true')
