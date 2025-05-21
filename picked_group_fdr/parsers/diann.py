from __future__ import annotations

import logging

import pandas as pd

from . import tsv

# for type hints only
from .. import scoring_strategy

logger = logging.getLogger(__name__)

"""report.parquet columns:
https://github.com/vdemichev/DiaNN?tab=readme-ov-file#main-output-reference

1 Run.Index          # experiment index, 0-based
2 Run                # experiment name
3 Channel            # empty for regular LFQ
4 Precursor.Id       # Modified.Sequence + Precursor.Charge
5 Modified.Sequence  # example: ^YVLYC(UniMod:4)AVISK
6 Stripped.Sequence  # example: YVLYCAVISK (no modifications)  
7 Precursor.Charge  
8 Precursor.Lib.Index  
9 Decoy              # 0 if target, 1 if decoy
10 Proteotypic  
11 Precursor.Mz  
12 Protein.Ids       # does not contain decoy prefix. example: C8Z8S1;C8ZDU4
13 Protein.Group     # does not contain decoy prefix. example: C8Z8S1;C8ZDU4
14 Protein.Names     # does not contain decoy prefix. example: C8Z8S1_YEAS8;C8ZDU4_YEAS8
15 Genes             # does not contain decoy prefix. example: EC1118_1G1_3290g;EC1118_1L7_1981g
16 RT  
17 iRT  
18 Predicted.RT  
19 Predicted.iRT  
20 IM  
21 iIM  
22 Predicted.IM  
23 Predicted.iIM  
24 Precursor.Quantity    # non-normalized MS2 quant, see https://github.com/vdemichev/DiaNN/discussions/951
25 Precursor.Normalised  # normalized MS2 quant
26 Ms1.Area              # non-normalized MS1 quant
27 Ms1.Normalised        # normalized MS1 quant
28 Ms1.Apex.Area  
29 Ms1.Apex.Mz.Delta  
30 Normalisation.Factor  
31 Quantity.Quality  
32 Empirical.Quality  
33 Normalisation.Noise  
34 Ms1.Profile.Corr  
35 Evidence  
36 Mass.Evidence  
37 Channel.Evidence  
38 Ms1.Total.Signal.Before  
39 Ms1.Total.Signal.After  
40 RT.Start  
41 RT.Stop  
42 FWHM  
43 PG.TopN  
44 PG.MaxLFQ  
45 Genes.TopN  
46 Genes.MaxLFQ  
47 Genes.MaxLFQ.Unique  
48 PG.MaxLFQ.Quality  
49 Genes.MaxLFQ.Quality  
50 Genes.MaxLFQ.Unique.Quality  
51 Q.Value  
52 PEP  
53 Global.Q.Value  
54 Lib.Q.Value  
55 Peptidoform.Q.Value  
56 Global.Peptidoform.Q.Value  
57 Lib.Peptidoform.Q.Value  
58 PTM.Site.Confidence  
59 Site.Occupancy.Probabilities  
60 Protein.Sites  
61 Lib.PTM.Site.Confidence  
62 Translated.Q.Value  
63 Channel.Q.Value  
64 PG.Q.Value  
65 PG.PEP  
66 GG.Q.Value  
67 Protein.Q.Value  
68 Global.PG.Q.Value  
69 Lib.PG.Q.Value  
70 Best.Fr.Mz  
71 Best.Fr.Mz.Delta  
"""


def parse_diann_report_file(
    evidence_file: str,
    get_proteins,
    score_column: str,
    for_quantification: bool = False,
    **kwargs,
):
    if evidence_file.endswith(".parquet"):
        df = pd.read_parquet(evidence_file)
    else:
        delimiter = tsv.get_delimiter(evidence_file)
        df = pd.read_csv(evidence_file, sep=delimiter)

    df.columns = df.columns.str.replace(".", "_")

    if "Decoy" not in df.columns:
        raise ValueError(
            "DIA-NN report file does not contain a 'Decoy' column. Make sure to use DIA-NN v2.0 or higher and add the --report-decoys flag as additional parameter."
        )

    logger.info("Parsing DIA-NN report file")
    for line_idx, row in enumerate(df.itertuples()):
        if line_idx % 500000 == 0:
            logger.info(f"    Reading line {line_idx}")

        score = row.PEP
        peptide = row.Modified_Sequence
        proteins = row.Protein_Ids.split(";")
        if row.Decoy == 1:
            proteins = [f"REV__{protein}" for protein in proteins]
        proteins = get_proteins(peptide, proteins)
        experiment = row.Run
        if not proteins:
            continue

        if not for_quantification:
            yield peptide, proteins, experiment, score
        else:
            tmt_intensities = []
            silac_intensities = []  # TODO: add SILAC support
            evidence_id = -1
            raw_file = row.Run
            fraction = -1  # TODO: add fractionation support
            yield (
                peptide,
                proteins,
                row.Precursor_Charge,
                raw_file,
                experiment,
                fraction,
                row.Ms1_Normalised,  # we use MS1 quant because MS2 quants are not comparable across multiple runs because different fragment ions can be selected.
                score,
                tmt_intensities,
                silac_intensities,
                evidence_id,
            )
