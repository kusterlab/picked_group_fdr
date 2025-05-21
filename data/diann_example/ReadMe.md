# DIA-NN example (LFQ)

This uses a subset of the LFQbench DIA dataset (PXD002952) as used in the Triqler-DIA paper (https://pubs.acs.org/doi/full/10.1021/acs.jproteome.2c00607):

```
We downloaded the LFQBench data set (6) from ProteomeXchange (PXD002952). Here, we used the TripleTOF6600 section of the study, which was harvested with a setup of 32 fixed windows MS2-windows. We also restricted ourselves to the low ratio difference samples, referred to as the HYE124 hybrid proteome samples in the original study. These consist of triplicates of sample A, composed of tryptically digested proteins from 65% w/w HeLa, 30% w/w yeast, and 5% w/w E. coli cells, and triplicates of sample B, composed of 65% w/w, 15% w/w yeast, and 20% w/w E. coli proteins. Samples from HYE110 and the TripleTOF5600 section of PXD002952 were omitted in this study.
```

To run this example, follow these steps:

## Option 1: reuse previous Sage results

**Download DIA-NN result files**

1. Download and unzip the [`diann.zip`](https://zenodo.org/records/15228406/files/diann.zip?download=1) file in this folder. It contains the necessary DIA-NN results files needed to run this example. For the settings used to generate these files, see step 3 in Option 2.
   ```
   wget -O diann.zip https://zenodo.org/records/15228406/files/diann.zip?download=1
   unzip diann.zip
   ```

**Running PickedGroupFDR**

2. Run the following command:
   ```
   python3 -u -m picked_group_fdr \
      --fasta ./diann/human_yeast_ecoli.fasta \
      --fasta_use_uniprot_id \
      --diann_reports ./diann/report.parquet \
      --protein_groups_out ./results/report.pg_matrix.tsv \
      --do_quant \
      --lfq_min_peptide_ratios 1 \
      --methods diann
   ```

A quantified protein group result in DIA-NN protein groups format will be available at `results/report.pg_matrix.tsv`. It is also possible to generate an output file in MaxQuant's `proteinGroups.txt` format or FragPipe's `combined_protein.tsv` format using the `--output_format` flag.


## Option 2: run DIA-NN yourself

**Download RAW files and convert**

1. Download the `HYE124_TTOF6600_32fix_lgillet_I150211_*.{wiff,wiff.scan}` files from the PRIDE repository [PXD002952](https://www.ebi.ac.uk/pride/archive/projects/PXD002952) into a new subfolder named `RAW`.

**Running DIA-NN**

2. Open the GUI for DIA-NN v2.x on Windows.
3. Add the downloaded .wiff files as Raw inputs.
4. Click the `Add FASTA` button and add the fasta file `human_yeast_ecoli.fasta` (included in the `diann.zip` file).
5. Tick the checkbox `FASTA digest for library-free search / library generation`
6. Create a new folder named `diann` in this folder.
7. Set the path for the Main output file and Output library to `diann/report.parquet` and `diann/report-lib.parquet`.
8. Add `--report-decoys` in the Additional options field.
9. Click the `Run` button. Using 6 cores, the analysis took about 9 hours.

**Running PickedGroupFDR**

10. Run the following command:
   ```
   python3 -u -m picked_group_fdr \
      --fasta ./diann/human_yeast_ecoli.fasta \
      --fasta_use_uniprot_id \
      --diann_reports ./diann/report.parquet \
      --protein_groups_out ./results/report.pg_matrix.tsv \
      --do_quant \
      --lfq_min_peptide_ratios 1 \
      --methods diann
   ```

A quantified protein group result in DIA-NN protein groups format will be available at `results/report.pg_matrix.tsv`. It is also possible to generate an output file in MaxQuant's `proteinGroups.txt` format or FragPipe's `combined_protein.tsv` format using the `--output_format` flag.
