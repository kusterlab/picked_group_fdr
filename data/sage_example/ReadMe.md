# Sage example (LFQ)

The integration of PickedGroupFDR with Sage is realized by taking PSM level identification information from `results.sage.tsv` and combining it with the quantification information in `lfq.tsv`.

The data consists of 9 RAW files with different pools of PrESTs spiked into an E. coli background.
The dataset is described in this publication: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6474350/

To run this example, follow these steps:

## Option 1: reuse previous Sage results

**Download Sage result files**

1. Download the `sage.zip` file from this folder. It contains the necessary Sage results files needed to run this example. For the settings used to generate these files, see step 3 in Option 2.

**Running PickedGroupFDR**

2. Run the `run_picked_group_fdr.sh` script. This runs the `picked_group_fdr` and `picked_group_fdr.pipeline.sage_quantification` modules.

A quantified protein group result in FragPipe format will be available at `results/combined_protein.tsv`.

As an additional file, it will produce a protein grouping result similar to MaxQuant's proteinGroups.txt format at `results/proteinGroups.txt` without quantification information.


## Option 2: run Sage yourself

**Download RAW files and convert**

1. Download the `mixture{A,B,AB}rep{1,2,3}.raw` files from the PRIDE repository [PXD008425](https://www.ebi.ac.uk/pride/archive/projects/PXD008425) into a new subfolder named `RAW`.
2. Convert the raw files into `mzML` format into a new folder named `mzML`, e.g. with [msconvert](https://proteowizard.sourceforge.io/download.html) or [ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser)

**Running Sage**

3. Run sage using the `config.json` in this folder and the fasta file `iprg2016_with_labels.fasta` included in the `sage.zip` file.
   ```
   sage --fasta iprg2016_with_labels.fasta --output_directory ./ config.json mzML/*
   ```

**Running PickedGroupFDR**

4. Run the `run_picked_group_fdr.sh` script. This runs the `picked_group_fdr` and `picked_group_fdr.pipeline.sage_quantification` modules.

A quantified protein group result in FragPipe format will be available at `results/combined_protein.tsv`.

As an additional file, it will produce a protein grouping result similar to MaxQuant's proteinGroups.txt format at `results/proteinGroups.txt` without quantification information.
