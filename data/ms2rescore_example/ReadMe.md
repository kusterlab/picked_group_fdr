# MS2Rescore example

The integration of PickedGroupFDR with MS2Rescore relies on its percolator/mokapot output format.

The data consists of 9 RAW files with different pools of PrESTs spiked into an E. coli background.
The dataset is described in this publication: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6474350/

To run this example, follow these steps:

## Option 1: reuse previous MS2Rescore results

**Download MS2Rescore result files**

1. Download and unzip the `ms2rescore.zip` file from this folder. It contains the necessary MS2Rescore results files for a Sage search needed to run this example. For the settings used to generate these files, see the `./data/sage_example` ReadMe and the example below.

**Running PickedGroupFDR**

2. Run the following command:
   ```
   python3 -u -m picked_group_fdr \
      --fasta ./sage/iprg2016_with_labels.fasta \
      --perc_evidence ./sage/results.sage.ms2rescore.mokapot.psms.txt ./sage/results.sage.ms2rescore.mokapot.decoy.psms.txt
      --protein_groups_out ./results/combined_protein.tsv \
      --output_format fragpipe \
      --methods picked_protein_group_no_remap
   ```

A quantified protein group result in FragPipe format will be available at `results/combined_protein.tsv`. It is also possible to generate an output file in MaxQuant's `proteinGroups.txt` format using the `--output_format maxquant` flag.

## Option 2: run MS2Rescore yourself

**Download Sage result files**

1. Download and unzip the `sage.zip` file from the `./data/sage_example` folder. It contains the necessary Sage results files needed to run this example. For the settings used to generate these files, see the `./data/sage_example` ReadMe.

**Download RAW files and convert**

2. Download the `mixture{A,B,AB}rep{1,2,3}.raw` files from the PRIDE repository [PXD008425](https://www.ebi.ac.uk/pride/archive/projects/PXD008425) into a new subfolder named `RAW`.
3. Convert the raw files into `mzML` format into a new folder named `mzML`, e.g. with [msconvert](https://proteowizard.sourceforge.io/download.html) or [ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser)

**Run MS2rescore**

4. Make sure MS2Rescore is installed and run the following command:
   ```
   ms2rescore --psm-file sage/results.sage.tsv \
           --psm-file-type 'sage' \
           --spectrum-path mzML/ \
           -f iprg2016_with_labels.fasta \
           -n 4 \
           -o sage/results.sage.ms2rescore
   ```

**Run PickedGroupFDR**

5. Run the following command:
   ```
   python3 -u -m picked_group_fdr \
      --fasta ./sage/iprg2016_with_labels.fasta \
      --perc_evidence ./sage/results.sage.ms2rescore.mokapot.psms.txt ./sage/results.sage.ms2rescore.mokapot.decoy.psms.txt
      --protein_groups_out ./results/combined_protein.tsv \
      --output_format fragpipe \
      --methods picked_protein_group_no_remap
   ```

A quantified protein group result in FragPipe format will be available at `results/combined_protein.tsv`. It is also possible to generate an output file in MaxQuant's `proteinGroups.txt` format using the `--output_format maxquant` flag.