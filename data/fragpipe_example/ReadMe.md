# FragPipe example (LFQ)

The integration of PickedGroupFDR with FragPipe is realized by updating the `psm.tsv` and `protein.tsv` files in each of the experiment folders of FragPipe. Quantifications of the protein groups then takes place using the [IonQuant](https://ionquant.nesvilab.org/) software.

The data consists of 9 RAW files with different pools of PrESTs spiked into an E. coli background.
The dataset is described in this publication: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6474350/

The `fragpipe.zip` file contains the necessary FragPipe results files needed to run this example.

To run this example, follow these steps:

**Download RAW files and FragPipe**

1. Download the `mixture{A,B,AB}rep{1,2,3}.raw` files from the PRIDE repository [PXD008425](https://www.ebi.ac.uk/pride/archive/projects/PXD008425) into a new subfolder named `RAW`.
2. Download and install FragPipe (version used for this example was v19.1).

**Running FragPipe (optional)**

The zip file already contains all necessary FragPipe results files needed to run this example. 
If you want to rerun the analysis with different parameters, follow these steps:

3. Use the fasta file `iprg2016_with_labels.fasta` included in the `fragpipe.zip` file.
4. In the `Validation` tab, in the `FDR Filter and Report` section, tick the `Print decoys` checkbox and add the following in the `Filter` field.
   ```
   --sequential --prot 1.0 --ion 1.0 --pept 1.0 --psm 1.0
   ```
   This sets the FDR to 100% on all levels, which allows PickedGroupFDR to achieve increased sensitivity over the regular FragPipe results.

**Running PickedGroupFDR**

5. Run the `run_picked_group_fdr.sh` script.

This will produce a protein grouping result similar to MaxQuant's proteinGroups.txt format at `results/proteinGroups.txt`.

**Running IonQuant**

6. Run the `run_ionquant.sh` script with the path to the FragPipe software directory as the first argument:
   ```
   ./run_ionquant.sh "/path/to/FragPipe-19.1"
   ```

The quantified protein group results will now be available at `results/combined_protein.tsv`.