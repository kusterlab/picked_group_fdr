# Picked Protein Group FDR

Scalable, accurate and sensitive protein group FDRs for large-scale mass spectrometry experiments

## Running Picked Protein Group FDR using the GUI

On Windows, you can download the `PickedGroupFDR_GUI_windows.zip` from the latest release, unzip it and open `PickedGroupFDR.exe` to start the GUI (no installation necessary).

Make sure that you have run the MaxQuant search with 100% protein-level FDR.

Alternatively, on all platforms, first install Picked Protein Group FDR as explained below. Then install `PyQt5` (`pip install PyQt5`) and run:

```shell
python gui.py
```

## Running Picked Protein Group FDR from the command line

1. install Picked Protein Group FDR as explained below.
2. make sure that you have run the MaxQuant search with 100% protein-level FDR.
3. the posterior error probabilities (PEP) of MaxQuant are not well-calibrated. Therefore, we first recalculate these with [Mokapot](https://mokapot.readthedocs.io/en/latest/) (=[Percolator](http://percolator.ms/) for Python):
   ```shell
   python3 -u -m picked_group_fdr.pipeline.andromeda2pin </path/to/mq_evidence_txt> --outputTab andromeda.tab --databases </path/to/fasta_file>
   python3 -u -m picked_group_fdr.pipeline.run_mokapot 0.01 0.01 percolator <num_threads>
   python3 -u -m picked_group_fdr.pipeline.update_evidence_from_pout --mq_evidence </path/to/mq_evidence_txt> --perc_results percolator/andromeda.mokapot.psms.txt percolator/andromeda.mokapot.decoys.psms.txt --mq_evidence_out percolator/evidence.txt
   ```
    Alternatively, you can use [Prosit](https://www.proteomicsdb.org/prosit/)'s Percolator results files directly:
   ```shell
   python3 -u -m picked_group_fdr.pipeline.update_evidence_from_pout --mq_evidence </path/to/mq_evidence_txt> --perc_results prosit_target.psms prosit_decoy.psms --mq_evidence_out percolator/evidence.txt --pout_input_type prosit
   ```
4. to obtain protein group level FDRs, run:
   ```shell
   python -m picked_group_fdr --mq_evidence percolator/evidence.txt --fasta </path/to/fasta_file>
   ```


## Installation

Picked Protein Group FDR is available on PyPI and can be installed with `pip`:

```shell
pip install picked_group_fdr
```

Alternatively, you can install directly from this repository:

```shell
git clone https://github.com/kusterlab/picked_group_fdr.git
pip install .
```


## Manuscript figures

To reproduce the figures from the manuscript:

1. make sure that you have `make` installed.
2. download the input files from zenodo: https://zenodo.org/record/7157677
3. specify the location of the input files and run the make command with the figure you want to reproduce, e.g.:
   ```shell
   export DATA_DIR=</path/to/zenodo/files>
   make Figure3a
   ```

Check the file `MakefileFigures` to see which figures are available.
