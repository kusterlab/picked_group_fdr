# Picked Protein Group FDR

Scalable, accurate and sensitive protein group FDRs for large-scale mass spectrometry experiments

## Running Picked Protein Group FDR using the GUI

On Windows, you can download the `PickedGroupFDR_GUI_windows.zip` from the latest release, unzip it and open `PickedGroupFDR.exe` to start the GUI (no installation necessary).

Alternatively, on all platforms, first install Picked Protein Group FDR as explained below. Then install `PyQt5` (`pip install PyQt5`) and run:

```shell
python gui.py
```

## Running Picked Protein Group FDR from the command line

First install Picked Protein Group FDR as explained below, then run:

```shell
python -m picked_group_fdr --mq_evidence </path/to/mq_evidence_txt> --fasta </path/to/fasta_file>
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
