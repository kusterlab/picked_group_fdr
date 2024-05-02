Manuscript figures
==================

To reproduce the figures from the manuscript:

#. make sure that you have ``make`` installed.
#. download the input files from zenodo: https://zenodo.org/record/7157677
#. specify the location of the input files and run the make command with the figure you want to reproduce, e.g.:

    .. code:: shell

        export DATA_DIR=</path/to/zenodo/files>
        make Figure3a

Check the file ``MakefileFigures`` to see which figures can be plotted.