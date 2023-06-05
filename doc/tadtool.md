### tadtool parameters


```
$ tadtool
usage: tadtool <command> [options]

Commands:
    plot                Main interactive TADtool plotting window
    tads                Call TADs with pre-defined parameters
    subset              Reduce a matrix to a smaller region.

Run tadtool <command> -h for help on a specific command.
tadtool: error: the following arguments are required: command

```



```
$ tadtool tads --help
usage: tadtool tads [-h] [-a ALGORITHM] [-n NORMALISATION_WINDOW] [-v] matrix regions window_size cutoff [output]

Call TADs with pre-defined parameters

positional arguments:
  matrix                Square Hi-C Matrix as tab-delimited or .npy file (created with numpy.save) or sparse matrix format
                        (each line: <row region index> <column region index> <matrix value>)
  regions               BED file (no header) with regions corresponding to the number of rows in the provided matrix.
                        Fourth column, if present, denotes name field, which is used as an identifier in sparse matrix
                        notation.
  window_size           Window size in base pairs
  cutoff                Cutoff for TAD-calling algorithm at given window size.
  output                Optional output file to save TADs.

optional arguments:
  -h, --help            show this help message and exit
  -a ALGORITHM, --algorithm ALGORITHM
                        TAD-calling algorithm. Options: insulation, ninsulation, directionality. Default: insulation.
  -n NORMALISATION_WINDOW, --normalisation-window NORMALISATION_WINDOW
                        Normalisation window in number of regions. Only affects ninsulation algorithm. If not specified,
                        window will be the whole chromosome.
  -v, --write-values    Write index values to file instead of TADs.

```



tadtool/examples/chr12_20-35Mb.matrix_sparse.txt

```
6   497 2.740072e-01
8   48  2.820830e-02
8   351 7.410965e-02
8   411 8.491875e-02
8   488 7.516951e-02
8   524 8.719267e-02
9   10  2.692810e-01
9   28  1.870377e-01
9   181 1.663241e-02
```


