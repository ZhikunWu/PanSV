
## cooler

```
$ cooler cload pairs -h
Usage: cooler cload pairs [OPTIONS] BINS PAIRS_PATH COOL_PATH

  Bin any text file or stream of pairs.

  Pairs data need not be sorted. Accepts compressed files. To pipe input from
  stdin, set PAIRS_PATH to '-'.

  BINS : One of the following

      <TEXT:INTEGER> : 1. Path to a chromsizes file, 2. Bin size in bp

      <TEXT> : Path to BED file defining the genomic bin segmentation.

  PAIRS_PATH : Path to contacts (i.e. read pairs) file.

  COOL_PATH : Output COOL file path or URI.

Options:
  --metadata TEXT                 Path to JSON file containing user metadata.
  --assembly TEXT                 Name of genome assembly (e.g. hg19, mm10)
  -c1, --chrom1 INTEGER           chrom1 field number (one-based)  [required]
  -p1, --pos1 INTEGER             pos1 field number (one-based)  [required]
  -c2, --chrom2 INTEGER           chrom2 field number (one-based)  [required]
  -p2, --pos2 INTEGER             pos2 field number (one-based)  [required]
  --chunksize INTEGER             Number of input lines to load at a time
  -0, --zero-based                Positions are zero-based
  --comment-char TEXT             Comment character that indicates lines to
                                  ignore.  [default: #]
  -N, --no-symmetric-upper        Create a complete square matrix without
                                  implicit symmetry. This allows for distinct
                                  upper- and lower-triangle values
  --input-copy-status [unique|duplex]
                                  Copy status of input data when using
                                  symmetric-upper storage. | `unique`:
                                  Incoming data comes from a unique half of a
                                  symmetric map, regardless of how the
                                  coordinates of a pair are ordered. `duplex`:
                                  Incoming data contains upper- and lower-
                                  triangle duplicates. All input records that
                                  map to the lower triangle will be discarded!
                                  | If you wish to treat lower- and upper-
                                  triangle input data as distinct, use the
                                  ``--no-symmetric-upper`` option.   [default:
                                  unique]
  --field TEXT                    Specify quantitative input fields to
                                  aggregate into value columns using the
                                  syntax ``--field <field-name>=<field-
                                  number>``. Optionally, append ``:`` followed
                                  by ``dtype=<dtype>`` to specify the data
                                  type (e.g. float), and/or ``agg=<agg>`` to
                                  specify an aggregation function different
                                  from sum (e.g. mean). Field numbers are
                                  1-based. Passing 'count' as the target name
                                  will override the default behavior of
                                  storing pair counts. Repeat the ``--field``
                                  option for each additional field.
  --temp-dir DIRECTORY            Create temporary files in a specified
                                  directory. Pass ``-`` to use the platform
                                  default temp dir.
  --no-delete-temp                Do not delete temporary files when finished.
  --max-merge INTEGER             Maximum number of chunks to merge before
                                  invoking recursive merging  [default: 200]
  --storage-options TEXT          Options to modify the data filter pipeline.
                                  Provide as a comma-separated list of key-
                                  value pairs of the form 'k1=v1,k2=v2,...'.
                                  See http://docs.h5py.org/en/stable/high/data
                                  set.html#filter-pipeline for more details.
  -h, --help                      Show this message and exit.

```




```
$ cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 /home/wuzhikun/Project/Centromere/Cooler/Arabidopsis_thaliana/Chromap/Arabidopsis_thaliana.ChromSize.txt:1000 /home/wuzhikun/Project/Centromere/Cooler/Arabidopsis_thaliana/Pairs/Arabidopsis_thaliana.sorted.pairs_deduplicate.gz  /home/wuzhikun/Project/Centromere/Cooler/Arabidopsis_thaliana/Pairs/Arabidopsis_thaliana.sorted.pairs_deduplicate_cooler_1000.cool

WARNING:py.warnings:/home/wuzhikun/anaconda3/envs/HiC-Pro/lib/python3.8/site-packages/cooler/util.py:733: FutureWarning: is_categorical is deprecated and will be removed in a future version. Use is_categorical_dtype instead.
  is_cat = pd.api.types.is_categorical(bins["chrom"])

INFO:cooler.create:Writing chunk 0: /home/wuzhikun/Project/Centromere/Cooler/Arabidopsis_thaliana/Pairs/tmpmmq86vr5.multi.cool::0
INFO:cooler.create:Creating cooler at "/home/wuzhikun/Project/Centromere/Cooler/Arabidopsis_thaliana/Pairs/tmpmmq86vr5.multi.cool::/0"
INFO:cooler.create:Writing chroms
INFO:cooler.create:Writing bins
INFO:cooler.create:Writing pixels
INFO:cooler.create:Writing indexes
INFO:cooler.create:Writing info
INFO:cooler.create:Merging into /home/wuzhikun/Project/Centromere/Cooler/Arabidopsis_thaliana/Pairs/Arabidopsis_thaliana.sorted.pairs_deduplicate_cooler_1000.cool
INFO:cooler.create:Creating cooler at "/home/wuzhikun/Project/Centromere/Cooler/Arabidopsis_thaliana/Pairs/Arabidopsis_thaliana.sorted.pairs_deduplicate_cooler_1000.cool::/"
INFO:cooler.create:Writing chroms
INFO:cooler.create:Writing bins
INFO:cooler.create:Writing pixels
INFO:cooler.reduce:nnzs: [130188]
INFO:cooler.reduce:current: [130188]
INFO:cooler.create:Writing indexes
INFO:cooler.create:Writing info

```


```
cooler zoomify  --nproc 5  --out /home/wuzhikun/Project/Centromere/Cooler/Arabidopsis_thaliana/Pairs/Arabidopsis_thaliana.sorted.pairs_deduplicate_cooler_1000.mcool  --resolutions 1000,2000  --balance     /home/wuzhikun/Project/Centromere/Cooler/Arabidopsis_thaliana/Pairs/Arabidopsis_thaliana.sorted.pairs_deduplicate_cooler_1000.cool

INFO:cooler.cli.zoomify:Recursively aggregating "/home/wuzhikun/Project/Centromere/Cooler/Arabidopsis_thaliana/Pairs/Arabidopsis_thaliana.sorted.pairs_deduplicate_cooler_1000.cool"
INFO:cooler.cli.zoomify:Writing to "/home/wuzhikun/Project/Centromere/Cooler/Arabidopsis_thaliana/Pairs/Arabidopsis_thaliana.sorted.pairs_deduplicate_cooler_1000.mcool"
INFO:cooler.reduce:Copying base matrices and producing 2 new zoom levels.
INFO:cooler.reduce:Bin size: 1000
INFO:cooler.reduce:Aggregating from 1000 to 2000.
WARNING:py.warnings:/home/wuzhikun/anaconda3/envs/HiC-Pro/lib/python3.8/site-packages/cooler/util.py:733: FutureWarning: is_categorical is deprecated and will be removed in a future version. Use is_categorical_dtype instead.
  is_cat = pd.api.types.is_categorical(bins["chrom"])

INFO:cooler.create:Creating cooler at "/home/wuzhikun/Project/Centromere/Cooler/Arabidopsis_thaliana/Pairs/Arabidopsis_thaliana.sorted.pairs_deduplicate_cooler_1000.mcool::/resolutions/2000"
INFO:cooler.create:Writing chroms
INFO:cooler.create:Writing bins
INFO:cooler.create:Writing pixels
INFO:cooler.reduce:0 130188
INFO:cooler.create:Writing indexes
INFO:cooler.create:Writing info
INFO:cooler.cli.zoomify:Balancing zoom level with bin size 1000
INFO:cooler.cli.balance:Balancing "/home/wuzhikun/Project/Centromere/Cooler/Arabidopsis_thaliana/Pairs/Arabidopsis_thaliana.sorted.pairs_deduplicate_cooler_1000.mcool::resolutions/1000"
WARNING:py.warnings:/home/wuzhikun/anaconda3/envs/HiC-Pro/lib/python3.8/site-packages/numpy/core/fromnumeric.py:3440: RuntimeWarning: Mean of empty slice.
  return _methods._mean(a, axis=axis, dtype=dtype,

WARNING:py.warnings:/home/wuzhikun/anaconda3/envs/HiC-Pro/lib/python3.8/site-packages/numpy/core/_methods.py:189: RuntimeWarning: invalid value encountered in double_scalars
  ret = ret.dtype.type(ret / rcount)

WARNING:py.warnings:/home/wuzhikun/anaconda3/envs/HiC-Pro/lib/python3.8/site-packages/cooler/balance.py:114: RuntimeWarning: Mean of empty slice.
  scale = nzmarg.mean()

INFO:cooler.cli.zoomify:Balancing zoom level with bin size 2000
INFO:cooler.cli.balance:Balancing "/home/wuzhikun/Project/Centromere/Cooler/Arabidopsis_thaliana/Pairs/Arabidopsis_thaliana.sorted.pairs_deduplicate_cooler_1000.mcool::resolutions/2000"
WARNING:py.warnings:/home/wuzhikun/anaconda3/envs/HiC-Pro/lib/python3.8/site-packages/numpy/core/fromnumeric.py:3440: RuntimeWarning: Mean of empty slice.
  return _methods._mean(a, axis=axis, dtype=dtype,

WARNING:py.warnings:/home/wuzhikun/anaconda3/envs/HiC-Pro/lib/python3.8/site-packages/numpy/core/_methods.py:189: RuntimeWarning: invalid value encountered in double_scalars
  ret = ret.dtype.type(ret / rcount)

WARNING:py.warnings:/home/wuzhikun/anaconda3/envs/HiC-Pro/lib/python3.8/site-packages/cooler/balance.py:114: RuntimeWarning: Mean of empty slice.
  scale = nzmarg.mean()

```

