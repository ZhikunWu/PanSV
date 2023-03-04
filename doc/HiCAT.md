## [HiCAT](https://github.com/xjtu-omics/HiCAT)


### install HiCAT
```
conda install -c conda-forge networkx
conda install -c conda-forge pandas
conda install -c conda-forge python-levenshtein
conda install -c bioconda stringdecomposer

```


### run HiCAT


```
(Assembly3) wuzhikun@mu03 17:13:47 ^_^ /home/wuzhikun/software/HiCAT 
$ python HiCAT.py -i ./testdata/cen21.fa -t ./testdata/AlphaSat.fa
Run stringdecomposer

python /home/wuzhikun/software/HiCAT/stringdecomposer/bin/stringdecomposer ./testdata/cen21.fa ./testdata/AlphaSat
.fa -o ./HiCAT_out2022-10-28 17:13:56,642 - StringDecomposer - INFO - cmd: ['/home/wuzhikun/software/HiCAT/stringdecomposer/bin/stri
ngdecomposer', './testdata/cen21.fa', './testdata/AlphaSat.fa', '-o', './HiCAT_out']2022-10-28 17:13:56,643 - StringDecomposer - INFO - Run /home/wuzhikun/software/HiCAT/stringdecomposer/stringdecom
poser/build/bin/dp with parameters ./testdata/cen21.fa ./testdata/AlphaSat.fa 1 5000 500 -1,-1,-1,1Scores: insertion=-1 deletion=-1 mismatch=-1 match=1
Prepared reads
100%: Aligned CP068257.1:11699867-12031015
2022-10-28 17:13:58,096 - StringDecomposer - INFO - Saved raw decomposition to ./HiCAT_out/final_decomposition_raw
.tsv2022-10-28 17:13:58,116 - StringDecomposer - INFO - Transforming raw alignments...
2022-10-28 17:13:58,344 - StringDecomposer - INFO - Transformation finished. Results can be found in ./HiCAT_out/f
inal_decomposition.tsv2022-10-28 17:13:58,344 - StringDecomposer - INFO - Thank you for using StringDecomposer!
Run HiCAT HOR

python /home/wuzhikun/software/HiCAT/HiCAT_HOR.py -d ./HiCAT_out/final_decomposition.tsv -b ./HiCAT_out/input_fast
a.1.fa -o ./HiCAT_out -s 0.94 -st 0.005 -m 40 -sp 5 -sn 10 -t 1start
build block sequence and read base sequence
calculate ed distance
ed distance thread: 1
1949
pre merge matrix
249
generation cover distribution for cluster
HOR thread: 1
get result
Time: 234

```


