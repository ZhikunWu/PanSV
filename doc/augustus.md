

### augustus

### [Predicting Genes with AUGUSTUS](https://bioinf.uni-greifswald.de/augustus/binaries/tutorial/prediction.html)


```
将 hints.gff 文件内容和 genome.fasta 内容进行分割。不对完整的序列进行切断。此程序将基因组序列按长度进行排序后，将序列写入到一个个分割的fasta文件中。每当写入到一个fasta文件的序列长度大于设置的值时，则将下一条序列下如下一个fasta文件。同时，也将相应的 hints 信息写入对应的 hint 文件。
$ perl split_hints_and_scaffolds_for_augustus.pl --minsize 1000000 --output split genome.fasta hints.gff

并行计算
$ for x in `ls split/*.fa | perl -pe 's/.*\///; s/.fa//' | sort -k 1.14n`
do
    echo "augustus --species=my_species --extrinsicCfgFile=extrinsic.cfg --alternatives-from-evidence=true --hintsfile=split/$x.hints --allow_hinted_splicesites=atac --alternatives-from-evidence=true --gff3=on --UTR=on split/$x.fa > split/$x.out" >> command_augustus.list
done
$ ParaFly -c command_augustus.list -CPU 12
 
合并结果
$ for x in `ls split/*.out | sort -k 1.20n`
do
    cat $x >> aug.out
done

$ join_aug_pred.pl aug.out > aug.gff3
```




### run using exist model


(Braker) wuzhikun@mu03 17:48:58 ^_^ /home/wuzhikun/Project/BAssembly/pipeline/august 
```
$ augustus --species=arabidopsis --gff3=on /home/wuzhikun/Project/BAssembly/Assembly/hifiOnt100kHic/hifiOnt100kHic.hic.p_ctg_200k.fasta.split/hifiOnt100kHic.hic.p_ctg_200k.id_ptg000010l.fasta  > hifiOnt100kHic.hic.p_ctg_200k.id_ptg000010l_august.gff3
```


out file:
```

$ grep -v '^#' hifiOnt100kHic.hic.p_ctg_200k.id_ptg000010l_august.gff3 | head 
ptg000010l  AUGUSTUS    gene    3489    5301    0.01    -   .   ID=g1
ptg000010l  AUGUSTUS    transcript  3489    5301    0.01    -   .   ID=g1.t1;Parent=g1
ptg000010l  AUGUSTUS    transcription_end_site  3489    3489    .   -   .   Parent=g1.t1
ptg000010l  AUGUSTUS    exon    3489    3789    .   -   .   Parent=g1.t1
ptg000010l  AUGUSTUS    stop_codon  3741    3743    .   -   0   Parent=g1.t1
ptg000010l  AUGUSTUS    intron  3790    4115    0.8 -   .   Parent=g1.t1
ptg000010l  AUGUSTUS    intron  4223    4369    1   -   .   Parent=g1.t1
ptg000010l  AUGUSTUS    intron  4500    4915    0.42    -   .   Parent=g1.t1
ptg000010l  AUGUSTUS    intron  5092    5180    0.53    -   .   Parent=g1.t1
ptg000010l  AUGUSTUS    CDS 3741    3789    0.77    -   1   ID=g1.t1.cds;Parent=g1.t1
```

### If you also want the protein sequences you can retrieve them with

```
getAnnoFasta.pl augustus.abinitio.gff
```

