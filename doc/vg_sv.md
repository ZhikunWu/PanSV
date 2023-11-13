
## [vg](https://github.com/vgteam/vg)

[Learning vg on toy examples](https://gtpb.github.io/CPANG19/pages/toy_examples)

### parameters

```
$ /home/wuzhikun/software/vg --version
error:[vg] command --version not found
vg: variation graph tool, version v1.50.1 "Monopoli"

usage: /home/wuzhikun/software/vg <command> [options]

main mapping and calling pipeline:
  -- autoindex     mapping tool-oriented index construction from interchange formats
  -- construct     graph construction
  -- rna           construct splicing graphs and pantranscriptomes
  -- index         index graphs or alignments for random access or mapping
  -- map           MEM-based read alignment
  -- giraffe       fast haplotype-aware short read alignment
  -- mpmap         splice-aware multipath alignment of short reads
  -- augment       augment a graph from an alignment
  -- pack          convert alignments to a compact coverage index
  -- call          call or genotype VCF variants
  -- help          show all subcommands

For more commands, type `vg help`.
For technical support, please visit: https://www.biostars.org/tag/vg/

```

### manual

```
# specify data and prefix
myfa="test.fa"
myvcf="test.vcf.gz"
myprefix="test"

# first autoindex the graph for alignment
vg autoindex --workflow giraffe --prefix ${myprefix} --ref-fasta ${myfa} --vcf ${myvcf} --threads 60 --tmp-dir ../tmp

# giraffe mapping short reads to graph
mkdir mygam
cat sample.lst | parallel -j 5 -k "vg giraffe -Z ${myprefix}.giraffe.gbz -m ${myprefix}.min -d ${myprefix}.dist -f ../GBS/{}.fastq.trim.fq -t 12 -N {} > mygam/{}.gam"

# count read support
mkdir mypack
cat sample.lst | parallel -j 5 -k "vg pack -x ${myprefix}.giraffe.gbz -t 12 -g mygam/{}.gam -o mypack/{}.pack -Q 5 -s 5"

# calculate snarl
vg snarls ${myprefix}.giraffe.gbz > ${myprefix}.snarls

# genotype the variants
mkdir myvcf
cat sample.lst | parallel -j 10 -k "vg call ${myprefix}.giraffe.gbz -z -a -k mypack/{}.pack -s {} -t 6 -r ${myprefix}.snarls > myvcf/{}.vcf"

```




### run
```
$ vg autoindex --workflow giraffe -r /home/wuzhikun/Project/Vigna/Final/Final/Vigna_unguiculata_assembly.fasta -v /home/wuzhikun/Project/PanVigna/SVCall/PanGraph/Samples_SV_merge_vgformat.vcf.gz -p /home/wuzhikun/Project/PanVigna/SVCall/PanGraph/Samples_SV_merge_giraffe
[vg autoindex] Executing command: vg autoindex --workflow giraffe -r /home/wuzhikun/Project/Vigna/Final/Final/Vigna_unguiculata_assembly.fasta -v /home/wuzhikun/Project/PanVigna/SVCall/PanGraph/Samples_SV_merge_vgformat.vcf.gz -p /home/wuzhikun/Project/PanVigna/SVCall/PanGraph/Samples_SV_merge_giraffe
[IndexRegistry]: Checking for phasing in VCF(s).
[IndexRegistry]: Chunking inputs for parallelism.
[IndexRegistry]: Chunking FASTA(s).
[IndexRegistry]: Chunking VCF(s).
[IndexRegistry]: Constructing VG graph from FASTA and VCF input.
Warning: could not set SEQ [canonicalize]. Vu05 1874540 Sniffles2.INS.204M4 C   <INS>   60  PASS    AC=8;COVERAGE=181,146,143,115,136;END=1874540;SPAN=2092;STDEV_LEN=172.461;STDEV_POS=0;STRAND=+-;SUPPORT=40;SUPP_VEC=10111100000000000100000001111000001;SVLEN=2092;SVTYPE=INS;PRECISE
Warning: could not set SEQ [canonicalize]. Vu08 7729171 Sniffles2.INS.849M7 G   <INS>   59  PASS    AC=1;COVERAGE=398,49,49,49,236;END=7729171;SPAN=4556;STDEV_LEN=57.983;STDEV_POS=0.707;STRAND=+-;SUPPORT=8;SUPP_VEC=00000000000000000000000000001001000;SVLEN=4556;SVTYPE=INS;IMPRECISE
Warning: could not set SEQ [canonicalize]. Vu01 9643044 Sniffles2.INS.60BM0 G   <INS>   60  PASS    AC=5;COVERAGE=118,110,110,110,180;END=9643044;SPAN=1565;STDEV_LEN=15.175;STDEV_POS=2.229;STRAND=+-;SUPPORT=37;SUPP_VEC=00000000010000000000001100001101000;SVLEN=1565;SVTYPE=INS;PRECISE
Warning: could not set SEQ [canonicalize]. Vu10 19910054    Sniffles2.INS.237FM9    A   <INS>   43  PASS    AC=1;COVERAGE=70,70,70,70,73;END=19910054;SPAN=2060;STDEV_LEN=397.394;STDEV_POS=92.631;STRAND=+-;SUPPORT=12;SUPP_VEC=00000000000000100000000000100000000;SVLEN=2060;SVTYPE=INS;IMPRECISE
Warning: could not set SEQ [canonicalize]. Vu01 25660163    Sniffles2.INS.16C3M0    A   <INS>   57  PASS    AC=1;COVERAGE=67,34,34,34,73;END=25660163;SPAN=1486;STDEV_LEN=16.058;STDEV_POS=7.026;STRAND=+-;SUPPORT=3;SUPP_VEC=00000000100000100000010000010100100;SVLEN=1486;SVTYPE=INS;PRECISE
Warning: could not set SEQ [canonicalize]. Vu03 29859260    Sniffles2.INS.1A21M2    C   <INS>   60  PASS    AC=1;COVERAGE=82,32,32,32,71;END=29859260;SPAN=11828;STDEV_LEN=0;STDEV_POS=0;STRAND=+-;SUPPORT=8;SUPP_VEC=00000000000000000000000000000010000;SVLEN=11828;SVTYPE=INS;IMPRECISE
Warning: could not set SEQ [canonicalize]. Vu11 30186523    Sniffles2.INS.2F0EMA    G   <INS>   53  PASS    AC=4;COVERAGE=68,63,60,60,51;END=30186523;SPAN=1406;STDEV_LEN=18.607;STDEV_POS=1.251;STRAND=+-;SUPPORT=4;SUPP_VEC=10101111111111111111011111011111111;SVLEN=1406;SVTYPE=INS;PRECISE
[IndexRegistry]: Constructing XG graph from VG graph.
[IndexRegistry]: Constructing a greedy path cover GBWT
[IndexRegistry] forked child 52763
[IndexRegistry]: Constructing GBZ.
[IndexRegistry]: Constructing distance index for Giraffe.
[IndexRegistry]: Constructing minimizer index.

```


files:
```
-rw-rw-r-- 1 295M Aug 15 13:45 Samples_SV_merge_giraffe.dist
-rw-rw-r-- 1 525M Aug 15 13:35 Samples_SV_merge_giraffe.giraffe.gbz
-rw-rw-r-- 1 5.0G Aug 15 13:45 Samples_SV_merge_giraffe.min
```


```
$ vg convert -x Samples_SV_merge_giraffe.giraffe.gbz  > Samples_SV_merge_giraffe.giraffe.xg
```



```
$ vg giraffe --gbz-name /home/wuzhikun/Project/PanVigna/SVCall/PanGraph/Samples_SV_merge_giraffe.giraffe.gbz --dist-name /home/wuzhikun/Project/PanVigna/SVCall/PanGraph/Samples_SV_merge_giraffe.dist --minimizer-name /home/wuzhikun/Project/PanVigna/SVCall/PanGraph/Samples_SV_merge_giraffe.min --xg-name /home/wuzhikun/Project/PanVigna/SVCall/PanGraph/Samples_SV_merge_giraffe.giraffe.xg  --fastq-in /home/wuzhikun/Project/PanVigna/raw/SRR22136022_NGS.R1.fastq.gz --fastq-in /home/wuzhikun/Project/PanVigna/raw/SRR22136022_NGS.R2.fastq.gz   --threads 10 --sample SRR22136022  -o BAM  >  SRR22136022.bam


```


```
vg pack -t $task.cpus -x $gbz -g $gam -Q5 -o ${prefix}.aln.pack
vg call -t $task.cpus $gbz -C 100 -k ${prefix}.aln.pack --min-support $params.min_support -a -r $snarls -z -s $sample_name > ${prefix}.vcf
bgzip ${prefix}.vcf
tabix ${prefix}.vcf.gz
```

