
## [hifiasm_meta](https://github.com/xfengnefx/hifiasm-meta)


### parameter
```
$ /home/wuzhikun/software/hifiasm_meta
Usage: hifiasm_meta 0.3-r073 (hifiasm code base 0.13-r308)
Options:
  Input/Output:
    -o STR      prefix of output files [hifiasm_meta.asm]
    -B STR      prefix of bin files, if it s different from -o [hifiasm_meta.asm]
    -i          ignore saved read correction and overlaps
    -t INT      number of threads [1]
    -z INT      length of adapters that should be removed [0]
    --version   show version number
  Read selection:
    --force-rs
                enable and force read selection.
    --lowq-10
                lower 10% runtime kmer frequency threshold. [50]
    --lowq-5
                lower 5% runtime kmer frequency threshold. [50]
    --lowq-3
                lower 3% runtime kmer frequency threshold. [10]
  Overlap/Error correction:
    -k INT      k-mer length (must be <64) [51]
    -w INT      minimizer window size [51]
    -f INT      number of bits for bloom filter; 0 to disable [37]
    -D FLOAT    drop k-mers occurring >FLOAT*coverage times [5.0]
    -N INT      consider up to max(-D*coverage,-N) overlaps for each oriented read [100]
    -r INT      round of correction [3]
  Assembly:
    -a INT      round of assembly cleaning [4]
    -n INT      remove tip unitigs composed of <=INT reads [3]
    -x FLOAT    max overlap drop ratio [0.8]
    -y FLOAT    min overlap drop ratio [0.2]
    --noch      disable contained reads sparing heuristics.
  Auxiliary:
    -e          ban assembly, i.e. terminate before generating string graph
    --write-paf dump overlaps (paf).
    --dump-all-ovlp
                dump all overlaps ever calculated in the final overlapping (paf).
    --write-ec
                dump error corrected reads (fasta).
Example: ./hifiasm_meta -o asm -t 32 asm.fq.gz 2>log
See `man ./hifiasm_meta.1` for detailed descriptions command-line options.
```


### run



```
$ /home/wuzhikun/software/hifiasm_meta -o organelle_assembly -t 20  /home/wuzhikun/Project/BAssembly/pipeline/organelle/MH388764.1_hifi_align.id.10k.fasta
```


```
********** checkpoint: p_ctg **********

[M::hamt_clean_graph] (peak RSS so far: 19.7 GB)
[M::hamt_clean_graph] removed single-read unconnected contigs, 1 on primary, 0 on alt.
[M::hamt_output_unitig_graph_advance] Writing GFA... 
[M::hamt_output_unitig_graph_advance] Writing GFA... 


********** checkpoint: post-assembly **********

[M::hamt_clean_graph] (peak RSS so far: 19.7 GB)
[M::hamt_ug_opportunistic_elementary_circuits] collected 0 circuits, used 0.00s
[M::hamt_ug_opportunistic_elementary_circuits] wrote all rescued circles, used 0.00s
[T::hamt_ug_opportunistic_elementary_circuits_helper_deduplicate_minhash] got the sequences, used 0.0s
[T::hamt_minhash_mashdist] sketched - 0.0s.
[T::hamt_minhash_mashdist] compared - 0.0s.
[T::hamt_ug_opportunistic_elementary_circuits_helper_deduplicate_minhash] collected mash distances for 0 seqs, used 0.0s
[M::hamt_ug_opportunistic_elementary_circuits_helper_deduplicate_minhash] had 0 paths, 0 remained (0 dropped by length diff, 0 by length abs),used 0.0s after sketching.
[M::hamt_ug_opportunistic_elementary_circuits] deduplicated rescued circles, used 0.00s
[M::hamt_ug_opportunistic_elementary_circuits] wrote deduplicated rescued circles, used 0.00s
[M::hamt_simple_binning] Will try to bin on 3 contigs (skipped 0 because blacklist).
Using random seed: 42
Perplexity too large for the number of data points!

```

```
gfatools gfa2fa organelle_assembly.p_ctg.gfa > organelle_assembly.p_ctg.fasta
```
