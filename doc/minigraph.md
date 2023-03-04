
## [minigraph](https://github.com/lh3/minigraph)


### install minigraph

```
git clone https://github.com/lh3/minigraph
cd minigraph && make
```


### parameters
```
$ minigraph
Usage: minigraph [options] <target.gfa> <query.fa> [...]
Options:
  Indexing:
    -k INT       k-mer size (no larger than 28) [17]
    -w INT       minizer window size [11]
  Mapping:
    -c           perform base alignment; RECOMMENDED
    -f FLOAT     ignore top FLOAT fraction of repetitive minimizers [0.0002]
    -U INT[,INT] choose the minimizer occurrence threshold within this interval [50,250]
    -j FLOAT     expected sequence divergence [0.1]
    -g NUM       stop chain enlongation if there are no minimizers in INT-bp [5000]
    -F NUM       max fragment length (effective with -xsr or in the fragment mode) [0]
    -r NUM[,NUM] bandwidth for the two rounds of chaining [500,20000]
    -n INT[,INT] minimal number of minimizers on a graph/linear chain [5,5]
    -m INT[,INT] minimal graph/linear chaining score [50,40]
    -p FLOAT     min secondary-to-primary score ratio [0.8]
    -N INT       retain at most INT secondary mappings [5]
    -D           skip self diagonal matches
  Graph generation:
    --ggen       perform incremental graph generation
    -q INT       min mapping quality [5]
    -l NUM       min alignment length [100000]
    -d NUM       min alignment length for depth calculation [20000]
    -L INT       min variant length [50]
    --call       call the graph path in each bubble and output BED
  Input/output:
    -t INT       number of threads [4]
    -o FILE      output mappings to FILE [stdout]
    -K NUM       minibatch size for mapping [500M]
    -S           output linear chains in * sName sLen nMz div sStart sEnd qStart qEnd
    --vc         output in the vertex coordinate
  Preset:
    -x STR       preset []
                 - lr: noisy long read mapping (the default)
                 - asm: asm-to-ref mapping
                 - sr: short reads
                 - ggs: incremental graph generation

```


