
## [pangene](https://github.com/lh3/pangene)


### install
```
git clone https://github.com/lh3/pangene
cd pangene && make
```



### parameters

```
$ /home/wuzhikun/software/pangene/pangene 
Usage: pangene [options] <in.paf> [...]
Options:
  Input preprocessing:
    -d CHAR       gene-protein delimiter [:]
    -X STR/@FILE  exclude genes in STR list or in @FILE []
    -I STR/@FILE  attempt to include genes in the output graph []
    -e FLOAT      drop an alignment if its identity <FLOAT [0.5]
    -l FLOAT      drop an alignment if <FLOAT fraction of the protein aligned [0.5]
    -m FLOAT      score adjustment coefficient [2]
  Graph construction:
    -f FLOAT      min overlap fraction [0.5]
    -p FLOAT      gene considered if dominant in FLOAT fraction of genes [0.05]
    -c INT        drop a gene if average occurrence is >INT [10]
    -g INT        drop a gene if its in- or out-degree >INT [15]
    -r INT        drop a gene if it connects >INT distant loci [3]
    -b FLOAT      demote a branching arc if weaker than the best by FLOAT [0.02]
    -B FLOAT      cut a branching arc if weaker by FLOAT [0.5]
    -y FLOAT      cut a distant branching arc if weaker by FLOAT [0.05]
    -T INT        apply branch cutting for INT times [15]
    -F            dont consider genes on different contigs as distant
    -a INT        prune an arc if it is supported by <INT genomes [1]
  Output:
    -w            Suppress walk lines (W-lines)
    --bed[=STR]   output 12-column BED where STR is walk, raw or flag [walk]
    --version     print version number

```



