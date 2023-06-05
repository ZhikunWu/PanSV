

## [SVision](https://github.com/xjtu-omics/SVision)

### install svison

```
## Get latest source code
git clone https://github.com/xjtu-omics/SVision.git
cd SVision

## Create conda environment and install SVision 
conda env create -f environment.yml
python setup.py install
```

```
SVision -o ./home/user/svision_out -b ./supports/HG00733.svision.demo.bam -m /home/user/svision_model/svision-cnn-model.ckpt -g /path/to/reference.fa -n HG00733 -s 5 --graph --qname
```




```

(svisionenv) wuzhikun@fat02 09:26:42 ^_^ /home/wuzhikun/Project/Population 
$ SVision --help
usage: SVision [-h] -o OUT_PATH -b BAM_PATH -m MODEL_PATH -g GENOME -n SAMPLE
               [-t THREAD_NUM] [-s MIN_SUPPORT] [-c CHROM] [--hash]
               [--cluster] [--qname] [--graph] [--contig]
               [--min_mapq MIN_MAPQ] [--min_sv_size MIN_SV_SIZE]
               [--max_sv_size MAX_SV_SIZE] [--window_size WINDOW_SIZE]
               [--patition_max_distance PATITION_MAX_DISTANCE]
               [--cluster_max_distance CLUSTER_MAX_DISTANCE]
               [--batch_size BATCH_SIZE] [--min_gt_depth MIN_GT_DEPTH]
               [--homo_thresh HOMO_THRESH] [--hete_thresh HETE_THRESH]
               [--k_size K_SIZE] [--min_accept MIN_ACCEPT]
               [--max_hash_len MAX_HASH_LEN]

SVision 1.3.6 
 
Short Usage: SVision [parameters] -o <output path> -b <input bam path> -g <reference> -m <model path>

optional arguments:
  -h, --help            show this help message and exit

Input/Output parameters:
  -o OUT_PATH           Absolute path to output
  -b BAM_PATH           Absolute path to bam file
  -m MODEL_PATH         Absolute path to CNN predict model
  -g GENOME             Absolute path to your reference genome (.fai required
                        in the directory)
  -n SAMPLE             Name of the BAM sample name

Optional parameters:
  -t THREAD_NUM         Thread numbers (default: 1)
  -s MIN_SUPPORT        Minimum support read number required for SV calling
                        (default: 5)
  -c CHROM              Specific region (chr1:xxx-xxx) or chromosome (chr1) to
                        detect
  --hash                Activate local realignment for unmapped sequences
                        (default: False)
  --cluster             Cluster calls that might occur together (default:
                        False)
  --qname               Report support names for each events (default: False)
  --graph               Report graph for events (default: False)
  --contig              Activate contig mode (default: False)

Collect parameters:
  --min_mapq MIN_MAPQ   Minimum mapping quality of reads to consider (default:
                        10)
  --min_sv_size MIN_SV_SIZE
                        Minimum SV size to detect (default: 50)
  --max_sv_size MAX_SV_SIZE
                        Maximum SV size to detect (default: 1000000)
  --window_size WINDOW_SIZE
                        The sliding window size in segment collection
                        (default: 10000000)

Cluster parameters:
  --patition_max_distance PATITION_MAX_DISTANCE
                        Maximum distance to partition signatures (default:
                        5000)
  --cluster_max_distance CLUSTER_MAX_DISTANCE
                        Clustering maximum distance for a partition (default:
                        0.3)

Predict parameters:
  --batch_size BATCH_SIZE
                        Batch size for the CNN prediction model (default: 128)

Genotype parameters:
  --min_gt_depth MIN_GT_DEPTH
                        Minimum reads required for genotyping (default: 4)
  --homo_thresh HOMO_THRESH
                        Minimum variant allele frequency to be called as
                        homozygous (default: 0.8)
  --hete_thresh HETE_THRESH
                        Minimum variant allele frequency to be called as
                        heterozygous (default: 0.2)

Hash table parameters:
  --k_size K_SIZE       Size of kmer (default: 10)
  --min_accept MIN_ACCEPT
                        Minimum match length for realignment (default: 50)
  --max_hash_len MAX_HASH_LEN
                        Maximum length of unmapped sequence length for
                        realignment (default: 1000)

```


```
SVision  -o /home/wuzhikun/Project/HuaXi2021/mapping/minimap2/M692-0_SV -b /home/wuzhikun/Project/HuaXi2021/mapping/minimap2/M692-0.bam -m /home/zhengjingjing/software/SVision/svision_model/svision-cnn-model.ckpt -g /home/wuzhikun/database/genome/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa -n M692-0 -s 3 --graph --qname -t 24 

```




[output vcf format](https://github.com/xjtu-omics/SVision/wiki/Output-format)


```
SVision adopts the standard VCF format with extra info columns. Some important info columns are listed as below:

The SV ID column is given in the format of a_b, where b indicates site a contains other type of SVs.

Filters used in the output.

Covered: The entire SV is spanned by long-reads, producing the most confident calls.

Uncovered: SV is partially spanned by long-reads, i.e. reads spanning one of the breakpoints.

Clustered: SV is partially spanned by long-reads, but can be spanned through reads clusters.

We add extra attributes in the INFO column of VCF format for SVision detected complex structural variants.

BRPKS: The CNN recognized internal structure of CSVs through tMOR.

GraphID: The graph index used to indicate the graph structure, which requires --graph and is obtained by calculating isomorphic graphs. The ID for simple SVs is -1.

GFA_FILE_PREFIX: File name of CSV corresponding GFA file.

GFA_S: Nodes contained in a CSV graph represented based on GFA format.

GFA_L: Links contained in a CSV graph represented based on GFA format

Example of a SVision CSV call from the demo data
```

