

## [hal](https://github.com/ComparativeGenomicsToolkit/hal)




```
$ halExtract
acToo few (required positional) arguments

USAGE:
halExtract [Options] <inHalPath> <outHalPath>

ARGUMENTS:
inHalPath:    input hal file
outHalPath:   output hal file

OPTIONS:
--cacheBytes <value>:        obsolete name for --hdf5CacheBytes [default = 1048576]
--cacheMDC <value>:          obsolete name for --hdf5CacheMDC  [default = 113]
--cacheRDC <value>:          obsolete name for --hdf5CacheRDC [default = 521]
--cacheW0 <value>:           obsolete name for --hdf5CacheW0 [default = 0.75]
--chunk <value>:             obsolete name for --hdf5Chunk  [default = 1000]
--deflate <value>:           obsolete name for --hdf5Compression [default = 2]
--format <value>:            choose the back-end storage format. [default = hdf5]
--hdf5CacheBytes <value>:    maximum size in bytes of regular hdf5 cache [default = 
                             1048576]
--hdf5CacheMDC <value>:      number of metadata slots in hdf5 cache [default = 113]
--hdf5CacheRDC <value>:      number of regular slots in hdf5 cache.  should be a 
                             prime number ~= 10 * DefaultCacheRDCBytes / chunk 
                             [default = 521]
--hdf5CacheW0 <value>:       w0 parameter for hdf5 cache [default = 0.75]
--hdf5Chunk <value>:         hdf5 chunk size [default = 1000]
--hdf5Compression <value>:   hdf5 compression factor [0:none - 9:max] [default = 2]
--hdf5InMemory:              load all data in memory (and disable hdf5 cache) 
                             [default = 0]
--help:                      display this help page [default = 0]
--inMemory:                  obsolete name for --hdf5InMemory [default = 0]
--mmapFileSize <value>:      mmap HAL file initial size (in gigabytes) [default = 64]
--outputFormat <value>:      format for output hal file (same as input file by 
                             default) [default = ]
--root <value>:              root of subtree to extract [default = ""]

```


```
wuzhikun@mu02 18:07:34 O_O /home/wuzhikun/software/hal-release-V2.1 
$ maf2hal
Too few (required positional) arguments

maf2hal v2.2: import maf into hal database.

USAGE:
maf2hal [Options] <mafFile> <halFile>

ARGUMENTS:
mafFile:   output maf file
halFile:   input hal file

OPTIONS:
--append:                    append maf as subtree to existing alignment. reference 
                             must alaready be present in hal dabase as a leaf. 
                             [default = 0]
--cacheBytes <value>:        obsolete name for --hdf5CacheBytes [default = 1048576]
--cacheMDC <value>:          obsolete name for --hdf5CacheMDC  [default = 113]
--cacheRDC <value>:          obsolete name for --hdf5CacheRDC [default = 521]
--cacheW0 <value>:           obsolete name for --hdf5CacheW0 [default = 0.75]
--chunk <value>:             obsolete name for --hdf5Chunk  [default = 1000]
--deflate <value>:           obsolete name for --hdf5Compression [default = 2]
--format <value>:            choose the back-end storage format. [default = hdf5]
--hdf5CacheBytes <value>:    maximum size in bytes of regular hdf5 cache [default = 
                             1048576]
--hdf5CacheMDC <value>:      number of metadata slots in hdf5 cache [default = 113]
--hdf5CacheRDC <value>:      number of regular slots in hdf5 cache.  should be a 
                             prime number ~= 10 * DefaultCacheRDCBytes / chunk 
                             [default = 521]
--hdf5CacheW0 <value>:       w0 parameter for hdf5 cache [default = 0.75]
--hdf5Chunk <value>:         hdf5 chunk size [default = 1000]
--hdf5Compression <value>:   hdf5 compression factor [0:none - 9:max] [default = 2]
--hdf5InMemory:              load all data in memory (and disable hdf5 cache) 
                             [default = 0]
--help:                      display this help page [default = 0]
--inMemory:                  obsolete name for --hdf5InMemory [default = 0]
--mmapFileSize <value>:      mmap HAL file initial size (in gigabytes) [default = 64]
--refGenome <value>:         name of reference genome in MAF (first found if empty) 
                             [default = ]
--targetGenomes <value>:     comma-separated (no spaces) list of target genomes 
                             (others are excluded) (vist all if empty) [default = ]

```

```
(Assembly3) wuzhikun@mu02 18:13:45 ^_^ /home/wuzhikun/software/hal-release-V2.1 
$ halSynteny
Too few (required positional) arguments

halSynteny v2.2: Convert alignments into synteny blocks

USAGE:
halSynteny [Options] <alignment> <outPslPath>

ARGUMENTS:
alignment:    input file in HAL or PSL format (PSL must specify --alignmentIsPsl)
outPslPath:   output psl file ffor synteny blocks

OPTIONS:
--alignmentIsPsl:              alignment is in PSL format [default = 0]
--cacheBytes <value>:          obsolete name for --hdf5CacheBytes [default = 1048576]
--cacheMDC <value>:            obsolete name for --hdf5CacheMDC  [default = 113]
--cacheRDC <value>:            obsolete name for --hdf5CacheRDC [default = 521]
--cacheW0 <value>:             obsolete name for --hdf5CacheW0 [default = 0.75]
--format <value>:              choose the back-end storage format. [default = hdf5]
--hdf5CacheBytes <value>:      maximum size in bytes of regular hdf5 cache [default =
                                1048576]
--hdf5CacheMDC <value>:        number of metadata slots in hdf5 cache [default = 113]
--hdf5CacheRDC <value>:        number of regular slots in hdf5 cache.  should be a 
                               prime number ~= 10 * DefaultCacheRDCBytes / chunk 
                               [default = 521]
--hdf5CacheW0 <value>:         w0 parameter for hdf5 cache [default = 0.75]
--hdf5InMemory:                load all data in memory (and disable hdf5 cache) 
                               [default = 0]
--help:                        display this help page [default = 0]
--inMemory:                    obsolete name for --hdf5InMemory [default = 0]
--maxAnchorDistance <value>:   upper bound on distance for syntenic psl blocks 
                               [default = 5000]
--minBlockSize <value>:        lower bound on synteny block length [default = 5000]
--queryChromosome <value>:     chromosome to infer synteny (default is whole genome) 
                               [default = ""]
--queryGenome <value>:         source genome [default = ""]
--targetGenome <value>:        reference genome name [default = ""]

```
