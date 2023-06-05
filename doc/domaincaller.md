

## [domaincaller](https://github.com/XiaoTaoWang/domaincaller)


### parameters
```
$ domaincaller --help


usage: domaincaller <--uri cool -O output> [options]

Detect minimum domains using adaptive DI

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Print version number and exit
  --uri URI             Cool URI. (default: None)
  -O OUTPUT, --output OUTPUT
  -D DI_OUTPUT, --DI-output DI_OUTPUT
  -W WEIGHT_COL, --weight-col WEIGHT_COL
                        Name of the column in .cool to be used to construct the normalized matrix. Specify "-W RAW" if you
                        want to run with the raw matrix. (default: weight)
  --exclude [EXCLUDE [EXCLUDE ...]]
                        List of chromosomes to exclude. (default: ['chrY', 'chrM'])
  -p CPU_CORE, --cpu-core CPU_CORE
                        Number of processes to launch. (default: 1)
  --removeCache         Remove cache data before exiting. (default: False)
  --logFile LOGFILE     Logging file name. (default: domaincaller.log)
(HiC-Pro) wuzhikun@mu03 17:02:21 ^_^ /home/wuzhikun/github/Primates3D 

```


```
(HiC-Pro) wuzhikun@mu03 17:20:10 O_O /home/wuzhikun/Project/Centromere/Cooler/HG02622/Pairs 
$ domaincaller  -O HG02622.nodups.maternal.tad --uri HG02622.nodups.maternal.cool
root                      INFO    @ 03/09/23 17:20:15: 
# ARGUMENT LIST:
# Output TAD file = HG02622.nodups.maternal.tad
# Output DI file = None
# Cool URI = HG02622.nodups.maternal.cool
# Column for matrix balancing = weight
# Excluded Chromosomes = ['chrY', 'chrM']
# Number of processes used = 1
# Remove cache data = False
# Log file name = domaincaller.log
numexpr.utils             INFO    @ 03/09/23 17:20:15: Note: NumExpr detected 24 cores but "NUMEXPR_MAX_THREADS" not set, so enforcing safe limit of 8.
numexpr.utils             INFO    @ 03/09/23 17:20:15: NumExpr defaulting to 8 threads.
root                      INFO    @ 03/09/23 17:20:16: Loading Hi-C matrix ...

```


```
(HiC-Pro) wuzhikun@fat01 16:22:42 ^_^ /home/wuzhikun/Project/Centromere/genome/Primates 
$ domaincaller  --uri test.hg38.1000000.mcool --output test.hg38.1000000.domain --cpu-core 20
root                      INFO    @ 02/04/23 16:23:18: 
# ARGUMENT LIST:
# Output TAD file = test.hg38.1000000.domain
# Output DI file = None
# Cool URI = test.hg38.1000000.mcool
# Window Size = 2000000
# Column for matrix balancing = weight
# Excluded Chromosomes = ['chrY', 'chrM']
# Number of processes used = 20
# Log file name = domaincaller.log
numexpr.utils             INFO    @ 02/04/23 16:23:19: Note: NumExpr detected 64 cores but "NUMEXPR_MAX_THREADS" not set, so enforcing safe limit of 8.
numexpr.utils             INFO    @ 02/04/23 16:23:19: NumExpr defaulting to 8 threads.

```





```
(HiC-Pro) wuzhikun@fat01 10:47:31 ^_^ /home/wuzhikun/Project/Centromere/Cooler/Arabidopsis_thaliana/Pairs 
$ domaincaller  --uri Arabidopsis_thaliana.nodups_copy.cool --output test.hg38.1000000.domain --cpu-core 20
root                      INFO    @ 02/06/23 10:53:09: 
# ARGUMENT LIST:
# Output TAD file = test.hg38.1000000.domain
# Output DI file = None
# Cool URI = Arabidopsis_thaliana.nodups_copy.cool
# Column for matrix balancing = weight
# Excluded Chromosomes = ['chrY', 'chrM']
# Number of processes used = 20
# Remove cache data = False
# Log file name = domaincaller.log
numexpr.utils             INFO    @ 02/06/23 10:53:09: Note: NumExpr detected 64 cores but "NUMEXPR_MAX_THREADS" not set, so enforcing safe limit of 8.
numexpr.utils             INFO    @ 02/06/23 10:53:09: NumExpr defaulting to 8 threads.
root                      INFO    @ 02/06/23 10:53:34: Loading Hi-C matrix ...
root                      INFO    @ 02/06/23 10:53:55: Done!
/home/wuzhikun/anaconda3/envs/HiC-Pro/lib/python3.8/site-packages/tadlib/domaincaller/chromLev.py:211: RuntimeWarning: invalid value encountered in long_scalars
  noiselevel = check.sum() / check.size

```



