## [cactus](https://github.com/ComparativeGenomicsToolkit/cactus)


### [cactus pangenome](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md)

### [cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/c6637f8904d84e7ef4bb09e5887c63b6fe63b158/doc/pangenome.md)


### install cactus

```
wegt https://github.com/ComparativeGenomicsToolkit/cactus/releases/download/v2.4.0/cactus-bin-v2.4.0.tar.gz
```


```
cd /home/wuzhikun/software/cactus-bin-v2.4.0/bin
virtualenv -p python3.8 cactus_env
echo "export PATH=/home/wuzhikun/software/cactus-bin-v2.4.0/bin:\$PATH" >> cactus_env/bin/activate
echo "export PYTHONPATH=/home/wuzhikun/software/cactus-bin-v2.4.0/lib:\$PYTHONPATH" >> cactus_env/bin/activate
source cactus_env/bin/activate
python3 -m pip install -U setuptools pip==21.3.1
python3 -m pip install -U -r /home/wuzhikun/software/cactus-bin-v2.4.0/toil-requirement.txt
python3 -m pip install -U /home/wuzhikun/software/cactus-bin-v2.4.0
```

```
$ source /home/wuzhikun/cactus_env/bin/activate


$ which cactus-minigraph 
~/cactus_env/bin/cactus-minigraph


$ ls /home/wuzhikun/cactus_env/bin/cactus*

/home/wuzhikun/cactus_env/bin/cactus
/home/wuzhikun/cactus_env/bin/cactus-align
/home/wuzhikun/cactus_env/bin/cactus-align-batch
/home/wuzhikun/cactus_env/bin/cactus-blast
/home/wuzhikun/cactus_env/bin/cactus-graphmap
/home/wuzhikun/cactus_env/bin/cactus-graphmap-join
/home/wuzhikun/cactus_env/bin/cactus-graphmap-split
/home/wuzhikun/cactus_env/bin/cactus-hal2maf
/home/wuzhikun/cactus_env/bin/cactus-minigraph
/home/wuzhikun/cactus_env/bin/cactus-prepare
/home/wuzhikun/cactus_env/bin/cactus-prepare-toil
/home/wuzhikun/cactus_env/bin/cactus-preprocess
/home/wuzhikun/cactus_env/bin/cactus-refmap
/home/wuzhikun/cactus_env/bin/cactus-terra-helper
/home/wuzhikun/cactus_env/bin/cactus-update-prepare

```

### Parameters

```
(cactus_env) wuzhikun@cu26 15:53:22 ^_^ /home/wuzhikun/Project/PanSV/Cactus 
$ /home/wuzhikun/cactus_env/bin/cactus-minigraph --help
usage: cactus-minigraph [-h] [--logCritical] [--logError] [--logWarning] [--logDebug] [--logInfo] [--logOff]
                        [--logLevel {Critical,Error,Warning,Debug,Info,critical,error,warning,debug,info,CRITICAL,ERROR,WARNING,DEBUG,INFO}]
                        [--logFile LOGFILE] [--rotatingLogging] [--workDir WORKDIR]
                        [--coordinationDir COORDINATION_DIR] [--noStdOutErr] [--stats]
                        [--clean {always,onError,never,onSuccess}]
                        [--cleanWorkDir {always,onError,never,onSuccess}] [--clusterStats [CLUSTERSTATS]]
                        [--restart]
                        [--batchSystem {aws_batch,parasol,single_machine,grid_engine,lsf,mesos,slurm,tes,torque,htcondor,kubernetes}]
                        [--disableHotDeployment] [--disableAutoDeployment] [--maxLocalJobs MAXLOCALJOBS]
                        [--manualMemArgs] [--runCwlInternalJobsOnWorkers] [--coalesceStatusCalls]
                        [--statePollingWait STATEPOLLINGWAIT] [--awsBatchRegion AWS_BATCH_REGION]
                        [--awsBatchQueue AWS_BATCH_QUEUE] [--awsBatchJobRoleArn AWS_BATCH_JOB_ROLE_ARN]
                        [--parasolCommand PARASOLCOMMAND] [--parasolMaxBatches PARASOLMAXBATCHES]
                        [--scale SCALE] [--dont_allocate_mem | --allocate_mem] [--tesEndpoint TES_ENDPOINT]
                        [--tesUser TES_USER] [--tesPassword TES_PASSWORD] [--tesBearerToken TES_BEARER_TOKEN]
                        [--linkImports | --noLinkImports] [--moveExports | --noMoveExports] [--disableCaching]
                        [--provisioner {aws,gce,None}] [--nodeTypes NODETYPES] [--minNodes MINNODES]
                        [--maxNodes MAXNODES] [--targetTime TARGETTIME] [--betaInertia BETAINERTIA]
                        [--scaleInterval SCALEINTERVAL] [--preemptableCompensation PREEMPTABLECOMPENSATION]
                        [--nodeStorage NODESTORAGE] [--nodeStorageOverrides NODESTORAGEOVERRIDES] [--metrics]
                        [--assumeZeroOverhead] [--maxServiceJobs MAXSERVICEJOBS]
                        [--maxPreemptableServiceJobs MAXPREEMPTABLESERVICEJOBS] [--deadlockWait DEADLOCKWAIT]
                        [--deadlockCheckInterval DEADLOCKCHECKINTERVAL] [--defaultMemory INT]
                        [--defaultCores FLOAT] [--defaultDisk INT]
                        [--defaultAccelerators ACCELERATOR[,ACCELERATOR...]] [--defaultPreemptable [BOOL]]
                        [--maxCores INT] [--maxMemory INT] [--maxDisk INT] [--retryCount RETRYCOUNT]
                        [--enableUnlimitedPreemptableRetries] [--doubleMem] [--maxJobDuration MAXJOBDURATION]
                        [--rescueJobsFrequency RESCUEJOBSFREQUENCY] [--maxLogFileSize MAXLOGFILESIZE]
                        [--writeLogs [WRITELOGS]] [--writeLogsGzip [WRITELOGSGZIP]] [--writeLogsFromAllJobs]
                        [--writeMessages WRITE_MESSAGES] [--realTimeLogging] [--disableChaining]
                        [--disableJobStoreChecksumVerification] [--sseKey SSEKEY] [--setEnv NAME=VALUE or NAME]
                        [--servicePollingInterval SERVICEPOLLINGINTERVAL] [--forceDockerAppliance]
                        [--statusWait STATUSWAIT] [--disableProgress] [--debugWorker]
                        [--disableWorkerOutputCapture] [--badWorker BADWORKER]
                        [--badWorkerFailInterval BADWORKERFAILINTERVAL] --reference REFERENCE
                        [--mapCores MAPCORES] [--configFile CONFIGFILE] [--latest]
                        [--containerImage CONTAINERIMAGE] [--binariesMode {docker,local,singularity}]
                        jobStore seqFile outputGFA

positional arguments:
  seqFile               Seq file (will be modified if necessary to include graph Fasta sequence)
  outputGFA             Output Minigraph GFA

optional arguments:
  -h, --help            show this help message and exit
  --reference REFERENCE
                        Reference genome name (added to minigraph first). Order in seqfile used otherwise
  --mapCores MAPCORES   Number of cores for minigraph. Overrides graphmap cpu in configuration
  --configFile CONFIGFILE
                        Specify cactus configuration file
  --latest              Use the latest version of the docker container rather than pulling one matching this
                        version of cactus
  --containerImage CONTAINERIMAGE
                        Use the the specified pre-built containter image rather than pulling one from quay.io
  --binariesMode {docker,local,singularity}
                        The way to run the Cactus binaries

```


### run cactus-minigraph

### Input files

```
$ seqkit grep -n -p  chr5 An-1.chr.all.v2.0.fasta > An-1.chr.all.v2.0_chr5.fasta
...

$ cat /home/wuzhikun/Project/PanSV/Cactus/sequence/seqfile.txt
An-1	/home/wuzhikun/Project/PanSV/Cactus/sequence/An-1.chr.all.v2.0_chr5.fasta
C24	/home/wuzhikun/Project/PanSV/Cactus/sequence/C24.chr.all.v2.0_chr5.fasta
Col	/home/wuzhikun/Project/PanSV/Cactus/sequence/Col-CEN_v1.2_chr5.fasta
Ler	/home/wuzhikun/Project/PanSV/Cactus/sequence/Ler.chr.all.v2.0_chr5.fasta


```



```
(cactus_env) (Assembly3) wuzhikun@fat01 09:23:10 ^_^ /home/wuzhikun/Project/PanSV/Cactus
$ /home/wuzhikun/cactus_env/bin/cactus-minigraph  /home/wuzhikun/Project/PanSV/Cactus/arabidopsis  /home/wuzhikun/Project/PanSV/Cactus/sequence/seqfile.txt  /home/wuzhikun/Project/PanSV/Cactus/arabidopsis_4lines.gfa  --reference Col --mapCores 40 


[2023-02-01T09:23:17+0800] [MainThread] [I] [toil.statsAndLogging] Enabling realtime logging in Toil
[2023-02-01T09:23:17+0800] [MainThread] [I] [toil.statsAndLogging] Cactus Command: /home/wuzhikun/cactus_env/bin/cactus-minigraph /home/wuzhikun/Project/PanSV/Cactus/arabidopsis /home/wuzhikun/Project/PanSV/Cactus/sequence/seqfile.txt /home/wuzhikun/Project/PanSV/Cactus/arabidopsis_4lines.gfa --reference Col --mapCores 40
[2023-02-01T09:23:17+0800] [MainThread] [I] [toil.statsAndLogging] Cactus Commit: 47f9079cc31a5533ffb76f038480fdec1b6f7c4f
[2023-02-01T09:23:17+0800] [MainThread] [I] [toil.statsAndLogging] Importing file:///home/wuzhikun/Project/PanSV/Cactus/sequence/An-1.chr.all.v2.0_chr5.fasta
[2023-02-01T09:23:23+0800] [MainThread] [I] [toil.statsAndLogging] Importing file:///home/wuzhikun/Project/PanSV/Cactus/sequence/C24.chr.all.v2.0_chr5.fasta
[2023-02-01T09:23:23+0800] [MainThread] [I] [toil.statsAndLogging] Importing file:///home/wuzhikun/Project/PanSV/Cactus/sequence/Col-CEN_v1.2_chr5.fasta
[2023-02-01T09:23:23+0800] [MainThread] [I] [toil.statsAndLogging] Importing file:///home/wuzhikun/Project/PanSV/Cactus/sequence/Ler.chr.all.v2.0_chr5.fasta
Traceback (most recent call last):
  File "/home/wuzhikun/cactus_env/bin/cactus-minigraph", line 8, in <module>
    sys.exit(main())
  File "/home/wuzhikun/cactus_env/lib/python3.8/site-packages/cactus/refmap/cactus_minigraph.py", line 122, in main
    gfa_id = toil.start(Job.wrapJobFn(minigraph_construct_workflow, config_node, input_seq_id_map, input_seq_order, options.outputGFA))
  File "/home/wuzhikun/cactus_env/lib/python3.8/site-packages/toil/common.py", line 996, in start
    self._batchSystem = self.createBatchSystem(self.config)
  File "/home/wuzhikun/cactus_env/lib/python3.8/site-packages/toil/common.py", line 1144, in createBatchSystem
    return batch_system(**kwargs)
  File "/home/wuzhikun/cactus_env/lib/python3.8/site-packages/toil/batchSystems/singleMachine.py", line 149, in __init__
    self.accelerator_identities = get_individual_local_accelerators()
  File "/home/wuzhikun/cactus_env/lib/python3.8/site-packages/toil/lib/accelerators.py", line 84, in get_individual_local_accelerators
    return [{'kind': 'gpu', 'brand': 'nvidia', 'api': 'cuda', 'count': 1} for _ in range(count_nvidia_gpus())]
  File "/home/wuzhikun/cactus_env/lib/python3.8/site-packages/toil/lib/accelerators.py", line 64, in count_nvidia_gpus
    subprocess.check_output(["nvidia-smi", "-q", "-x"])
  File "/home/wuzhikun/anaconda3/envs/Assembly3/lib/python3.8/subprocess.py", line 411, in check_output
    return run(*popenargs, stdout=PIPE, timeout=timeout, check=True,
  File "/home/wuzhikun/anaconda3/envs/Assembly3/lib/python3.8/subprocess.py", line 489, in run
    with Popen(*popenargs, **kwargs) as process:
  File "/home/wuzhikun/anaconda3/envs/Assembly3/lib/python3.8/subprocess.py", line 854, in __init__
    self._execute_child(args, executable, preexec_fn, close_fds,
  File "/home/wuzhikun/anaconda3/envs/Assembly3/lib/python3.8/subprocess.py", line 1702, in _execute_child
    raise child_exception_type(errno_num, err_msg, err_filename)
PermissionError: [Errno 13] Permission denied: 'nvidia-smi'
```


It seens need GPU


```
$ qsub -I -l nodes=cu26 -q gpu
$ nvidia-smi
Wed Feb  1 10:22:22 2023       
+-----------------------------------------------------------------------------+
| NVIDIA-SMI 515.43.04    Driver Version: 515.43.04    CUDA Version: 11.7     |
|-------------------------------+----------------------+----------------------+
| GPU  Name        Persistence-M| Bus-Id        Disp.A | Volatile Uncorr. ECC |
| Fan  Temp  Perf  Pwr:Usage/Cap|         Memory-Usage | GPU-Util  Compute M. |
|                               |                      |               MIG M. |
|===============================+======================+======================|
|   0  Tesla V100-PCIE...  On   | 00000000:84:00.0 Off |                    0 |
| N/A   33C    P0    27W / 250W |      0MiB / 32768MiB |      0%      Default |
|                               |                      |                  N/A |
+-------------------------------+----------------------+----------------------+
                                                                               
+-----------------------------------------------------------------------------+
| Processes:                                                                  |
|  GPU   GI   CI        PID   Type   Process name                  GPU Memory |
|        ID   ID                                                   Usage      |
|=============================================================================|
|  No running processes found                                                 |
+-----------------------------------------------------------------------------+
```

add parameters

```
--defaultAccelerators gpu
```


### pipeline

#### minigraph

```
(cactus_env) (Assembly3) wuzhikun@fat01 09:27:42 ^_^ /home/wuzhikun/Project/PanSV/Cactus 
$ minigraph -x ggs -t 40  -o /home/wuzhikun/Project/PanSV/Cactus/arabidopsis.gfa  /home/wuzhikun/Project/PanSV/Cactus/sequence/Col-CEN_v1.2_chr5.fasta /home/wuzhikun/Project/PanSV/Cactus/sequence/An-1.chr.all.v2.0_chr5.fasta /home/wuzhikun/Project/PanSV/Cactus/sequence/C24.chr.all.v2.0_chr5.fasta /home/wuzhikun/Project/PanSV/Cactus/sequence/Ler.chr.all.v2.0_chr5.fasta  

[WARNING] it is recommended to add -c for graph generation
[M::main::0.102*0.97] loaded the graph from "/home/wuzhikun/Project/PanSV/Cactus/sequence/Col-CEN_v1.2_chr5.fasta"
[M::mg_index::1.563*1.69] indexed the graph
[M::mg_opt_update::1.668*1.65] occ_max1=63; lc_max_occ=2
[M::ggen_map::1.788*1.60] loaded file "/home/wuzhikun/Project/PanSV/Cactus/sequence/An-1.chr.all.v2.0_chr5.fasta"
[M::ggen_map::20.108*1.05] mapped 1 sequence(s) to the graph
[M::mg_ggsimple::21.136*1.05] inserted 2014 events, including 1 inversions
[M::mg_index::22.614*1.12] indexed the graph
[M::mg_opt_update::22.725*1.11] occ_max1=63; lc_max_occ=2
[M::ggen_map::22.817*1.11] loaded file "/home/wuzhikun/Project/PanSV/Cactus/sequence/C24.chr.all.v2.0_chr5.fasta"
[M::ggen_map::121.232*1.02] mapped 1 sequence(s) to the graph
[M::mg_ggsimple::122.302*1.02] inserted 1503 events, including 1 inversions
[M::mg_index::123.956*1.04] indexed the graph
[M::mg_opt_update::124.059*1.04] occ_max1=63; lc_max_occ=2
[M::ggen_map::124.173*1.04] loaded file "/home/wuzhikun/Project/PanSV/Cactus/sequence/Ler.chr.all.v2.0_chr5.fasta"
[M::ggen_map::228.575*1.02] mapped 1 sequence(s) to the graph
[M::mg_ggsimple::229.366*1.02] inserted 1054 events, including 0 inversions
[M::main] Version: 0.20-r559
[M::main] CMD: minigraph -x ggs -t 40 -o /home/wuzhikun/Project/PanSV/Cactus/arabidopsis.gfa /home/wuzhikun/Project/PanSV/Cactus/sequence/Col-CEN_v1.2_chr5.fasta /home/wuzhikun/Project/PanSV/Cactus/sequence/An-1.chr.all.v2.0_chr5.fasta /home/wuzhikun/Project/PanSV/Cactus/sequence/C24.chr.all.v2.0_chr5.fasta /home/wuzhikun/Project/PanSV/Cactus/sequence/Ler.chr.all.v2.0_chr5.fasta
[M::main] Real time: 229.450 sec; CPU: 234.330 sec; Peak RSS: 0.880 GB


```


#### cactus-graphmap


```
(cactus_env) wuzhikun@cu26 10:35:58 ^_^ /home/wuzhikun/Project/PanSV/Cactus 
$ /home/wuzhikun/cactus_env/bin/cactus-graphmap  /home/wuzhikun/Project/PanSV/Cactus/arabidopsis_graphmap  /home/wuzhikun/Project/PanSV/Cactus/sequence/seqfile.txt  /home/wuzhikun/Project/PanSV/Cactus/arabidopsis.gfa  /home/wuzhikun/Project/PanSV/Cactus/arabidopsis.paf  --outputFasta  /home/wuzhikun/Project/PanSV/Cactus/arabidopsis.gfa.fa --reference Col --mapCores 40 --realTimeLogging  --defaultAccelerators gpu
```

output files:
```
-rw-rw-r-- 1 wuzhikun wuzhikun  7.3M Feb  1 10:31 arabidopsis.paf
-rw-rw-r-- 1 wuzhikun wuzhikun 1012K Feb  1 10:31 arabidopsis.gaf.gz
-rw-rw-r-- 1 wuzhikun wuzhikun   34M Feb  1 10:31 arabidopsis.gfa.fa

```


#### cactus-preprocess

need mask file

```
$ mkdir mask
$ cat  /home/wuzhikun/Project/PanSV/Cactus/sequence/seqfile.txt | sed 's/sequence/mask/g' > /home/wuzhikun/Project/PanSV/Cactus/sequence/seqfile_mask.txt


$ cat /home/wuzhikun/Project/PanSV/Cactus/sequence/seqfile_mask.txt
(An-1:1.0,C24:1.0,Col:1.0,Ler:1.0,_MINIGRAPH_:1.0);
An-1    /home/wuzhikun/Project/PanSV/Cactus/mask/An-1.chr.all.v2.0_chr5.fasta
C24 /home/wuzhikun/Project/PanSV/Cactus/mask/C24.chr.all.v2.0_chr5.fasta
Col /home/wuzhikun/Project/PanSV/Cactus/mask/Col-CEN_v1.2_chr5.fasta
Ler /home/wuzhikun/Project/PanSV/Cactus/mask/Ler.chr.all.v2.0_chr5.fasta
_MINIGRAPH_ file:///home/wuzhikun/Project/PanSV/Cactus/arabidopsis.gfa.fa
```

run cactus-preprocess

```
(cactus_env) wuzhikun@cu26 11:09:17 O_O /home/wuzhikun/Project/PanSV/Cactus 
$ /home/wuzhikun/cactus_env/bin/cactus-preprocess   --logFile preprocess.log    --maskFile /home/wuzhikun/Project/PanSV/Cactus/arabidopsis.paf  --minLength 100000  --maskMode brnn  --brnnCores 20  /home/wuzhikun/Project/PanSV/Cactus/arabidopsis_preprocess  /home/wuzhikun/Project/PanSV/Cactus/sequence/seqfile.txt /home/wuzhikun/Project/PanSV/Cactus/sequence/seqfile_mask.txt

```



output files

```
drwxrwxr-x 1 wuzhikun wuzhikun  4.0K Feb  1 11:35 mask
-rw-rw-r-- 1 wuzhikun wuzhikun   13K Feb  1 11:35 preprocess.log


$ tree mask/
mask/
├── An-1.chr.all.v2.0_chr5.fasta
├── An-1.chr.all.v2.0_chr5.fasta.bed
├── An-1.chr.all.v2.0_chr5.fasta.mask.bed
├── C24.chr.all.v2.0_chr5.fasta
├── C24.chr.all.v2.0_chr5.fasta.bed
├── C24.chr.all.v2.0_chr5.fasta.mask.bed
├── Col-CEN_v1.2_chr5.fasta
├── Col-CEN_v1.2_chr5.fasta.bed
├── Col-CEN_v1.2_chr5.fasta.mask.bed
├── Ler.chr.all.v2.0_chr5.fasta
├── Ler.chr.all.v2.0_chr5.fasta.bed
└── Ler.chr.all.v2.0_chr5.fasta.mask.bed

```


#### cactus-graphmap


```
/home/wuzhikun/cactus_env/bin/cactus-graphmap  --maskFilter 100000 --logFile graphmap2.log --outputFasta /home/wuzhikun/Project/PanSV/Cactus/arabidopsis.align.gfa.fa   /home/wuzhikun/Project/PanSV/Cactus/arabidopsis_align /home/wuzhikun/Project/PanSV/Cactus/sequence/seqfile_mask.txt /home/wuzhikun/Project/PanSV/Cactus/arabidopsis.gfa /home/wuzhikun/Project/PanSV/Cactus/arabidopsis.align.paf 
```

output files

```
-rw-rw-r-- 1 wuzhikun wuzhikun  7.3M Feb  1 13:14 arabidopsis.align.paf
-rw-rw-r-- 1 wuzhikun wuzhikun 1012K Feb  1 13:14 arabidopsis.align.gaf.gz
-rw-rw-r-- 1 wuzhikun wuzhikun   34M Feb  1 13:14 arabidopsis.align.gfa.fa
-rw-rw-r-- 1 wuzhikun wuzhikun   17K Feb  1 13:14 graphmap2.log

```

#### cactus-align

```
(cactus_env) wuzhikun@cu26 13:49:51 ^_^ /home/wuzhikun/Project/PanSV/Cactus 
$ /home/wuzhikun/cactus_env/bin/cactus-align   --reference Col --consCores 20 --gpu --pangenome --outVG  /home/wuzhikun/Project/PanSV/Cactus/cactus_align2 /home/wuzhikun/Project/PanSV/Cactus/sequence/seqfile_mask.txt  /home/wuzhikun/Project/PanSV/Cactus/arabidopsis.align.paf  /home/wuzhikun/Project/PanSV/Cactus/arabidopsis.align2.hal
```



output files:
```
-rw-rw-r-- 1 wuzhikun wuzhikun   65M Feb  1 13:45 arabidopsis.align2.hal
-rw-rw-r-- 1 wuzhikun wuzhikun   70M Feb  1 13:45 arabidopsis.align2.vg

```

#### cactus-graphmap-join


```
(cactus_env) wuzhikun@cu26 14:21:51 ^_^ /home/wuzhikun/Project/PanSV/Cactus 
$ /home/wuzhikun/cactus_env/bin/cactus-graphmap-join /home/wuzhikun/Project/PanSV/Cactus/graphmap_join  --vg /home/wuzhikun/Project/PanSV/Cactus/arabidopsis.align2.vg --hal /home/wuzhikun/Project/PanSV/Cactus/arabidopsis.align2.hal  --outDir /home/wuzhikun/Project/PanSV/Cactus/arabidopsis_pg --outName arabidopsis --reference Col --vcf --giraffe
```



output files

```
drwxrwxr-x 1 wuzhikun wuzhikun  4.0K Feb  1 14:17 arabidopsis_pg


$ tree arabidopsis_pg/
arabidopsis_pg/
├── arabidopsis.d2.dist
├── arabidopsis.d2.gbz
├── arabidopsis.d2.min
├── arabidopsis.full.hal
├── arabidopsis.gfa.gz
├── arabidopsis.raw.vcf.gz
├── arabidopsis.raw.vcf.gz.tbi
├── arabidopsis.stats.tgz
├── arabidopsis.vcf.gz
└── arabidopsis.vcf.gz.tbi

```


