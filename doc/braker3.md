## [braker](https://github.com/Gaius-Augustus/BRAKER)


### install 


```
mamba install -c bioconda braker3
```

install information
``
The config/ directory from AUGUSTUS can be accessed with the variable AUGUSTUS_CONFIG_PATH.
BRAKER3 requires this directory to be in a writable location, so if that is not the case, copy this directory to a writable location, e.g.:
cp -r /home/wuzhikun/anaconda3/envs/braker3/config/ /absolute_path_to_user_writable_directory/
export AUGUSTUS_CONFIG_PATH=/absolute_path_to_user_writable_directory/config

Due to license and distribution restrictions, GeneMark-ETP and ProtHint should be additionally installed for BRAKER3 to fully work.
These packages can be either installed as part of the BRAKER3 environment, or the PATH variable should be configured to point to them.
The GeneMark key should be located in /home/wuzhikun/.gm_key and GENEMARK_PATH should include the path to the GeneMark executables gmes_petap.pl or gmetp.pl.
```


export AUGUSTUS_CONFIG_PATH=/home/wuzhikun/anaconda3/envs/braker3/config

--addUTR=on


[GeneMark-ETP](https://github.com/gatech-genemark/GeneMark-ETP)

[GeneMark key](http://topaz.gatech.edu/GeneMark/license_download.cgi)

[ProtHint](https://github.com/gatech-genemark/ProtHint)

[GUSHR](https://github.com/Gaius-Augustus/GUSHR)



### run braker3


```
(braker3) wuzhikun@cu17 15:38:40 ^_^ /home/wuzhikun/Project/Vigna/braker3 
$ braker.pl --threads=10 --genome /home/wuzhikun/Project/Vigna/Repeats/Final/repeatMasker_soft/Vigna_unguiculata_assembly.fasta.masked --prot_seq /home/wuzhikun/database/genome/Vigna_radiata/GCF_000741045.1_Vradiata_ver6_protein.faa --bam /home/wuzhikun/Project/Vigna/RNA/NGS/mapping/hisat/pod/ngs_RNA_hisat.bam --skipOptimize --gff3 –-workingdir /home/wuzhikun/Project/Vigna/braker3 --GUSHR_PATH /home/wuzhikun/anaconda3/envs/braker3/bin/GUSHR/gushr.py  --GENEMARK_PATH=/home/wuzhikun/anaconda3/envs/braker3/bin/GeneMark-ETP/bin  --PROTHINT_PATH=/home/wuzhikun/anaconda3/envs/braker3/bin/ProtHint/bin 

# Fri Aug  4 15:39:09 2023:Both protein and RNA-Seq data in input detected. BRAKER will be executed in ETP mode (BRAKER3).
#*********
# Fri Aug  4 15:39:13 2023: Log information is stored in file /home/wuzhikun/Project/Vigna/braker3/braker/braker.log
#*********
# WARNING: Detected whitespace in fasta header of file /home/wuzhikun/database/genome/Vigna_radiata/GCF_000741045.1_Vradiata_ver6_protein.faa. This may later on cause problems! The pipeline will create a new file without spaces or "|" characters and a genome_header.map file to look up the old and new headers. This message will be suppressed from now on!
#*********
Undefined subroutine &main::DumpFile called at /home/wuzhikun/anaconda3/envs/braker3/bin/braker.pl line 5440.
```






```
export PYTHONPATH=/home/wuzhikun/anaconda3/envs/braker3/lib:$PYTHONPATH
export PATH=/home/wuzhikun/anaconda3/envs/braker3/bin:$PATH
export PERL5LIB=/home/wuzhikun/anaconda3/envs/braker3/lib:$PERL5LIB


(braker3) wuzhikun@cu17 15:54:36 ^_^ /home/wuzhikun/Project/Vigna/braker3 
$ braker.pl --threads=10 --genome /home/wuzhikun/Project/Vigna/Repeats/Final/repeatMasker_soft/Vigna_unguiculata_assembly.fasta.masked --prot_seq /home/wuzhikun/Project/Vigna/braker3/TAIR10_pep_20101214.fa --bam /home/wuzhikun/Project/Vigna/RNA/NGS/mapping/hisat/pod/ngs_RNA_hisat.bam --skipOptimize --gff3 –-workingdir /home/wuzhikun/Project/Vigna/braker3 --GUSHR_PATH /home/wuzhikun/anaconda3/envs/braker3/bin/GUSHR  --GENEMARK_PATH=/home/wuzhikun/anaconda3/envs/braker3/bin/GeneMark-ETP/bin  --PROTHINT_PATH=/home/wuzhikun/anaconda3/envs/braker3/bin/ProtHint/bin 

```



