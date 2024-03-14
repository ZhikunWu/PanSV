

### [基于PASA进行基因预测](https://www.cnblogs.com/zhanmaomao/p/12456073.html)

### 配置PASA config

```

## 配置
cd pasa_conf
cp pasa.CONFIG.template conf.txt
vi conf.txt

## 需要修改如下内容：
MYSQL_RW_USER=shehb
MYSQL_RW_PASSWORD=123456
MYSQL_RO_USER=pasa
MYSQL_RO_PASSWORD=123456
MYSQLSERVER=localhost  此处不能填写IP
PASA_ADMIN_EMAIL=邮箱
BASE_PASA_URL=http://pasa-dev.tigr.org/cgi-bin/
```

### 修改pasa.alignAssembly.Template.txt

```
cd pasa_conf
cp pasa.alignAssembly.Template.txt alignAssembly.config
vi alignAssembly.config

DATABASE=/tem/mydb.sqlite
validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=80
validate_alignments_in_db.dbi:--MIN_AVG_PER_ID=80
```



### modify


```
$ grep -n  dev /home/wuzhikun/software/PASApipeline-v2.5.1/scripts/process_PBLAT_alignments.pl
91:my $ooc_cmd = "$blat_path $genome_db $transcript_db -q=rna -dots=100 -maxIntron=$MAX_INTRON -threads=$CPU  -makeOoc=11.ooc /dev/null";

```

change **/dev/null** to **~/dev/null**



### library
```
export PERL5LIB=/home/wuzhikun/anaconda3/envs/Perl5/lib

(Perl5) wuzhikun@fat02 11:16:15 ^_^ /home/wuzhikun/Project/PanSV/pipeline/PASA 
$ Launch_PASA_pipeline.pl -c /home/wuzhikun/Project/PanSV/pipeline/PASA/alignAssembly.config  -C -R -g /home/wuzhikun/Project/PanSV/Repeat/repeatMasker/CPG02/CPG02.scaffold.fasta.masked -t /home/wuzhikun/Project/PanSV/GenePrediction/RNA/TransMerge/CPG02/nonreduandant_transcript.fasta  --ALIGNERS blat  --CPU 20

```


### run PASA
```

$ Launch_PASA_pipeline.pl -c /home/wuzhikun/Project/PanSV/pipeline/PASA/alignAssembly.config  -C -R -g /home/wuzhikun/Project/PanSV/Repeat/repeatMasker/CPG02/CPG02.scaffold.fasta.masked -t /home/wuzhikun/Project/PanSV/GenePrediction/RNA/TransMerge/CPG02/nonreduandant_transcript.fasta  --ALIGNERS blat  --CPU 20
```



