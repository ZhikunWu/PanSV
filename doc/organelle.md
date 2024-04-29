
### [GetOrganelle](https://github.com/Kinggerm/GetOrganelleDB?tab=readme-ov-file)

### download database
```
curl -L https://github.com/Kinggerm/GetOrganelleDB/releases/download/0.0.1/v0.0.1.tar.gz | tar zx
get_organelle_config.py -a embplant_pt,embplant_mt --use-local ./0.0.1
```

### run
```
(mitohifi) wuzhikun@cu12 11:12:12 ^_^ /home/wuzhikun/Project/BAssembly/pipeline/organelle 
$ nohup get_organelle_from_reads.py  -1 /home/wuzhikun/Project/BAssembly/clean/Bo.NGS.R1.fastq -2 /home/wuzhikun/Project/BAssembly/clean/Bo.NGS.R2.fastq -t 16 -o Boleracea_plastome  -F embplant_pt -R 10 > organelle.log 2>&1 &

```


out files:
```
drwxrwxr-x 1 wuzhikun wuzhikun 4.0K Apr 11 11:27 seed
-rw-rw-r-- 1 wuzhikun wuzhikun  78M Apr 11 11:33 extended_1_paired.fq
-rw-rw-r-- 1 wuzhikun wuzhikun  78M Apr 11 11:33 extended_2_paired.fq
-rw-rw-r-- 1 wuzhikun wuzhikun 1.9M Apr 11 11:33 extended_1_unpaired.fq
-rw-rw-r-- 1 wuzhikun wuzhikun 1.4M Apr 11 11:33 extended_2_unpaired.fq
drwxrwxr-x 1 wuzhikun wuzhikun 4.0K Apr 11 11:35 extended_spades
-rw-rw-r-- 1 wuzhikun wuzhikun 150K Apr 11 11:35 embplant_pt.K115.scaffolds.graph1.1.path_sequence.fasta
-rw-rw-r-- 1 wuzhikun wuzhikun 150K Apr 11 11:35 embplant_pt.K115.scaffolds.graph1.2.path_sequence.fasta
-rw-rw-r-- 1 wuzhikun wuzhikun 126K Apr 11 11:35 embplant_pt.K115.contigs.graph1.selected_graph.gfa
-rw-rw-r-- 1 wuzhikun wuzhikun 323K Apr 11 11:35 extended_K115.assembly_graph.fastg
-rw-rw-r-- 1 wuzhikun wuzhikun 323K Apr 11 11:35 extended_K115.assembly_graph.fastg.extend-embplant_pt-embplant_mt.fastg
-rw-rw-r-- 1 wuzhikun wuzhikun 6.5K Apr 11 11:35 extended_K115.assembly_graph.fastg.extend-embplant_pt-embplant_mt.csv
-rw-rw-r-- 1 wuzhikun wuzhikun 8.9K Apr 11 11:35 get_org.log.txt

```


