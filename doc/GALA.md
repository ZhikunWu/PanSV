
## [GALA](https://github.com/ganlab/GALA)

Long-reads Gap-free Chromosome-scale Assembler

### install
```
git clone https://github.com/ganlab/GALA.git
```





```
(Assembly) wuzhikun@cu22 19:41:51 ^_^ /home/wuzhikun/Project/Vigna/GALA 
$ /home/wuzhikun/software/GALA/gala /home/wuzhikun/Project/Vigna/GALA/draft_names_paths.txt  fa /home/wuzhikun/data/Vigna/HiFi/m64181_210607_093522.ccs.fasta   pacbio-corrected -f Vigna -o /home/wuzhikun/Project/Vigna/GALA 
```


### run step by step

```
(Assembly) wuzhikun@cu22 08:40:52 ^_^ /home/wuzhikun/Project/Vigna/GALA/scaffolds 

$ /home/wuzhikun/software/GALA/comp draft_names_paths.txt
```
```
$ sh draft_comp.sh 
```
out files:
```
comparison
├── draft_01vsdraft_02.paf
├── draft_01vsdraft_03.paf
├── draft_02vsdraft_01.paf
├── draft_02vsdraft_03.paf
├── draft_03vsdraft_01.paf
└── draft_03vsdraft_02.paf

```


$ /home/wuzhikun/software/GALA/ccm comparison 3
```

out files:
```
scaffolds
├── scaffolds_draft_01.scaff
├── scaffolds_draft_02.scaff
└── scaffolds_draft_03.scaff

```


