

## [EVidenceModeler(EVM)流程做基因组注释简单小例子](https://www.jianshu.com/p/abb6e6629832)



### parameters
```
(Anno) wuzhikun@fat01 17:25:08 ^_^ /home/wuzhikun/anaconda3/envs/Anno/bin 
$ EVidenceModeler --help

################# Evidence Modeler ##############################
#
#  parameters:
#
#  Required:
# 
#  --sample_id <str>           sample_id (for naming outputs)
#  --genome <str>              genome sequence in fasta format
#  --weights <str>             weights for evidence types file. See documentation for formatting.
#  --gene_predictions <str>    gene predictions gff3 file
#
#  # partitioning info (required too) 
#  --segmentSize <str>          * :length of a single sequence for running EVM
#  --overlapSize  <str>         * :length of sequence overlap between segmented sequences
#
#
#  Optional but recommended:
#  --protein_alignments <str>         protein alignments gff3 file
#  --transcript_alignments <str>      transcript alignments gff3 file
#
#  Optional and miscellaneous
#
#  --repeats <str>              gff3 file with repeats masked from genome file
#
#  
#  --terminalExons <str>        supplementary file of additional terminal exons to consider (from PASA long-orfs)
#
#  --stop_codons <str>            list of stop codons (default: TAA,TGA,TAG)
#                                 *for Tetrahymena, set --stop_codons TGA
#
#  --min_intron_length <int>      minimum length for an intron (default 20 bp)
#  --exec_dir <str>               directory that EVM cds to before running.
#
#  --CPU <int>                   number of parallel computes (default: 4)
#
#  --search_long_introns  <int>  when set, reexamines long introns (can find nested genes, but also can result in FPs) (default: 0 (off))
#
#
#  --re_search_intergenic <int>  when set, reexamines intergenic regions of minimum length (can add FPs) (default: 0  (off))
#  --terminal_intergenic_re_search <int>   reexamines genomic regions outside of the span of all predicted genes (default: 10000)
#
# flags:
#
#  --forwardStrandOnly   runs only on the forward strand
#  --reverseStrandOnly   runs only on the reverse strand
#
#  -S                    verbose flag
#  --debug               debug mode, writes lots of extra files.
#  --report_ELM          report the eliminated EVM preds too.
#
#  --version             report version (EVidenceModeler-v2.0.0) and exit.
#
#################################################################################################################################
```



time ~/anaconda3/envs/EVM/opt/evidencemodeler-2.1.0/EvmUtils/evidence_modeler.pl \
--genome ../../repeat/chr1.fa.masked \
--weights /data/myan/raw_data/practice/pan.genome/at.nc/version2.5.2019-10-09/genome.annotation/braker2/evm/weights.txt \
--gene_predictions evm_abinitio.gff3 \
--protein_alignments evm_pro.gff3 \
--transcript_alignments transcripts.fasta.transdecoder.genome.gff3 > evm.out


```
(Anno) wuzhikun@fat01 17:45:33 O_O /home/wuzhikun/Project/BAssembly/pipeline/august 
$ EVidenceModeler --sample_id test --segmentSize 100000 --overlapSize 10000  --genome /home/wuzhikun/Project/BAssembly/Assembly/hifiOnt100kHic/hifiOnt100kHic.hic.p_ctg_200k.fasta.split/hifiOnt100kHic.hic.p_ctg_200k.id_ptg000010l.fasta --weights weight.txt --gene_predictions hifiOnt100kHic.hic.p_ctg_200k.id_ptg000010l_august.gff3 > evm.out

```

out files:
```
-rw-rw-rw- 1 wuzhikun wuzhikun 4.6K Apr 12 17:45 test.partitions.listing
drwxrwxrwx 1 wuzhikun wuzhikun 4.0K Apr 12 17:45 test.partitions
-rw-rw-r-- 1 wuzhikun wuzhikun  25K Apr 12 17:45 test.partitions.evm_cmds
-rw-rw-r-- 1 wuzhikun wuzhikun  25K Apr 12 17:45 test.partitions.evm_cmds.completed
-rw-rw-r-- 1 wuzhikun wuzhikun  16K Apr 12 17:45 evm.out
-rw-rw-r-- 1 wuzhikun wuzhikun 1.1M Apr 12 17:46 test.EVM.gff3
-rw-rw-r-- 1 wuzhikun wuzhikun 436K Apr 12 17:46 test.EVM.pep
-rw-rw-r-- 1 wuzhikun wuzhikun 1.2M Apr 12 17:46 test.EVM.cds
-rw-rw-r-- 1 wuzhikun wuzhikun 2.4K Apr 12 17:46 __test-EVM_chckpts.cmds_log
-rw-rw-r-- 1 wuzhikun wuzhikun 131K Apr 12 17:46 test.EVM.bed
drwxrwxr-x 1 wuzhikun wuzhikun 4.0K Apr 12 17:46 __test-EVM_chckpts

```



```

nohup EVidenceModeler --sample_id CPG11 --segmentSize 100000 --overlapSize 10000  --genome /home/wuzhikun/Project/PanSV/Scaffold/CPG11.genome.fasta  --weights /home/wuzhikun/Project/BAssembly/pipeline/evidencemodeler/Trans_weight.txt --gene_predictions  /home/wuzhikun/Project/PanSV/GenePrediction/augustus/CPG11_augustus.gff --protein_alignments  /home/wuzhikun/Project/PanSV/GenePrediction/homo/miniprot/CPG11_align.gff --transcript_alignments  /home/wuzhikun/Project/PanSV/GenePrediction/RNA/Stringtie/CPG11/CPG11_RNA_stringtie.gtf  --repeats /home/wuzhikun/Project/PanSV/Repeats/repeatMasker/CPG11/CPG11.genome.fasta.out.gff --CPU 16  > /home/wuzhikun/Project/BAssembly/pipeline/evidencemodeler/CPG11.evm.out 2>CPG11.evm.out.log  &
```
