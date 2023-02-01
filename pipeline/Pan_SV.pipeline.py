#!/usr/bin/env python
# coding: utf-8
import os
import sys
import yaml


### plot workflow
### snakemake -s /home/wuzhikun/github/TrioWGS/pipeline/TrioWGS.pipeline.py  --configfile /home/wuzhikun/github/TrioWGS/pipeline/TrioWGS.pipeline.yaml -j 48 --dag | dot -Tsvg > TrioWGS.svg



IN_PATH = config['project_path']
THREADS = config['THREADS']
SAMPLES = config['SAMPLES']
SCAFFOLDS = config['SCAFFOLDS']
ThreadFold = config['ThreadFold']

# ENV_BIN = config["ENV_BIN"]
SCRIPT_DIR = config["SCRIPT_DIR"]
PIPE_DIR = config["PIPE_DIR"]
RULE_DIR = config["RULE_DIR"]
SRC_DIR = config["SRC_DIR"]


# include: RULE_DIR + '/BaseFunction.rule.py'
include: RULE_DIR + '/QualityControl.rule.py'
# include: RULE_DIR + '/Mapping.rule.py'
include: RULE_DIR + '/MultiAssembly.rule.py'
include: RULE_DIR + '/GenePrediction.rule.py'
include: RULE_DIR + '/StructuralVariant.rule.py'

rule all:
    input:
        ################################# Quality control ############################
        # expand(IN_PATH + '/FastQC/clean/{sample}_stats.xls', sample=SAMPLES),
        # ################################ mapping ####################
        # expand(IN_PATH + '/mapping/{sample}/{sample}_align.bam', sample=SAMPLES),
        # expand(IN_PATH + '/mapping/{sample}/{sample}_bam_stats.xls', sample=SAMPLES),
        # expand(IN_PATH + "/SyRI/{scaffold}/SyRI/{scaffold}_aln_out.txt", scaffold=SCAFFOLDS),
        # expand(IN_PATH + "/clean/{sample}_ONT.fastq.gz", sample=SAMPLES),
        # expand(IN_PATH + '/clean/{sample}.RNA.R1.fq.gz', sample=SAMPLES),
        # expand(IN_PATH + "/Assembly/NextDenovo/{sample}/03.ctg_graph/nd.asm.fasta", sample=SAMPLES),
        # expand(IN_PATH + "/Assembly/Polish/{sample}/Racon/{sample}_ONT_polish_racon3.fasta", sample=SAMPLES),
        expand(IN_PATH + "/Assembly/Polish/{sample}/Pilon/{sample}.SRS_pilon.fasta", sample=SAMPLES),
        # expand(IN_PATH + "/Assembly/Scaffold/{sample}/ragtag.scaffold.fasta", sample=SAMPLES),
        expand(IN_PATH + "/Evaluation/QV/{sample}/{sample}.qv", sample=SAMPLES),
        # expand(IN_PATH + "/Evaluation/BUSCO/{sample}/run_embryophyta_odb10/short_summary.txt", sample=SAMPLES),
        # expand(IN_PATH + "/SVCall/cuteSV/{sample}.cutesv.vcf", sample=SAMPLES),
        # expand(IN_PATH + "/GenePrediction/RNA/Stringtie/{sample}/{sample}_RNA_stringtie.gtf", sample=SAMPLES),
        # expand(IN_PATH + "/GenePrediction/RNA/trinity/{sample}/trinity/Trinity.fasta", sample=SAMPLES),
        # expand(IN_PATH + "/Evaluation/Coverage/{sample}/{sample}_ONT_WGS.bw", sample=SAMPLES),
        # expand(IN_PATH +  "Evaluation/Coverage/{sample}/{sample}_ONT.per-base.bed.gz", sample=SAMPLES),
        # expand(IN_PATH + "/Repeat/repeatMasker/{sample}/{sample}.fasta.out", sample=SAMPLES),
        expand(IN_PATH + "/GenePrediction/RNA/TransMerge/{sample}/nonreduandant_transcript.fasta", sample=SAMPLES),
        expand(IN_PATH + "/GenePrediction/snap/{sample}/{sample}_transcript.fa", sample=SAMPLES),
        # expand(IN_PATH + "/GenePrediction/augustus/{sample}/{sample}_prediction.gff", sample=SAMPLES),
        IN_PATH + "/GenePrediction/Protein/Merged_protein_nonredundant.fasta",