#!/usr/bin/env python
# coding: utf-8
import os
import sys
import yaml



### snakemake -s /home/wuzhikun/github/PanSV/pipeline/PanVigna.pipeline.py  --configfile /home/wuzhikun/github/PanSV/pipeline/PanVigan.pipeline.yaml -j 48 --nolock --dryrun



IN_PATH = config['project_path']
THREADS = config['THREADS']
SAMPLES = config['SAMPLES']
ThreadFold = config['ThreadFold']
CHRS = config["CHRS"]

# ENV_BIN = config["ENV_BIN"]
SCRIPT_DIR = config["SCRIPT_DIR"]
PIPE_DIR = config["PIPE_DIR"]
RULE_DIR = config["RULE_DIR"]
#SRC_DIR = config["SRC_DIR"]



include: RULE_DIR + "/NanoporeQC.rule.py"
include: RULE_DIR + "/MultiAssembly.rule.py"
include: RULE_DIR + "/AssemblyEvaluation.rule.py"
include: RULE_DIR + "/StructuralVariant.rule.py"
include: RULE_DIR + "/PangenomeGraph.rule.py"
include: RULE_DIR + "/NGSVariants.rule.py"
include: RULE_DIR + "/Methylation.rule.py"

rule all:
    input:
        #expand(IN_PATH + "/QualityControl/raw/{sample}_raw_stats.txt", sample=SAMPLES),
        #expand(IN_PATH + "/QualityControl/clean/{sample}_clean_stats.txt", sample=SAMPLES),
        #expand(IN_PATH + "/Assembly/NextDenovo/{sample}/03.ctg_graph/nd.asm.fasta", sample=SAMPLES)
        #expand(IN_PATH + "/QualityValue/NextDenovo/{sample}/{sample}.qv", sample=SAMPLES),
        #expand(IN_PATH + "/QualityValue/ONTPolish/{sample}/{sample}.qv", sample=SAMPLES),
        #expand(IN_PATH + "/QualityValue/NGSPolish/{sample}/{sample}.qv", sample=SAMPLES),
        #expand(IN_PATH + "/QualityControl/NGS/raw/{sample}_NGS.raw.stats.txt", sample=SAMPLES),
        #expand(IN_PATH + "/QualityControl/NGS/clean/{sample}_NGS.clean.stats.txt", sample=SAMPLES),
        #expand(IN_PATH + "/Centromere/SRF/{sample}_srf_aln_cluster.bed", sample=SAMPLES),
        #expand(IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_nextPolish3.fasta", sample=SAMPLES),
        #expand(IN_PATH + "/Assembly/Polish/{sample}/{sample}_contig_rank.txt", sample=SAMPLES),
        #IN_PATH + "/Assembly/Polish/Samples_contigs_rank.xls",
        #expand(IN_PATH + "/Assembly/RagTag/{sample}/ragtag.scaffold.stats.xls", sample=SAMPLES),
        #IN_PATH + "/Assembly/RagTag/All_sample_ragtag.scaffold.stats.xls",
        #IN_PATH + "/Assembly/Polish/Samples_contigs_rank.pdf",
        #expand(IN_PATH + "/SVCall/SV/{sample}.sniffles.vcf", sample=SAMPLES),
        #expand(IN_PATH + "/SVCall/SV/{sample}.NanoSV.vcf", sample=SAMPLES),
        #expand(IN_PATH + "/SVCall/SV/{sample}.cutesv.vcf", sample=SAMPLES),
        #expand(IN_PATH + "/Cactus/Assembly/{sample}.assembly.fasta", sample=SAMPLES),
        #IN_PATH + "/Cactus/minigraph/Assembly.minigraph.gfa",
        #expand(IN_PATH + "/QualityControl/NGS/raw/{sample}_NGS.raw.stats.txt", sample=SAMPLES),
        #expand(IN_PATH + "/QualityControl/NGS/clean/{sample}_NGS.clean.stats.txt", sample=SAMPLES),
        #expand(IN_PATH + '/NGS/mapping/{sample}/{sample}.sorted.bam', sample=SAMPLES),
        #expand(IN_PATH + '/NGS/mapping/{sample}/{sample}.sorted.bam', sample=SAMPLES),
        #expand(IN_PATH + '/NGS/VCF/{sample}/{sample}_{Chr}.g.vcf.gz', sample=SAMPLES, Chr=CHRS),
        expand(IN_PATH + "/SVCall/mapping/{sample}_ONT.bam", sample=SAMPLES),
        expand(IN_PATH + "/SVCall/Sniffles2/snf/{sample}.snf", sample=SAMPLES),
        expand(IN_PATH + '/SVCall/mapping/{sample}_ONT_bam_stats.xls', sample=SAMPLES),
        #expand(IN_PATH + "/TR/straglr/{sample}/{chr}_straglr_TR.tsv", sample=SAMPLES, chr=CHRS),
        #expand(IN_PATH + "/TR/straglr/{sample}_all_chr_straglr_TR.txt", sample=SAMPLES),
        #expand(IN_PATH + "/TR/straglr/{sample}_all_chr_straglr_TR.bed", sample=SAMPLES),



