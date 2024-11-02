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
THREADS = config["THREADS"]
#REVISED = config["REVISED"]
TISSUES = config["TISSUES"]
POPS = [1,2,3,4,5,6,7,8,9,10]


# ENV_BIN = config["ENV_BIN"]
SCRIPT_DIR = config["SCRIPT_DIR"]
PIPE_DIR = config["PIPE_DIR"]
RULE_DIR = config["RULE_DIR"]
#SRC_DIR = config["SRC_DIR"]
SRC_DIR = SCRIPT_DIR
#GROUPS = config["GROUPS"]

include: RULE_DIR + "/NanoporeQC.rule.py"
#include: RULE_DIR + "/MultiAssembly.rule.py"
#include: RULE_DIR + "/AssemblyEvaluation.rule.py"
#include: RULE_DIR + "/StructuralVariant.rule.py"
#include: RULE_DIR + "/PangenomeGraph.rule.py"
#include: RULE_DIR + "/NGSVariants.rule.py"
#include: RULE_DIR + "/NGSVariantsJoint.rule.py"
#include: RULE_DIR + "/Methylation.rule.py"
#include: RULE_DIR + "/GenePredictionFengchan.rule.py"
include: RULE_DIR + "/GenePrediction2.rule.py"
include: RULE_DIR + "/ProteinAnno.rule.py"
#include: RULE_DIR + "/PopulationVariant.rule.py"
#include: RULE_DIR + "/Vigna_RNA.rule.py"
include: RULE_DIR + "/SVMerge.rule.py"
include: RULE_DIR + "/PanVigna.rule.py"

rule all:
    input:
        expand(IN_PATH + "/GenePrediction/evidencemodeler/{sample}/{sample}.EVM.gff3", sample=SAMPLES),
        #expand(IN_PATH + "/GenePrediction/Gene/{sample}/{sample}.pep.fasta", sample=SAMPLES),
        #IN_PATH + "/NGS/Population/Chr_joint_variant_GT.vcf",
        #IN_PATH + "/NGS/Population/Chr_joint_variant_GT.bed",
        #IN_PATH + "/NGS/Population/PCA/Samples_PCA.eigenval",
        #IN_PATH + "/NGS/Population/PCA/Samples_PCA.txt",
        #IN_PATH + "/NGS/Population/PCA/Sample_group_PCA12.pdf",
        #IN_PATH + "/SVCall/Giraffe3/Sample884.giraffe.merge.vcf",
        #IN_PATH + "/NGS/Population/Sample34_joint_variant.bed",
        #IN_PATH + "/NGS/Population/PCA_sample34/Samples34_PCA.txt",
        #IN_PATH + "/NGS/Population/PCA_sample34/Samples34_group_PCA12.pdf",
        #IN_PATH + "/NGS/Population/Chr_joint_variant_GT.LDfilt.bed",
        #IN_PATH + "/NGS/Population/Chr_joint_variant_maf05.bed",
        #IN_PATH + "/NGS/Joint/VCF/Chrom_Fengchan6_joint_variant.stat.txt",
        #IN_PATH + "/NGS/Population/Chr_joint_variant.stat.txt",
