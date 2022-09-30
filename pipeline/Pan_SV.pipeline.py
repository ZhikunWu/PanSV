#!/usr/bin/env python
# coding: utf-8
import os
import sys
import yaml


### plot workflow
### snakemake -s /home/wuzhikun/github/TrioWGS/pipeline/TrioWGS.pipeline.py  --configfile /home/wuzhikun/github/TrioWGS/pipeline/TrioWGS.pipeline.yaml -j 48 --dag | dot -Tsvg > TrioWGS.svg



IN_PATH = config['project_path']
THREADS = config['THREADS']
threadFold = config['threadFold']
SAMPLES = config['SAMPLES']


# ENV_BIN = config["ENV_BIN"]
SCRIPT_DIR = config["SCRIPT_DIR"]
PIPE_DIR = config["PIPE_DIR"]
RULE_DIR = config["RULE_DIR"]
SRC_DIR = config["SRC_DIR"]


# include: RULE_DIR + '/BaseFunction.rule.py'
include: RULE_DIR + '/QualityControl.rule.py'
include: RULE_DIR + '/Mapping.rule.py'



rule all:
    input:
        ################################# Quality control ############################
        expand(IN_PATH + '/FastQC/clean/{sample}_stats.xls', sample=SAMPLES),
        ################################ mapping ####################
        expand(IN_PATH + '/mapping/{sample}/{sample}_align.bam', sample=SAMPLES),
        expand(IN_PATH + '/mapping/{sample}/{sample}_bam_stats.xls', sample=SAMPLES),

