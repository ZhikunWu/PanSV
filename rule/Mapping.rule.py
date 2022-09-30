

###################################################### mapping reads to genome ########################################################
rule BWAIndex:
    input:
        ref = config['REF'],
    output:
        bwt = config['REF'] + '.bwt',
        # fai = config['REF'] + '.fai',
    run:
        if not os.path.exists(output.bwt):
            cmd = 'bwa index %s' % input.ref
            os.system(cmd)
        # if not os.path.exists(output.fai):
        #     cmd2 = 'samtools faidx %s' % input.ref
        #     os.system(cmd2)


rule BWAAlign:
    input:
        R1 = IN_PATH + '/clean/fastp/{sample}.clean.R1.fq.gz',
        R2 = IN_PATH + '/clean/fastp/{sample}.clean.R2.fq.gz',
        bwt = config['REF'] + '.bwt',
    output:
        sam = temp(IN_PATH + '/mapping/{sample}/{sample}.sam'),
    threads:
        THREADS
    params:
        REF = config['REF'],
        # rg = "@RG\tID:01\tPL:ILLUMINA\tPU:{sample}\tSM:{sample}",
    log:
        IN_PATH + "/log/{sample}.bwa_align.log",
    run:
        shell('bwa mem -M -t {threads} {params.REF} {input.R1} {input.R2} > {output.sam} 2>{log}')




rule SAM2BAM:
    input:
        sam = IN_PATH + '/mapping/{sample}/{sample}.sam',
    output:
        bam = IN_PATH + '/mapping/{sample}/{sample}_align.bam',
    threads:
        THREADS 
    log:
        IN_PATH + "/log/{sample}.bam_sort.log",
    run:
        # shell("sambamba view -t {threads} -S -f  bam {input.sam} | sambamba sort -o {output.bam}  > {log} 2>&1")
        shell("samtools view -b -S {input.sam} | samtools sort -@ {threads} -o {output.bam} -m 5G  -  > {log} 2>&1")
        shell("samtools index {output.bam}")
###############################################################



############################ duplicate ###############################
rule MarkDuplicates:
    input:
        bam = IN_PATH + '/mapping/{sample}/{sample}_align.bam',
    output:
        bam = temp(IN_PATH + '/mapping/{sample}/{sample}_deduplicated.bam'),
        # bam = IN_PATH + '/mapping/{sample}/{sample}_deduplicated.bam',
        metrics = IN_PATH + '/mapping/{sample}/{sample}_dup_metrics.txt',
    params:
        picard = config["picard"],
        Memory = config["Memory"],
        # gatk4 = config["gatk4"],
    threads:
        THREADS
    log:
        IN_PATH + "/log/{sample}.markDuplicate.log",
    benchmark:
        IN_PATH + "/benchmark/MarkDuplicates_{sample}.txt"
    run:
        ### -Xmx30g
        shell('java {params.Memory} -jar {params.picard}  MarkDuplicates INPUT={input.bam} OUTPUT={output.bam} METRICS_FILE={output.metrics} REMOVE_DUPLICATES=true ASSUME_SORTED=true >{log} 2>&1')
        # shell("{params.gakt4} MarkDuplicates I={input.bam} O={output.bam} METRICS_FILE={output.metrics} REMOVE_DUPLICATES=true ASSUME_SORTED=true >{log} 2>&1")


rule AddGroup:
    input:
        bam = IN_PATH + '/mapping/{sample}/{sample}_deduplicated.bam',
    output:
        bam = temp(IN_PATH + '/mapping/{sample}/{sample}.sorted.bam'),
        bai = IN_PATH + '/mapping/{sample}/{sample}.sorted.bam.bai',
    log:
        IN_PATH + "/log/{sample}.addGroup.log",
    params:
        picard = config["picard"],
        Memory = config["Memory"],
        # gatk4 = config["gatk4"],
    threads:
        THREADS
    benchmark:
        IN_PATH + "/benchmark/AddGroup_{sample}.txt"
    run:

        # shell('java {params.Memory} -jar {params.picard}  AddOrReplaceReadGroups -INPUT {input.bam} -OUTPUT {output.bam} -SORT_ORDER coordinate -RGID {wildcards.sample} -RGPL Illumina -RGLB {wildcards.sample} -RGPU {wildcards.sample} -RGSM {wildcards.sample} >{log} 2>&1')
        shell('java {params.Memory} -jar {params.picard}  AddOrReplaceReadGroups INPUT={input.bam} OUTPUT={output.bam} SORT_ORDER=coordinate RGID={wildcards.sample} RGPL=Illumina RGLB={wildcards.sample} RGPU={wildcards.sample} RGSM={wildcards.sample} >{log} 2>&1')
        # shell("{params.gatk4} AddOrReplaceReadGroups I={input.bam} O={output.bam} RGID={wildcards.sample} RGPL=Illumina RGLB={wildcards.sample} RGPU={wildcards.sample} RGSM={wildcards.sample} SORT_ORDER=coordinate > {log} 2>&1")
        # shell("sambamba index --nthreads {threads} {output.bam} >{log} 2>&1")
        shell("samtools index -@ {threads} {output.bam} >{log} 2>&1")



##############################################################################################################################






######################################### mapping statistics #########################################

rule BAMStats:
    input:
        bam = IN_PATH + '/mapping/{sample}/{sample}_align.bam',
    output:
        stat = IN_PATH + '/mapping/{sample}/{sample}_bam_stats.xls',
        # summary = IN_PATH + "/mapping/{sample}/{sample}_bam_summary.xls",
    # params:
    #     SamtoolsBamStats = SRC_DIR + "/SamtoolsBamStats.py",
    threads:
        THREADS
    log:
        IN_PATH + "/log/BAMStats_{sample}.log"
    run:
        shell("samtools stats --threads {threads} {input.bam} > {output.stat} 2>{log}")
        # shell("python {params.SamtoolsBamStats} --stat {output.stat} --out {output.summary} --sample {wildcards.sample} 2>>{log}")


rule BAMStatsMerge:
    input:
        summary = expand(IN_PATH + "/mapping/{sample}/{sample}_bam_summary.xls", sample=SAMPLES),
    output:
        sumAll = IN_PATH + "/mapping/All_bam_stats_summary.xls",
    run:
        bamSums = input.summary
        merge_multiple_sample_stats(bamSums, output.sumAll)




# rule DepthCoverageStats:
#     input:
#         # bam = rules.InDelAlign.output.bam,
#         bam = IN_PATH + '/mapping/{sample}/{sample}.sorted.bam',
#     output:
#         depth = temp(IN_PATH + '/mapping/{sample}/{sample}_realigned.depth.txt'),
#         temp = IN_PATH + '/mapping/{sample}/{sample}_depth_coverage_temp.xls',
#         out = IN_PATH + '/mapping/{sample}/{sample}_depth_coverage.xls',
#     threads:
#         THREADS
#     params:
#         DepthCovStats = SRC_DIR + '/DepthCovStats.py',
#         REF = config['REF'],
#     log:
#         IN_PATH + "/log/{sample}_DepthCoverageStats.log",
#     run:
#         shell('samtools depth {input.bam} > {output.depth}')
#         shell('python {params.DepthCovStats} --input {output.depth} --fasta {params.REF} --out {output.temp} >{log} 2>&1')
#         shell("grep -v '^G' {output.temp} |  grep -v '^K' | grep -v '^M' > {output.out}")



# rule DepthCoveragePlot:
#     input:
#         depth = rules.DepthCoverageStats.output.out,
#     output:
#         sortedChr = IN_PATH + '/mapping/{sample}/{sample}_depth_coverage_sortedChr.xls',
#         ### the name is fixed and ends with ".Depth.Coverage.pdf"
#         pdf = IN_PATH + '/mapping/{sample}/{sample}.Depth.Coverage.pdf',
#     threads:
#         THREADS
#     params:
#         SortChrDigitPart = SRC_DIR + '/SortChrDigitPart.py',
#         DepthCoverPlot = SRC_DIR + '/DepthCoverPlot.py',
#         odir = IN_PATH + '/mapping/{sample}',
#         maxDepth = config['maxDepth'],
#     log:
#         IN_PATH + "/log/{sample}.DepthCoveragePlot.log",    
#     run:
#         shell('python {params.SortChrDigitPart} --input {input.depth} --out {output.sortedChr} >{log} 2>&1')
#         shell('python {params.DepthCoverPlot} --sample {wildcards.sample} --data {input.depth} --odir {params.odir} --depth {params.maxDepth} >{log} 2>&1')



# rule ReadDepthWin:
#     input:
#         depth = rules.DepthCoverageStats.output.depth,
#     output:
#         depth = IN_PATH + '/mapping/{sample}/{sample}_region_depth.txt',
#     threads:
#         THREADS
#     params:
#         SortReadDepthWin = SRC_DIR + '/SortReadDepthWin.py',
#         window = config['sliding'],
#     log:
#         IN_PATH + "/log/{sample}.ReadDepthWin.log",  
#     run:
#         shell('python {params.SortReadDepthWin} --depth {input.depth} --out {output.depth} --window {params.window} >{log} 2>&1')


# rule ReadDepthWinPlot:
#     input:
#         depth = rules.ReadDepthWin.output.depth,
#     output:
#         pdf = IN_PATH + '/mapping/{sample}/{sample}_region_depth.pdf',
#     threads:
#         THREADS
#     params:
#         ReadDepthPlot = SCRIPT_DIR + '/ReadDepthPlot.R',
#         window_size = config['sliding'],
#         chromosomes = config['chromosomes'],
#         ylim = config['maxDepth'],
#         width = config['width'],
#         height = config['height'],
#     log:
#         IN_PATH + "/log/{sample}.ReadDepthWinPlot.log",
#     run:
#         # shell('source activate Rmeta && Rscript {params.ReadDepthPlot} --depth {input.depth} --window_size {params.window_size} --chromosomes {params.chromosomes} --ylim {params.ylim} --width {params.width} --height {params.height} --pdf {output.pdf} > {log} 2>&1')
#         cmd = 'source activate Rmeta && Rscript %s --depth %s --window_size %s --chromosomes %s --ylim %s --width %s --height %s --pdf %s > %s 2>&1' % (params.ReadDepthPlot, input.depth, params.window_size, params.chromosomes, params.ylim, params.width, params.height, output.pdf, log)
#         os.system(cmd)
############################################################################################################





# ################################### pileup ####################################
# rule mpileup:
#     input:
#         bam = IN_PATH + '/mapping/{sample}/{sample}.sorted.bam',
#     output:
#         mpileup = IN_PATH + '/mapping/{sample}/{sample}_bam.mpileup',
#     params:
#         REF = config["REF"],
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/{sample}_mpileup.log",
#     run:
#         ### -L 10000
#         shell("samtools mpileup -B -q 1 -d 5000  -f {params.REF} {input.bam} > {output.mpileup} 2>{log}")



# rule mpileupTrio:
#     input:
#         bam = expand(IN_PATH + '/mapping/{sample}/{sample}.sorted.bam', sample=SAMPLES),
#     output:
#         mpileup = IN_PATH + '/mapping/trio/{trio}_bam_trio.mpileup',
#     params:
#         REF = config["REF"],
#         metafile = config["metafile"],
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/{trio}_mpileup.log",
#     run:
#         ### -L 10000
#         parents = MetaInfor(params.metafile).group_parents(wildcards.trio)
#         proband = MetaInfor(params.metafile).group_proband(wildcards.trio)
#         PileParents = "  ".join(sorted(target_file(parents, input.bam)))
#         PileProband = target_file(proband, input.bam)[0]
#         shell("samtools mpileup -B -q 1 -d 5000  -f {params.REF} {PileParents} {PileProband} > {output.mpileup} 2>{log}")



# rule VarScanDenovo:
#     input:
#         mpileup = expand(IN_PATH + '/mapping/{sample}/{sample}_bam.mpileup', sample=SAMPLES),
#     output:
#         snp = IN_PATH + '/VarScan/denovo/{trio}/{trio}_denovo.snp.vcf',
#         indel = IN_PATH + '/VarScan/denovo/{trio}/{trio}_denovo.indel.vcf',
#     params:
#         metafile = config["metafile"],
#         VarScan = config["VarScan"],
#         outPrefix = IN_PATH + '/VarScan/denovo/{trio}/{trio}_denovo',
#     log:
#         IN_PATH + "/log/{trio}_VarScanDenovo.log",
#     run:
#         parents = MetaInfor(params.metafile).group_parents(wildcards.trio)
#         proband = MetaInfor(params.metafile).group_proband(wildcards.trio)
#         PileParents = "  ".join(sorted(target_file(parents, input.mpileup)))
#         PileProband = target_file(proband, input.mpileup)[0]
#         shell("java -Xmx20g -jar {params.VarScan} trio {PileParents} {PileProband} --output-name {params.outPrefix} --min-coverage 5 --min-var-freq 0.2 --p-value 0.05  --adj-var-freq 0.05 --adj-p-value 0.15 >{log} 2>&1")       


# #################################################################################

