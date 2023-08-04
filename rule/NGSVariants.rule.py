

###################################################### mapping reads to genome ########################################################
# rule BWAIndex:
#     input:
#         ref = config['REF'],
#     output:
#         bwt = config['REF'] + '.bwt',
#         # fai = config['REF'] + '.fai',
#     run:
#         if not os.path.exists(output.bwt):
#             cmd = 'bwa index %s' % input.ref
#             os.system(cmd)
#         # if not os.path.exists(output.fai):
#         #     cmd2 = 'samtools faidx %s' % input.ref
#         #     os.system(cmd2)


rule BWAAlign:
    input:
        R1 = IN_PATH + '/clean/{sample}_NGS.R1.fastq.gz',
        R2 = IN_PATH + '/clean/{sample}_NGS.R2.fastq.gz',
    output:
        sam = temp(IN_PATH + '/NGS/mapping/{sample}/{sample}.sam'),
    threads:
        THREADS
    params:
        RefGenome = config['RefGenome'],
        # rg = "@RG\tID:01\tPL:ILLUMINA\tPU:{sample}\tSM:{sample}",
    log:
        IN_PATH + "/log/{sample}.bwa_align.log",
    run:
        shell('bwa mem -M -t {threads} {params.RefGenome} {input.R1} {input.R2} > {output.sam} 2>{log}')


rule SAM2BAM:
    input:
        sam = IN_PATH + '/NGS/mapping/{sample}/{sample}.sam',
    output:
        tempbam = temp(IN_PATH + '/NGS/mapping/{sample}/{sample}_temp.bam'),
        bam = temp(IN_PATH + '/NGS/mapping/{sample}/{sample}_align.bam'),
    threads:
        THREADS
    log:
        IN_PATH + "/log/{sample}.bam_sort.log",
    run:
        shell('sambamba view --nthreads {threads} --sam-input --format bam  -o {output.tempbam} {input.sam} >{log} 2>&1' )
        shell('samtools sort --threads {threads} -o {output.bam} {output.tempbam} 2>>{log}')



rule MarkDuplicates:
    input:
        bam = IN_PATH + '/NGS/mapping/{sample}/{sample}_align.bam',
    output:
        bam = temp(IN_PATH + '/NGS/mapping/{sample}/{sample}_deduplicated.bam'),
        metrics = IN_PATH + '/NGS/mapping/{sample}/{sample}_dup_metrics.txt',
    params:
        picard = config["picard"],
    threads:
        THREADS
    log:
        IN_PATH + "/log/{sample}.markDuplicate.log",
    run:
        ### -Xmx30g
        shell('java -Xmx10g -jar {params.picard}  MarkDuplicates INPUT={input.bam} OUTPUT={output.bam} METRICS_FILE={output.metrics} REMOVE_DUPLICATES=true ASSUME_SORTED=true >{log} 2>&1')


rule AddGroup:
    input:
        bam = IN_PATH + '/NGS/mapping/{sample}/{sample}_deduplicated.bam',
    output:
        bam = IN_PATH + '/NGS/mapping/{sample}/{sample}.sorted.bam',
    log:
        IN_PATH + "/log/{sample}.addGroup.log",
    params:
        picard = config["picard"],
    threads:
        THREADS
    run:
        shell('java -Xmx10g -jar {params.picard}  AddOrReplaceReadGroups INPUT={input.bam} OUTPUT={output.bam} SORT_ORDER=coordinate RGID={wildcards.sample} RGPL=Illumina RGLB={wildcards.sample} RGPU={wildcards.sample} RGSM={wildcards.sample} >{log} 2>&1')
        shell("samtools index -@ {threads} {output.bam}")

#################################################


########################### GATK #################
rule HaplotypeCaller:
    input:
        bam = IN_PATH + '/NGS/mapping/{sample}/{sample}.sorted.bam',
    output:
        vcf = temp(IN_PATH + '/NGS/VCF/{sample}/{sample}_{Chr}.g.vcf'),
        # idx = IN_PATH + '/VCF/{sample}_{Chr}.gvcf.idx',
    params:
        GATK4 = config['GATK4'],
        REF = config['RefGenome'],
    threads:
        THREADS
    log:
        IN_PATH + "/log/{sample}_{Chr}.HaplotypeCaller.log",
    run:
        ### -L 20:10,000,000-10,200,000
        ### -L hg38.exome.regions.bed
        ### --dbsnp dbSNP.vcf
        ### --variant_index_type LINEAR
        ### --variant_index_parameter 128000
        ##### -nt / --num_threads controls the number of data threads sent to the processor (acting at the machine level)
        ##### -nct / --num_cpu_threads_per_data_thread controls the number of CPU threads allocated to each data thread (acting at the core level)
        ###  -variant_index_type LINEAR -variant_index_parameter 128000 
        shell('java -Xmx10g -jar {params.GATK4} HaplotypeCaller -R {params.REF} -L {wildcards.Chr} -I {input.bam} -O {output.vcf} -ERC GVCF --native-pair-hmm-threads {threads}  >{log} 2>&1') 
        ### -gt_mode DISCOVERY -stand_call_conf 30  -stand_emit_conf 30



rule bgzip:
    input:
        vcf = IN_PATH + '/NGS/VCF/{sample}/{sample}_{Chr}.g.vcf',
    output:
        vcf = IN_PATH + '/NGS/VCF/{sample}/{sample}_{Chr}.g.vcf.gz',
    run:
        shell("bgzip {input.vcf}")
##################################################



############# stats ####################
rule BAMStats:
    input:
        bam = IN_PATH + '/NGS/mapping/{sample}/{sample}.sorted.bam',
    output:
        stat = IN_PATH + '/NGS/mapping/{sample}/{sample}_bam_stats.xls',
    threads:
        THREADS
    log:
        IN_PATH + "/log/BAMStats_{sample}.log"
    run:
        shell("samtools stats --threads {threads} {input.bam} > {output.stat} 2>{log}")


########################################




# rule HaplotypeCallerMerge:
#     input:
#         vcf = expand(IN_PATH + '/VCF/{sample}_{Chr}.g.vcf', sample=SAMPLES, Chr=CHRS),
#         contig = IN_PATH + '/VCF/contig_{sample}.g.vcf',
#     output:
#         vcf = IN_PATH + '/VCF/{sample}.all.g.vcf',
#     threads:
#         THREADS
#     run:
#         vcfs = input.vcf
#         bases = ["_".join(v.split("/")[-1].split("_")[:-1]) for v in vcfs]
#         targetVCF = []
#         for v in range(len(bases)):
#             if wildcards.sample == bases[v]:
#                 targetVCF.append(vcfs[v])
#         for s in range(len(targetVCF)):
#             if s == 0:
#                 cmd = "cat %s > %s " % (targetVCF[s], output.vcf)
#                 os.system(cmd)
#             else:
#                 cmd = "grep -v '^#' %s >> %s " % (targetVCF[s], output.vcf)
#                 os.system(cmd)
#         cmd =  "grep -v '^#' %s >> %s " %  (input.contig, output.vcf)
#         os.system(cmd)


# rule HiFireplaceGATK2:
#     input:
#         RefGenome = IN_PATH + "/scaffold/MashMap/ont2/ont_merge_contig_replaced.fasta",
#         bam = IN_PATH + "/Evaluation/mapping/HiFireplace/Vigna_unguiculata_WGS_GATK.bam",
#     output:
#         vcf = IN_PATH + "/Evaluation/mapping/HiFireplace/Vigna_unguiculata_variant.vcf",
#     params:
#         GATK4 = config["GATK4"],
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/GATK_HiFireplace.log"  
#     run:
#         ### WGS
#         # shell("java  -Xmx50g -jar {params.GATK4} HaplotypeCaller  -I {input.bam}    -R {params.RefGenome}  -O {output.vcf} -L {wildcards.chr} -nct 1 -nt {threads} > {log} 2>&1")
#         ### java -jar /home/wuzhikun/anaconda3/envs/WGS/bin/picard.jar  CreateSequenceDictionary  R=/home/wuzhikun/Project/Vigna/HiC/ALLHiC/mapping/groups.asm.fasta  O=/home/wuzhikun/Project/Vigna/HiC/ALLHiC/mapping/groups.asm.dict
#         ### samtools faidx /home/wuzhikun/Project/Vigna/HiC/ALLHiC/mapping/groups.asm.fasta
#         cmd = "source activate WGS && java -Xmx50g -jar %s HaplotypeCaller  -I %s -R %s  -O %s  > %s 2>&1" % (params.GATK4, input.bam, input.RefGenome, output.vcf,  log)
#         print(cmd)
#         os.system(cmd)


# rule HiFireplaceSNVFilt2:
#     input:
#         vcf = IN_PATH + "/Evaluation/mapping/HiFireplace/Vigna_unguiculata_variant.vcf",
#     output:
#         snp = IN_PATH + "/Evaluation/mapping/HiFireplace/Vigna_unguiculata_variant_snv.vcf",
#         snp_filt = IN_PATH + "/Evaluation/mapping/HiFireplace/Vigna_unguiculata_variant_snv_filt.vcf",
#         indel = IN_PATH + "/Evaluation/mapping/HiFireplace/Vigna_unguiculata_variant_indel.vcf",
#         indel_filt = IN_PATH + "/Evaluation/mapping/HiFireplace/Vigna_unguiculata_variant_indel_filt.vcf",
#     params:
#         GATK3 = config['GATK3'],
#         REF = IN_PATH + "/scaffold/MashMap/ont2/ont_merge_contig_replaced.fasta",
#         Memory = config["Memory"],
#         snpFiltExpression = "'QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0'",
#         indelFiltExpression = "'QD < 2.0 ||  FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0'",
#     log:
#         IN_PATH + "/log/SNVFilt2_HiFireplace.log"  
#     run:
#         ### –filterExpression “QD <2.0 || ReadPosRankSum <−8.0 || FS >60.0 || MQ <40.0 || SOR >3.0 || MQRankSum <−10.0 || QUAL <30” –logging_level ERROR–missingValuesInExpressionsShouldEvaluateAsFailing
#         shell("java {params.Memory} -jar {params.GATK3} -T SelectVariants -R {params.REF} --variant {input.vcf} -selectType SNP -o {output.snp} >{log} 2>&1")
#         shell("java {params.Memory} -jar {params.GATK3} -T VariantFiltration -R {params.REF} --variant {output.snp} -o {output.snp_filt} --filterExpression {params.snpFiltExpression} --filterName 'SNPFilt' >>{log}")
#         shell("java {params.Memory} -jar {params.GATK3} -T SelectVariants -R {params.REF} --variant {input.vcf} -selectType INDEL -o {output.indel} >>{log}")
#         shell("java {params.Memory} -jar {params.GATK3} -T VariantFiltration -R {params.REF} --variant {output.indel} -o {output.indel_filt} --filterExpression {params.indelFiltExpression} --filterName 'InDelFilt' >>{log}")

#######################################################



##############################################################################################################################


############################### IndexDictionary ######################################
# rule IndexDictionary:
#     input:
#         ref = config['RefGenome'],
#     output:
#         fai = config['RefGenome']+ '.fai',
#         dic = config['RefGenome'].rstrip(".fasta") + ".dict",
#     params:
#         picard = config["picard"],
#     log:
#         IN_PATH + "/log/IndexDictionary.log",
#     run:
#         shell('samtools faidx {input.ref} >{log} 2>&1')
#         if input.ref.endswith('.fasta'):
#             base = input.ref.rstrip('.fasta')
#         else:
#             print('Please check the reference to make sure it ends with "fa", else there is something wrong with "GATK".')
#         base_dict = base + '.dict'
#         if not os.path.exists(base_dict):
#             cmd = 'java -Xmx30g -jar %s  CreateSequenceDictionary REFERENCE=%s > %s 2>&1' % (params.picard, input.ref, log)
#             os.system(cmd)
#########################################################################################  



# ######################################## Realign InDels and BaseRecalibrator #############################################



# rule RealignTarget:
#     input:
#         dic = rules.IndexDictionary.output.dic,
#         bai = IN_PATH + '/mapping/{sample}/{sample}.sorted.bam.bai',
#         bam = IN_PATH + '/mapping/{sample}/{sample}.sorted.bam',    
#     output:
#         #have one of the supported extensions (.bed, .list, .picard, .interval_list, or .intervals)
#         intervals = IN_PATH + '/mapping/{sample}/{sample}_realigned.intervals',
#     threads:
#         THREADS * threadFold
#     params:
#         GATK3 = config['GATK3'],
#         REF = config['REF'],
#         indels = config['indels'],
#         Memory = config["Memory"],
#     log:
#         IN_PATH + "/log/{sample}.realignTarget.log",
#     run:
#         shell('java {params.Memory} -jar {params.GATK3} -T RealignerTargetCreator -I {input.bam} -R {params.REF} -o {output.intervals}  -known {params.indels} -nt {threads} >{log} 2>&1') ### -known {params.indels} 



# rule InDelAlign:
#     input:
#         bam = IN_PATH + '/mapping/{sample}/{sample}.sorted.bam',
#         intervals = IN_PATH + '/mapping/{sample}/{sample}_realigned.intervals',
#     output:
#         bam = IN_PATH + '/mapping/{sample}/{sample}_realigned.bam', 
#     params:
#         GATK3 = config['GATK3'],
#         REF = config['REF'],
#         Memory = config["Memory"],
#         indels = config['indels'],
#     threads:
#         THREADS * threadFold
#     log:
#         IN_PATH + "/log/{sample}.InDelAlign.log",   
#     run:
#         shell('java {params.Memory} -jar {params.GATK3} -T IndelRealigner  -I {input.bam} -R {params.REF} -targetIntervals {input.intervals} -known {params.indels} -o {output.bam} >{log} 2>&1')
# ######################################################################################################






# ######################################### mapping statistics #########################################

# rule BAMStats:
#     input:
#         bam = IN_PATH + '/mapping/{sample}/{sample}.sorted.bam',
#     output:
#         stat = IN_PATH + '/mapping/{sample}/{sample}_bam_stats.xls',
#         summary = IN_PATH + "/mapping/{sample}/{sample}_bam_summary.xls",
#     params:
#         SamtoolsBamStats = SRC_DIR + "/SamtoolsBamStats.py",
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/BAMStats_{sample}.log"
#     run:
#         shell("samtools stats --threads {threads} {input.bam} > {output.stat} 2>{log}")
#         shell("python {params.SamtoolsBamStats} --stat {output.stat} --out {output.summary} --sample {wildcards.sample} 2>>{log}")


# rule BAMStatsMerge:
#     input:
#         summary = expand(IN_PATH + "/mapping/{sample}/{sample}_bam_summary.xls", sample=SAMPLES),
#     output:
#         sumAll = IN_PATH + "/mapping/All_bam_stats_summary.xls",
#     run:
#         bamSums = input.summary
#         merge_multiple_sample_stats(bamSums, output.sumAll)




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



