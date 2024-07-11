################## white & green ################
rule fastpRNA:
    input:
        R1 = IN_PATH + "/raw/{sample}_RNA.R1.fastq.gz",
        R2 = IN_PATH + "/raw/{sample}_RNA.R2.fastq.gz",
    output:
        R1 = IN_PATH + '/clean/{sample}.RNA.R1.fq.gz',
        R2 = IN_PATH + '/clean/{sample}.RNA.R2.fq.gz',
    threads:
        THREADS
    params:
        length_required = config["length_required"],
        qualified_quality_phred = config["qualified_quality_phred"],
        unqualified_percent_limit = config["unqualified_percent_limit"],
        cut_window_size = config["cut_window_size"],
        cut_mean_quality = config["cut_mean_quality"],
    log:
        IN_PATH + "/log/trim/{sample}.log", 
    run:
        shell("fastp --in1 {input.R1} --in2 {input.R2} --out1 {output.R1} --out2 {output.R2} --thread {threads} --compression 2 --length_required {params.length_required} --qualified_quality_phred {params.qualified_quality_phred}  --unqualified_percent_limit {params.unqualified_percent_limit} --cut_front --cut_tail --cut_window_size {params.cut_window_size} --cut_mean_quality {params.cut_mean_quality} >{log} 2>&1")


rule STAR:
    input:
        R1 = IN_PATH + '/clean/{sample}.RNA.R1.fq.gz',
        R2 = IN_PATH + '/clean/{sample}.RNA.R2.fq.gz',
        SA = IN_PATH + "/GenePrediction/STARIndex/Fengchan6/SA",
    output:
        bam = IN_PATH + '/RNA/NGS/mapping/{sample}Aligned.sortedByCoord.out.bam',
    params:
        IndexDir = IN_PATH + "/GenePrediction/STARIndex/Fengchan6",
        outDir = IN_PATH + '/RNA/NGS/mapping/{sample}',
    threads:
        THREADS
    log:
        IN_PATH + "/log/STAR/{sample}.log"
    run:
        if not os.path.exists(params.outDir):
            os.makedirs(params.outDir)
        shell(' STAR --runThreadN {threads} --genomeDir {params.IndexDir}  --readFilesIn {input.R1}  {input.R2}  --readFilesCommand gunzip -c --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {params.outDir}   > {log} 2>&1')
        shell("samtools index {output.bam}")



rule BAMStats2:
    input:
        bam = IN_PATH + '/RNA/NGS/mapping/{sample}Aligned.sortedByCoord.out.bam',
        bai = IN_PATH + '/RNA/NGS/mapping/{sample}Aligned.sortedByCoord.out.bam.bai',
    output:
        stat = IN_PATH + '/RNA/NGS/mapping/{sample}.bam_stats.xls',
    threads:
        THREADS
    log:
        IN_PATH + "/log/BAMStats2_{sample}.log"
    run:
        shell("samtools stats --threads {threads} {input.bam} > {output.stat} 2>{log}")

rule BAMStats3:
    input:
        bam = IN_PATH + "/RNA/NGS/mapping/{sample}Aligned.sortedByCoord.out.bam",
        bai = IN_PATH + '/RNA/NGS/mapping/{sample}Aligned.sortedByCoord.out.bam.bai',
    output:
        stat = IN_PATH + '/RNA/NGS/mapping/{sample}_bam_flag_stats.xls',
    threads:
        THREADS
    log:
        IN_PATH + "/log/BAMStats3_{sample}.log"
    run:
        shell("samtools flagstat --threads {threads} -O tsv {input.bam} > {output.stat} 2>{log}")

rule featureCount:
    input:
        IN_PATH + "/RNA/NGS/mapping/{sample}Aligned.sortedByCoord.out.bam",
    output:
        IN_PATH + '/RNA/NGS/Counts/{sample}_Reads_count.xls'
    threads:
        THREADS
    log:
        IN_PATH + '/log/{sample}_featurecount.log'
    params:
        feature_type = "exon", # "gene",  #config['feature_type'],
        feature_attribute = "Parent", # "ID", # config['feature_attribute'],
        GFF = IN_PATH + "/GenePrediction/Gene/Fengchan6/Fengchan6.gene.gff3",
    run:
        shell('featureCounts -T {threads} -d 50 -D 1000 -C -s 0 -t {params.feature_type} -g {params.feature_attribute} --primary -O -M --fraction -a {params.GFF} -o {output} {input} > {log} 2>&1' )


rule merge_featurecount:
    input:
        expand(IN_PATH + '/RNA/NGS/Counts/{sample}_Reads_count.xls', sample=SAMPLES)
    output:
        IN_PATH + '/RNA/NGS/Counts/Vigna_samples_reads_counts.xls',
    params:
        merge_featurecount = SCRIPT_DIR + '/merge_featurecount.py'
    log:
        IN_PATH + '/log/All_reads_counts.log',
    run:
        INPUT = ','.join(input)
        shell('python {params.merge_featurecount} -i {INPUT} -o {output}')


rule edgeR_DEG:
    input:
        IN_PATH + '/RNA/NGS/Counts/Vigna_samples_reads_counts.xls',
    output:
        down = IN_PATH + '/RNA/NGS/DEG/{group}/sig_Down_genes.xls',
        up = IN_PATH + '/RNA/NGS/DEG/{group}/sig_UP_genes.xls',
        diff = IN_PATH + '/RNA/NGS/DEG/{group}/RPKM_of_diff_genes.xls',
    params:
        Sample_pairs = config['GROUPS'],
        edgeR_DEG = '/home/wuzhikun/github/zkwu/kcmRNA/script/edgeR_DEG.R',
        edgeR_pvalue = 0.01,
        edgeR_logFC = 1,
        edgeR_fdr = "none", #fdr #none
        outDir = IN_PATH + '/RNA/NGS/DEG/{group}',
    log:
        IN_PATH + '/log/DEG/edgeR_DEG_{group}.log'
    run:
        if not os.path.exists(params.outDir):
            os.makedirs(params.outDir)
        pair = os.path.basename(params.outDir)
        print(pair)
        pairs = pair.split('___')
        sample1 = pairs[0]
        sample2 = pairs[1]
        if '__' in sample1 or '__' in sample2:
            CMD = 'Rscript  %s\t%s\t%s\t%s\t%s\t%f\t%f\t%s\t%s' % (params.edgeR_DEG, input, sample1, sample2, params.edgeR_fdr, params.edgeR_pvalue, params.edgeR_logFC, params.outDir, 'yes' )
            print(CMD)
            os.system(CMD)
        else:
            CMD = 'Rscript  %s\t%s\t%s\t%s\t%s\t%f\t%f\t%s\t%s' % (params.edgeR_DEG, input, sample1, sample2, params.edgeR_fdr, params.edgeR_pvalue, params.edgeR_logFC, params.outDir, 'no' )
            print(CMD)
            os.system(CMD)


##################################################
