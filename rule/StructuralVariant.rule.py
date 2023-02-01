
########### quality value #######
rule merylCount:
    input:
        # R1 = '/home/wuzhikun/PublicData/PRJNA836573/SRR19215719_clean.{pair}.fastq',
        R1 = IN_PATH + '/clean/{sample}.DNA.R1.fastq.gz',
        R2 = IN_PATH + '/clean/{sample}.DNA.R2.fastq.gz',
    output:
        R1 = directory(IN_PATH + "/Evaluation/QV/{sample}/{sample}_R1"),
        R2 = directory(IN_PATH + "/Evaluation/QV/{sample}/{sample}_R2"),
    threads:
        THREADS
    log:
        IN_PATH + "/log/merylCount_{sample}.log",
    run:
        shell("meryl count k=19 output {output.R1} {input.R1} > {log} 2>&1")
        shell("meryl count k=19 output {output.R2} {input.R2} > {log} 2>&1")

rule merylSum:
    input:
        R1 = IN_PATH + "/Evaluation/QV/{sample}/{sample}_R1",
        R2 = IN_PATH + "/Evaluation/QV/{sample}/{sample}_R2",
    output:
        merge = directory(IN_PATH + "/Evaluation/QV/{sample}/{sample}_merge"),
    threads:
        THREADS
    log:
        IN_PATH + "/log/merylSum_{sample}.log",
    run:
        Directory = " ".join(sorted([input.R1, input.R2]))
        shell("meryl union-sum output {output.merge} {Directory} > {log} 2>&1")


rule QV:
    input:
        merge = IN_PATH + "/Evaluation/QV/{sample}/{sample}_merge",
        assembly = rules.RagTag.output.scaffold,
    output:
        qv = IN_PATH + "/Evaluation/QV/{sample}/{sample}.qv",
    params:
        # qv = IN_PATH + "/Evaluation/QV/{sample}/{sample}",
        outDir = IN_PATH + "/Evaluation/QV/{sample}",
    threads:
        THREADS
    log:
        IN_PATH + "/log/QV_{sample}.log",
    run:
        shell("cd {params.outDir} && merqury.sh {input.merge}  {input.assembly}  {wildcards.sample} > {log} 2>&1")


################################



############### Evaluation #############

rule ONTAlign:
    input:
        fastq = IN_PATH + "/clean/{sample}_ONT.fastq.gz",
        assembly = rules.RagTag.output.scaffold,
    output:
        sam = temp(IN_PATH + "/SVCall/mapping/{sample}_ONT.sam"),
        bam = IN_PATH + "/SVCall/mapping/{sample}_ONT.bam",
    threads:
        THREADS
    log:
        IN_PATH + "/log/ONTAlign_{sample}.log"
    run:
        shell("minimap2 --MD -a -x map-ont -t {threads} {input.assembly} {input.fastq} > {output.sam} 2>{log}")
        shell("samtools view -@ {threads} -b {output.sam} | samtools sort -  -@ {threads} -o {output.bam} > {log} 2>&1")
        shell("samtools index {output.bam}")



rule ONTCov:
    input:
        bam = rules.ONTAlign.output.bam,
    output:
        bw = IN_PATH + "/Evaluation/Coverage/{sample}/{sample}_ONT_WGS.bw",
    log:
        IN_PATH + "/log/ONTCov_{sample}.log"
    run:
        shell("bamCoverage -b {input.bam} -o {output.bw} > {log} 2>&1")


rule ONTDepth:
    input:
        bam = rules.ONTAlign.output.bam,
    output:
        depth = IN_PATH +  "/Evaluation/Coverage/{sample}/{sample}_ONT.per-base.bed.gz",
    threads:
        THREADS
    params:
        outPrefix = IN_PATH + '/Evaluation/Coverage/{sample}/{sample}_ONT',
    log:
        IN_PATH + "/log/ONTDepth_{sample}.log"
    run:
        # shell('mosdepth --threads {threads} --fast-mode --flag 256  {params.outPrefix} {input.bam} > {log} 2>&1')
        cmd = "source activate NanoSV && mosdepth --threads %s --fast-mode --flag 256  %s %s > %s 2>&1" % (threads, params.outPrefix, input.bam, log)
        print(cmd)
        os.system(cmd)

########################################################




######################### SV ###########################
rule cuteSV:
    input:
        bam = IN_PATH + "/SVCall/mapping/{sample}_ONT.bam",
    output:
        vcf = IN_PATH + "/SVCall/cuteSV/{sample}.cutesv.vcf",
    threads:
        THREADS 
    params:
        RefGenome = config["RefGenome"],
        tempDir = IN_PATH + "/SVCall/cuteSV/temp/{sample}",
    log:
        IN_PATH + "/log/cuteSV_{sample}.log"
    run:
	    # For ONT data:
		# --max_cluster_bias_INS		100
		# --diff_ratio_merging_INS	0.3
		# --max_cluster_bias_DEL	100
		# --diff_ratio_merging_DEL	0.3
        if not os.path.exists(params.tempDir):
            os.makedirs(params.tempDir)
        # shell("cuteSV  --threads {threads} --sample  {wildcards.sample} --report_readid --min_support 3  --min_size 50 --max_size -1 --genotype {input.bam} {params.RefGenome}  {output.vcf} {params.tempDir} > {log} 2>&1")
        cmd = "source activate nanovar && cuteSV --max_cluster_bias_INS 100   --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --threads %s --sample  %s --report_readid --min_support 5  --min_size 50 --min_siglength 50  --max_size -1 --genotype %s  %s %s %s > %s  2>&1" % (threads, wildcards.sample, input.bam, params.RefGenome, output.vcf, params.tempDir, log)
        print(cmd)
        os.system(cmd)


rule MergeSV:
    input:
        vcf = expand(IN_PATH + "/SVCall/cuteSV/{sample}.cutesv.vcf", sample=SAMPLES),
    output:
        fileList = IN_PATH + "/SVCall/Merge/samples_file_list_DEL.txt",
        vcf = IN_PATH + "/SVCall/Merge/Sample_SVs_jasmine.vcf",
    threads:
        THREADS * ThreadFold
    log:
        IN_PATH + "/log/MergeSV.log",
    run:
        Files = "\n".join(sorted(input.vcf))
        out_h = open(output.fileList, "w")
        out_h.write("%s\n" % Files)
        out_h.close()
        # shell("jasmine --output_genotypes --ignore_strand --keep_var_ids threads={threads} file_list={output.fileList} out_file={output.vcf} > {log} 2>&1")
        cmd = "source activate nanovar && jasmine --output_genotypes --ignore_strand --keep_var_ids threads=%s file_list=%s out_file=%s > %s 2>&1" % (threads, output.fileList, output.vcf, log)
        os.system(cmd)

################################
