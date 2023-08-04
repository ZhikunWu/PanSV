



############### Evaluation #############

# rule ONTAlign:
#     input:
#         fastq = IN_PATH + "/clean/{sample}_ONT.fastq.gz",
#         assembly = rules.RagTag.output.scaffold,
#     output:
#         sam = temp(IN_PATH + "/SVCall/mapping/{sample}_ONT.sam"),
#         bam = IN_PATH + "/SVCall/mapping/{sample}_ONT.bam",
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/ONTAlign_{sample}.log"
#     run:
#         shell("minimap2 --MD -a -x map-ont -t {threads} {input.assembly} {input.fastq} > {output.sam} 2>{log}")
#         shell("samtools view -@ {threads} -b {output.sam} | samtools sort -  -@ {threads} -o {output.bam} > {log} 2>&1")
#         shell("samtools index {output.bam}")



# rule ONTCov:
#     input:
#         bam = rules.ONTAlign.output.bam,
#     output:
#         bw = IN_PATH + "/Evaluation/Coverage/{sample}/{sample}_ONT_WGS.bw",
#     log:
#         IN_PATH + "/log/ONTCov_{sample}.log"
#     run:
#         shell("bamCoverage -b {input.bam} -o {output.bw} > {log} 2>&1")


# rule ONTDepth:
#     input:
#         bam = rules.ONTAlign.output.bam,
#     output:
#         depth = IN_PATH +  "/Evaluation/Coverage/{sample}/{sample}_ONT.per-base.bed.gz",
#     threads:
#         THREADS
#     params:
#         outPrefix = IN_PATH + '/Evaluation/Coverage/{sample}/{sample}_ONT',
#     log:
#         IN_PATH + "/log/ONTDepth_{sample}.log"
#     run:
#         # shell('mosdepth --threads {threads} --fast-mode --flag 256  {params.outPrefix} {input.bam} > {log} 2>&1')
#         cmd = "source activate NanoSV && mosdepth --threads %s --fast-mode --flag 256  %s %s > %s 2>&1" % (threads, params.outPrefix, input.bam, log)
#         print(cmd)
#         os.system(cmd)

########################################################




######################### SV ###########################
rule ONTAlignRef:
    input:
        fastq = IN_PATH + "/clean/{sample}_ONT.fastq.gz",
        assembly = "/home/wuzhikun/Project/Vigna/Final/Final/Vigna_unguiculata_assembly.fasta",
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

rule BAMStats2:
    input:
        bam = IN_PATH + "/SVCall/mapping/{sample}_ONT.bam",
    output:
        stat = IN_PATH + '/SVCall/mapping/{sample}_ONT_bam_stats.xls',
    threads:
        THREADS
    log:
        IN_PATH + "/log/BAMStats2_{sample}.log"
    run:
        shell("samtools stats --threads {threads} {input.bam} > {output.stat} 2>{log}")


rule Sniffles1:
    input:
        bam = IN_PATH + "/SVCall/mapping/{sample}_ONT.bam",
    output:
        vcf = IN_PATH + "/SVCall/SV/{sample}.sniffles.vcf",
    threads:
        THREADS 
    log:
        IN_PATH + "/log/Sniffles1_{sample}.log"
    run:
        shell("sniffles --mapped_reads {input.bam} --vcf {output.vcf} --threads {threads}  --min_support 3 --min_length 50 --minmapping_qual 20 --num_reads_report -1 --min_seq_size 500  --genotype --report_BND --report-seq  >{log} 2>&1")


rule NanoSV1:
    input:
        bam = IN_PATH + "/SVCall/mapping/{sample}_ONT.bam",
    output:
        vcf = IN_PATH + "/SVCall/SV/{sample}.NanoSV.vcf",
    threads:
        THREADS 
    params:
        nanosvConfig = "/home/wuzhikun/anaconda3/envs/NanoSV/lib/python3.6/site-packages/nanosv/config.ini"
    log:
        IN_PATH + "/log/NanoSV_{sample}.log"
    run:
        ### depth_support = False
        cmd = "source activate NanoSV && NanoSV --threads %s -c %s -o %s %s > %s 2>&1" % (threads, params.nanosvConfig, output.vcf, input.bam, log)
        print(cmd)
        os.system(cmd)



rule cuteSV1:
    input:
        bam = IN_PATH + "/SVCall/mapping/{sample}_ONT.bam",
    output:
        vcf = IN_PATH + "/SVCall/SV/{sample}.cutesv.vcf",
    threads:
        THREADS 
    params:
        RefGenome = "/home/wuzhikun/Project/Vigna/Final/Final/Vigna_unguiculata_assembly.fasta",
        tempDir = IN_PATH + "/SVCall/SV/cuteSVTemp/{sample}",
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
        cmd = "source activate nanovar && cuteSV --max_cluster_bias_INS 100   --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --threads %s --sample  %s --report_readid --min_support 5  --min_size 50 --min_siglength 50  --max_size -1 --genotype %s  %s %s %s > %s  2>&1" % (threads, wildcards.sample, input.bam, params.RefGenome, output.vcf, params.tempDir, log)
        print(cmd)
        os.system(cmd)
##########################################################



##########################################################
rule snf:
    input:
        bam = IN_PATH + "/SVCall/mapping/{sample}_ONT.bam",
    output:
        snf = IN_PATH + "/SVCall/Sniffles2/snf/{sample}.snf",
        vcf = IN_PATH + "/SVCall/Sniffles2/snf/{sample}.vcf",
    threads:
        THREADS
    params:
        RefGenome = "/home/wuzhikun/Project/Vigna/Final/Final/Vigna_unguiculata_assembly.fasta",
        annotation = "/home/wuzhikun/Project/PanVigna/Vigna_SRF/Vigna_satellite_assembly_aln_3column.bed",
    log:
        IN_PATH + "/log/snf_{sample}.log"
    run:
        ### --tandem-repeats {params.annotation}
        shell("~/anaconda3/envs/sniffles2/bin/sniffles --threads {threads} --minsupport 5 --minsvlen 50 --sample-id {wildcards.sample} --reference {params.RefGenome} --tandem-repeats {params.annotation} --input {input.bam} --vcf {output.vcf} --snf {output.snf}  --long-ins-length 50000  > {log} 2>&1")


rule MergeVCF:
    input:
        snf = expand(IN_PATH + "/SVCall/Sniffles2/snf/{sample}.snf", sample=SAMPLES),
    output:
        vcf = IN_PATH + "/SVCall/Sniffles2/Samples_SV_merge.vcf",
    threads:
        THREADS
    params:
        RefGenome = "/home/wuzhikun/Project/Vigna/Final/Final/Vigna_unguiculata_assembly.fasta",
        annotation = "/home/wuzhikun/Project/PanVigna/Vigna_SRF/Vigna_satellite_assembly_aln_3column.bed",
    log:
        IN_PATH + "/log/MergeVCF.log"
    run:
        ### --tandem-repeats {params.annotation}
        Samples = " ".join(sorted(input.snf))
        shell("~/anaconda3/envs/sniffles2/bin/sniffles --input {Samples} --vcf {output.vcf} --threads {threads}  --reference {params.RefGenome} --tandem-repeats {params.annotation} --minsupport 5 --minsvlen 50  --long-ins-length 50000 > {log} 2>&1")




# rule MergeSV:
#     input:
#         vcf = expand(IN_PATH + "/SVCall/cuteSV/{sample}.cutesv.vcf", sample=SAMPLES),
#     output:
#         fileList = IN_PATH + "/SVCall/Merge/samples_file_list_DEL.txt",
#         vcf = IN_PATH + "/SVCall/Merge/Sample_SVs_jasmine.vcf",
#     threads:
#         THREADS * ThreadFold
#     log:
#         IN_PATH + "/log/MergeSV.log",
#     run:
#         Files = "\n".join(sorted(input.vcf))
#         out_h = open(output.fileList, "w")
#         out_h.write("%s\n" % Files)
#         out_h.close()
#         # shell("jasmine --output_genotypes --ignore_strand --keep_var_ids threads={threads} file_list={output.fileList} out_file={output.vcf} > {log} 2>&1")
#         cmd = "source activate nanovar && jasmine --output_genotypes --ignore_strand --keep_var_ids threads=%s file_list=%s out_file=%s > %s 2>&1" % (threads, output.fileList, output.vcf, log)
#         os.system(cmd)




rule straglr:
    input:
        bam = IN_PATH + "/SVCall/mapping/{sample}_ONT.bam",
    output:
        TR = IN_PATH + "/TR/straglr/{sample}/{chr}_straglr_TR.tsv",
        bed = IN_PATH + "/TR/straglr/{sample}/{chr}_straglr_TR.bed",
    threads:
        THREADS
    params:
        RefGenome = "/home/wuzhikun/Project/Vigna/Final/Final/Vigna_unguiculata_assembly.fasta",
        outdir = IN_PATH + "/TR/straglr/{sample}",
        tmpdir = IN_PATH + "/TR/Temp/{sample}",
        TR = IN_PATH + "/TR/straglr/{sample}/{chr}_straglr_TR",
    log:
        IN_PATH + "/log/straglr_{sample}_{chr}.log",
    run:
        ### --exclude hg38.exclude.bed
        ### source actiuvate HiC-Pro
        if not os.path.exists(params.outdir):
            os.makedirs(params.outdir)
        if not os.path.exists(params.tmpdir):
            os.makedirs(params.tmpdir)
        # shell("source activate HiC-Pro && cd {params.outdir} && straglr.py {input.bam}  {params.RefGenome} {params.TR}  --min_str_len 2 --max_str_len 100 --min_ins_size 50 --min_support 5  --chroms {wildcards.chr} --nprocs {threads} --tmpdir {params.tmpdir}  --genotype_in_size > {log} 2>&1")
        cmd = "source activate HiC-Pro && cd %s && straglr.py %s  %s %s  --min_str_len 2 --max_str_len 100 --min_ins_size 50 --min_support 5  --chroms %s --nprocs %s --tmpdir %s  --genotype_in_size > %s 2>&1" % (params.outdir, input.bam, params.RefGenome, params.TR, wildcards.chr, threads, params.tmpdir, log)
        print(cmd)
        os.system(cmd)



rule mergeRT:
    input:
        TR = expand(IN_PATH + "/TR/straglr/{sample}/{chr}_straglr_TR.tsv", sample=SAMPLES, chr=CHRS),
    output:
        TR = IN_PATH + "/TR/straglr/{sample}_all_chr_straglr_TR.txt",
    run:
        Target = []
        files = sorted(input.TR)
        for i in files:
            s = i.split("/")[-2]
            if s == wildcards.sample:
                Target.append(i)
        print(Target)
        for j in range(len(Target)):
            if j == 0:
                cmd = "cat %s > %s" % (Target[j], output.TR)
            else:
                cmd = "grep -v  '^#' %s >> %s" % (Target[j], output.TR)
            os.system(cmd)



rule mergeRT2:
    input:
        bed = expand(IN_PATH + "/TR/straglr/{sample}/{chr}_straglr_TR.bed", sample=SAMPLES, chr=CHRS),
    output:
        bed = IN_PATH + "/TR/straglr/{sample}_all_chr_straglr_TR.bed",
    run:
        Target = []
        files = sorted(input.bed)
        for i in files:
            s = i.split("/")[-2]
            if s == wildcards.sample:
                Target.append(i)
        print(Target)
        for j in range(len(Target)):
            if j == 0:
                cmd = "cat %s > %s" % (Target[j], output.bed)
            else:
                cmd = "grep -v  '^#' %s >> %s" % (Target[j], output.bed)
            os.system(cmd)

################################
