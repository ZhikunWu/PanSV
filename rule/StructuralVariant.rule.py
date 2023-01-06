



############### SV #############

rule ONTAlign:
    input:
        fastq = IN_PATH + "/clean/{sample}_ONT.fastq.gz"
    output:
        sam = temp(IN_PATH + "/SVCall/mapping/{sample}_ONT.sam"),
        bam = IN_PATH + "/SVCall/mapping/{sample}_ONT.bam",
    threads:
        THREADS
    params:
        RefGenome = config["RefGenome"],
    log:
        IN_PATH + "/log/ONTAlign_{sample}.log"
    run:
        shell("minimap2 --MD -a -x map-ont -t {threads} {params.RefGenome} {input.fastq} > {output.sam} 2>{log}")
        shell("samtools view -@ {threads} -b {output.sam} | samtools sort -  -@ {threads} -o {output.bam} > {log} 2>&1")
        shell("samtools index {output.bam}")



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
        cmd = "source activate nanovar && cuteSV --max_cluster_bias_INS 100   --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --threads %s --sample  %s --report_readid --min_support 3  --min_size 50 --min_siglength 50  --max_size -1 --genotype %s  %s %s %s > %s  2>&1" % (threads, wildcards.sample, input.bam, params.RefGenome, output.vcf, params.tempDir, log)
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
