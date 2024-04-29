# #################### Nanopore #################
# rule NanoQCRaw:
#     input:
#         fastq = IN_PATH + "/raw/{sample}_ONT.fastq.gz",
#     output:
#         # temp1 = temp(IN_PATH + "/raw/{sample}_sub.fastq"),
#         qc = IN_PATH + "/QualityControl/NanoQC/{sample}/nanoQC.html",
#     threads:
#         THREADS
#     params:
#         outdir = IN_PATH + "/QualityControl/NanoQC/{sample}",
#         # sub_freq = config["sub_freq"],
#     log:
#         IN_PATH + "/log/nanoQCRaw_{sample}.log"
#     run:
#         # shell("seqtk sample -s100 {input.fastq} {params.sub_freq}  > {output.temp1} 2>>{log}")
#         # shell("nanoQC --outdir {params.outdir} {input.fastq} 2>{log}")
#         cmd = "nanoQC --outdir %s %s 2>%s" % (params.outdir, input.fastq, log)
#         print(cmd)
#         os.system(cmd)



# rule RawStats:
#     input:
#         fastq = IN_PATH + "/raw/{sample}_ONT.fastq.gz",
#     output:
#         stats = IN_PATH + "/QualityControl/raw/{sample}_raw_stats.txt",
#     run:
#         ### nanovar2
#         # shell("source activate nanovar2 && nanoq -i {input.fastq} -s > {output.summary}")
#         cmd = "source activate nanovar2 && nanoq -i %s -s > %s" % (input.fastq, output.stats)
#         print(cmd)
#         os.system(cmd)



rule SeqFilt:
    input:
        fastq = IN_PATH + "/raw/{sample}_ONT.fastq.gz",
    output:
        fastq = IN_PATH + "/clean/{sample}_ONT.fastq.gz",
    run:
        cmd = "source activate nanovar2 && nanoq -i %s -S 70 -E 60 -l 1000 -q 7 -O g -c 6 > %s" % (input.fastq, output.fastq)
        print(cmd)
        os.system(cmd)



# rule CleanStats:
#     input:
#         fastq = IN_PATH + "/clean/{sample}_ONT.fastq.gz",
#     output:
#         stats = IN_PATH + "/QualityControl/clean/{sample}_clean_stats.txt",
#     run:
#         ### nanovar2
#         # shell("source activate nanovar2 && nanoq -i {input.fastq} -s > {output.summary}")
#         cmd = "source activate nanovar2 && nanoq -i %s -s > %s" % (input.fastq, output.stats)
#         print(cmd)
#         os.system(cmd)



# def sequence_stats_merge(filelist, out_file):
#     out_h = open(out_file, "w")
#     out_h.write("Sample\tReads\tBases\tN50_length\tLongest\tShortest\tMean_length\tMedian_length\t Mean_quality\tMedian_quality\n")
#     files = filelist.split(",")
#     for f in files:
#         name = f.strip().split("/")[-1].split("_")[0]
#         in_h = open(f, "r")
#         record = in_h.readline().strip()
#         out_h.write("%s\t%s\n" % (name, record))
#         in_h.close()
#     out_h.close()


# rule StatsMerge:
#     input:
#         raw = expand(IN_PATH + "/QualityControl/raw/{sample}_raw_stats.txt", sample=SAMPLES),
#         clean = expand(IN_PATH + "/QualityControl/clean/{sample}_clean_stats.txt", sample=SAMPLES),
#     output:
#         raw = IN_PATH + "/QualityControl/raw/All_samples_raw_stats.xls",
#         clean = IN_PATH + "/QualityControl/clean/All_samples_clean_stats.xls",
#     run:
#         RAW = ",".join(input.raw)
#         CLEAN = ",".join(input.clean)
#         sequence_stats_merge(RAW, output.raw)
#         sequence_stats_merge(CLEAN, output.clean)
# ###############################################################



# #################### NGS ################

rule fastp:
    input:
        R1 = IN_PATH + "/raw/{sample}_NGS.R1.fastq.gz",
        R2 = IN_PATH + "/raw/{sample}_NGS.R2.fastq.gz",
        # R1 = IN_PATH + "/raw/{sample}_R1.fq.gz",
        # R2 = IN_PATH + "/raw/{sample}_R2.fq.gz",
    output:
        R1 = IN_PATH + "/clean/{sample}_NGS.R1.fastq.gz",
        R2 = IN_PATH + "/clean/{sample}_NGS.R2.fastq.gz",
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


rule RawNGSStats:
    input:
        R1 = IN_PATH + "/raw/{sample}_NGS.R1.fastq.gz",
        R2 = IN_PATH + "/raw/{sample}_NGS.R2.fastq.gz",
        # R1 = IN_PATH + "/raw/{sample}_R1.fq.gz",
        # R2 = IN_PATH + "/raw/{sample}_R2.fq.gz",
    output:
        stats = IN_PATH + "/QualityControl/NGS/raw/{sample}_NGS.raw.stats.txt",
    threads:
        THREADS
    run:
        shell("seqkit stats -aT {input.R1} {input.R2} > {output.stats}")



rule MergeStat:
    input:
        stats = expand(IN_PATH + "/QualityControl/NGS/raw/{sample}_NGS.raw.stats.txt", sample=SAMPLES),
    output:
        stats = IN_PATH + "/QualityControl/NGS/Stats/S_Samples_NGS_all.raw.stats.xls",
    run:
        Stats = sorted(input.stats)
        print(len(Stats))
        for i in range(len(Stats)):
            s = Stats[i]
            if i == 0:
                cmd = "cat %s > %s" % (s, output.stats)
            else:
                cmd = "sed '1d' %s >> %s" % (s, output.stats)
            os.system(cmd)
            



rule NGSStats:
    input:
        R1 = IN_PATH + "/clean/{sample}_NGS.R1.fastq.gz",
        R2 = IN_PATH + "/clean/{sample}_NGS.R2.fastq.gz",
    output:
        stats = IN_PATH + "/QualityControl/NGS/clean/{sample}_NGS.clean.stats.txt",
    threads:
        THREADS
    run:
        shell("seqkit stats -aT {input.R1} {input.R2} > {output.stats}")


rule RNAStats:
    input:
        R1 = IN_PATH + "/clean/{sample}.RNA.R1.fq.gz",
        R2 = IN_PATH + "/clean/{sample}.RNA.R2.fq.gz",
    output:
        stats = IN_PATH + "/QualityControl/NGS/clean/{sample}_RNA.clean.stats.txt",
    threads:
        THREADS
    run:
        shell("seqkit stats -aT {input.R1} {input.R2} > {output.stats}")





rule MergeStatClean:
    input:
        stats = expand(IN_PATH + "/QualityControl/NGS/clean/{sample}_NGS.clean.stats.txt", sample=SAMPLES),
    output:
        stats = IN_PATH + "/QualityControl/NGS/Stats/S_Samples_NGS_all.clean.stats.xls",
    run:
        Stats = sorted(input.stats)
        for i in range(len(Stats)):
            s = Stats[i]
            if i == 0:
                cmd = "cat %s > %s" % (s, output.stats)
            else:
                cmd = "sed '1d' %s >> %s" % (s, output.stats)
            os.system(cmd)



rule ReadMd5:
    input:
        R1 = IN_PATH + "/raw/{sample}_NGS.R1.fastq.gz",
        R2 = IN_PATH + "/raw/{sample}_NGS.R2.fastq.gz",
    output:
        md5 = IN_PATH + "/QualityControl/md5/{sample}_raw_read.md5",
    run:
        shell("md5sum {input.R1} {input.R2} > {output.md5}")

#############################################################



