rule NanoQCRaw:
    input:
        fastq = IN_PATH + "/raw/{sample}.fastq.gz",
    output:
        temp1 = temp(IN_PATH + "/raw/{sample}_sub.fastq"),
        qc = IN_PATH + "/QualityControl/raw/NanoQC/{sample}/nanoQC.html",
    threads:
        THREADS
    params:
        outdir = IN_PATH + "/QualityControl/raw/NanoQC/{sample}",
        sub_freq = config["sub_freq"],
    log:
        IN_PATH + "/log/nanoQC_raw_{sample}.log"
    run:
        shell("seqtk sample -s100 {input.fastq} {params.sub_freq}  > {output.temp1} 2>>{log}")
        shell("nanoQC --outdir {params.outdir} {output.temp1} 2>>{log}")



rule RawStats:
    input:
        fastq = IN_PATH + "/raw/{sample}.fastq.gz",
    output:
        stats = IN_PATH + "/QualityControl/raw/{sample}_raw_stats.txt",
    run:
        ### nanovar2
        # shell("source activate nanovar2 && nanoq -i {input.fastq} -s > {output.summary}")
        cmd = "source activate nanovar2 && nanoq -i %s -s > %s" % (input.fastq, output.stats)
        print(cmd)
        os.system(cmd)


rule SeqFilt:
    input:
        fastq = IN_PATH + "/raw/{sample}.fastq.gz",
    output:
        fastq = IN_PATH + "/clean/{sample}.fastq.gz",
    run:
        cmd = "source activate nanovar2 && nanoq -i %s -S 70 -E 20 -l 1000 -q 7 -O g -c 6 > %s" % (input.fastq, output.fastq)
        print(cmd)
        os.system(cmd)



rule CleanStats:
    input:
        fastq = IN_PATH + "/clean/{sample}.fastq.gz",
    output:
        stats = IN_PATH + "/QualityControl/clean/{sample}_clean_stats.txt",
    run:
        ### nanovar2
        # shell("source activate nanovar2 && nanoq -i {input.fastq} -s > {output.summary}")
        cmd = "source activate nanovar2 && nanoq -i %s -s > %s" % (input.fastq, output.stats)
        print(cmd)
        os.system(cmd)