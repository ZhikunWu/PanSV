
######################### DNA Quality control  ####################################################
rule fastp:
    input:
        R1 = IN_PATH + '/raw/{sample}.html_1.fastq.gz',
        R2 = IN_PATH + '/raw/{sample}.html_2.fastq.gz',
    output:
        R1 = IN_PATH + '/clean/fastp/{sample}.clean.R1.fq.gz',
        R2 = IN_PATH + '/clean/fastp/{sample}.clean.R2.fq.gz',
    threads:
        THREADS
    params:
        length_required = config["length_required"],
        qualified_quality_phred = config["qualified_quality_phred"],
        unqualified_percent_limit = config["unqualified_percent_limit"],
        cut_window_size = config["cut_window_size"],
        cut_mean_quality = config["cut_mean_quality"],
    log:
        IN_PATH + "/log/trim/{sample}_fastp.log", 
    run:
        ### adapter parameter: --adaptersequence (read1), --adapter_sequence_r2(read2)
        ### do not trim adapter: --disable_adapter_trimming
        shell("fastp --in1 {input.R1} --in2 {input.R2} --out1 {output.R1} --out2 {output.R2} --thread {threads} --compression 2 --length_required {params.length_required} --qualified_quality_phred {params.qualified_quality_phred}  --unqualified_percent_limit {params.unqualified_percent_limit} --cut_front --cut_tail --cut_window_size {params.cut_window_size} --cut_mean_quality {params.cut_mean_quality} >{log} 2>&1")



rule fastqc:
    input:
        clean_R1 = rules.fastp.output.R1,
        clean_R2 = rules.fastp.output.R2,
    output:
        clean_R1 = IN_PATH + '/FastQC/clean/{sample}/{sample}.clean.R1_fastqc.zip',  
        clean_R2 = IN_PATH + '/FastQC/clean/{sample}/{sample}.clean.R2_fastqc.zip',
    threads:
        THREADS
    log:
        clean = IN_PATH + "/log/fastqc/{sample}_clean.log",
    params:
        clean_dir = IN_PATH + '/FastQC/clean/{sample}',
    run:
        shell('fastqc --threads {threads} --extract -f fastq {input.clean_R1} {input.clean_R2} -o {params.clean_dir} > {log.clean} 2>&1 ')





rule QCStats:
    input:
        clean_R1 = rules.fastp.output.R1,
        clean_R2 = rules.fastp.output.R2,
    output:
        clean_R1 = IN_PATH + '/FastQC/clean/{sample}_R1_stats.xls',
        clean_R2 = IN_PATH + '/FastQC/clean/{sample}_R2_stats.xls',
        clean = IN_PATH + '/FastQC/clean/{sample}_stats.xls',
    threads:
        THREADS
    log:
        IN_PATH + "/log/{sample}_QCStats.log",
    params:
        FastxStatMerge = SRC_DIR + '/FastxStatMerge.py',
        TrimStats = SRC_DIR + "/TrimStats.py",
    run:
        shell("seqstats {input.clean_R1} > {output.clean_R1}")
        shell("seqstats {input.clean_R2} > {output.clean_R2}")
        shell("python {params.FastxStatMerge} --stat {output.clean_R1} --stat2 {output.clean_R2} --sample {wildcards.sample} --out {output.clean} 2>>{log}")




##############################################################################################


################## RNA quality control #############
rule fastpRNA:
    input:
        R1 = IN_PATH + '/raw/{sample}.RNA.R1.fastq.gz',
        R2 = IN_PATH + '/raw/{sample}.RNA.R2.fastq.gz',
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
        IN_PATH + "/log/trim/{sample}_fastp.log", 
    run:
        ### adapter parameter: --adaptersequence (read1), --adapter_sequence_r2(read2)
        ### do not trim adapter: --disable_adapter_trimming
        shell("fastp --in1 {input.R1} --in2 {input.R2} --out1 {output.R1} --out2 {output.R2} --thread {threads} --compression 3 --length_required {params.length_required} --qualified_quality_phred {params.qualified_quality_phred}  --unqualified_percent_limit {params.unqualified_percent_limit} --cut_front --cut_tail --cut_window_size {params.cut_window_size} --cut_mean_quality {params.cut_mean_quality} >{log} 2>&1")
###################################################




################## Long read quality control ################
rule NanoQCRaw:
    input:
        fastq = IN_PATH + "/raw/{sample}_ONT.fastq.gz",
    output:
        temp1 = temp(IN_PATH + "/raw/{sample}_sub.fastq"),
        qc = IN_PATH + "/QualityControl/raw/NanoQC/{sample}/nanoQC.html",
    threads:
        THREADS
    params:
        outdir = IN_PATH + "/QualityControl/raw/NanoQC/{sample}",
    log:
        IN_PATH + "/log/nanoQC_raw_{sample}.log"
    run:
        shell("seqtk sample -s100 {input.fastq} 0.01  > {output.temp1} 2>>{log}")
        shell("nanoQC --outdir {params.outdir} {output.temp1} 2>>{log}")


rule NanoFilt:
    input:
        fastq = IN_PATH + "/raw/{sample}_ONT.fastq.gz",
    output:
        fastq = IN_PATH + "/clean/{sample}_ONT.fastq.gz",
    threads:
        THREADS
    log:
        IN_PATH + "/log/NanoFilt_{sample}.log"
    params:
        readtype = config["readtype"],
        minQuality = config["minQuality"],
        minLength = config["minLength"],
        headcrop = config["headcrop"],
        tailcrop = config["tailcrop"],
    run:
        # ### Filt out reads of phage lambda
        # shell("gunzip -c {input.fastq} | NanoLyse --reference {params.phage_lambda} | NanoFilt --quality  {params.minQuality} --length {params.minLength} --headcrop {params.headcrop} --tailcrop {params.tailcrop} | gzip > {output.fastq} 2>{log}")
        # shell("gunzip -c {input.fastq} | NanoFilt --readtype {params.readtype} --quality  {params.minQuality} --length {params.minLength} --headcrop {params.headcrop} --tailcrop {params.tailcrop} | gzip > {output.fastq} 2>{log}")
        cmd = "source activate NanoSV &&  gunzip -c %s | NanoFilt --readtype %s --quality  %s --length %s --headcrop %s --tailcrop %s | gzip > %s 2>%s" % (input.fastq, params.readtype, params.minQuality, params.minLength, params.headcrop, params.tailcrop, output.fastq, log)
        print(cmd)
        os.system(cmd)
        # shell("NanoFilt --readtype {params.readtype} --quality  {params.minQuality} --length {params.minLength} --headcrop {params.headcrop} --tailcrop {params.tailcrop} < {input.fastq} | gzip > {output.fastq} 2>{log}")
################################################################