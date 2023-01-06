####################### HISAT alignment #####################
rule HisatIndex:
    input:
        assembly = rules.RagTag.output.assembly,
    output:
        ht = IN_PATH + "/GenePrediction/RNA/HisatIndex/{sample}_hisat.1.ht2",
    params:
        outDir = IN_PATH + "/GenePrediction/RNA/HisatIndex",
        prefix = IN_PATH + "/GenePrediction/RNA/HisatIndex/{sample}_hisat",
    threads:
        THREADS
    log:
        IN_PATH + "/log/HisatIndex_{sample}.log", 
    run:
        if not os.path.exists(params.outDir):
            os.makedirs(params.outDir)
        shell("hisat2-build  {input.assembly}  {params.prefix} > {log} 2>&1")

################################################################


########################### align each #################
rule hisatEach:
    input:
        ht = rules.HisatIndex.output.ht,
        R1 = IN_PATH + '/clean/{sample}.RNA.R1.fq.gz',
        R2 = IN_PATH + '/clean/{sample}.RNA.R2.fq.gz',
    output:
        sam = temp(IN_PATH + "/GenePrediction/RNA/hisat/{sample}/ngs_RNA_hisat.sam"),
    params:
        prefix = rules.HisatIndex.params.prefix,
    threads:
        THREADS
    log:
        IN_PATH + "/log/hisatEach_{sample}.log", 
    run:
        shell("hisat2 --dta -q --sensitive --threads {threads} -x {params.prefix} -1 {input.R1} -2 {input.R2} -S {output.sam} > {log} 2>&1")


rule hisatBAM:
    input:
        sam = rules.hisatEach.output.sam,
    output:
        bam = IN_PATH + "/GenePrediction/RNA/hisat/{sample}/ngs_RNA_hisat.bam"
    threads:
        THREADS
    log:
        IN_PATH + "/log/hisatbam_{sample}.log", 
    run:
        shell("samtools view -@ {threads} -b {input.sam} | samtools sort -  -@ {threads} -o {output.bam} > {log} 2>&1")
        shell("samtools index {output.bam}")


rule stringtie:
    input:
        bam = rules.hisatBAM.output.bam,
    output:
        gtf = IN_PATH + "/GenePrediction/RNA/Stringtie/{sample}/{sample}_RNA_stringtie.gtf",
    threads:
        THREADS
    log:
        IN_PATH + "/log/stringtie_{sample}.log",
    run:
        shell("stringtie {input.bam} -p {threads} -o {output.gtf} > {log} 2>&1")

##############################################################################





##################### TE ################
rule EDTA:
    input:
        assembly = rules.RagTag.output.assembly,
    output:
        log = IN_PATH + "/log/EDTA_{sample}.log",
    params:
        EDTA = config["EDTA"],
        CDSSeq = config["CDSSeq"],
        outDir = IN_PATH + "/Repeats/EDTA/{sample}",
    threads:
        THREADS
    log:
        IN_PATH + "/log/EDTA_{sample}.log",
    run:
        shell("cd {params.outDir} && perl {params.EDTA} --step all  --genome  {input.assembly}  --species others --cds {params.CDSSeq} --overwrite 1 --sensitive 1 --anno 1 --evaluate 1 --threads {threads} > {log} 2>&1")

########################################