
IN_PATH = "/home/wuzhikun/Project/Faba"
THREADS = 20

SAMPLES = ["Cajanus_cajan"]



rule all:
    input:
        expand(IN_PATH + "/HiC/yahs/{sample}/{sample}_scaffolds_final.fa", sample=SAMPLES),
        



rule chromapIndex:
    input:
        fasta = IN_PATH + "/Assembly/{sample}/03.ctg_graph/nd.asm.fasta",
    output:
        index = IN_PATH + "/Assembly/{sample}/03.ctg_graph/nd.asm.fasta.index",
    log:
        IN_PATH + "/log/chromapIndex_{sample}.log"
    run:
        shell("chromap -i -r {input.fasta}  -o {output.index} > {log} 2>&1")


rule chromap:
    input:
        index = IN_PATH + "/Assembly/{sample}/03.ctg_graph/nd.asm.fasta.index",
        fasta = IN_PATH + "/Assembly/{sample}/03.ctg_graph/nd.asm.fasta",
        R1 = IN_PATH + "/clean/{sample}.HiC.R1.fastq.gz",
        R2 = IN_PATH + "/clean/{sample}.HiC.R1.fastq.gz",
    output:
        sam = IN_PATH + "/HiC/yahs/{sample}/{sample}_HiC.sam",
    threads:
        THREADS
    log:
        IN_PATH + "/log/chromap_{sample}.log"
    run:
        shell("chromap --preset hic -x {input.index} -r {input.fasta} -1 {input.R1} -2 {input.R2} --remove-pcr-duplicates --SAM -o {output.sam} -t {threads} > {log} 2>&1")


rule sam2bam:
    input:
        sam = IN_PATH + "/HiC/yahs/{sample}/{sample}_HiC.sam",
    output:
        bam = IN_PATH + "/HiC/yahs/{sample}/{sample}_HiC.bam",
    threads:
        THREADS
    run:
        shell("samtools view -bh {input.sam} | samtools sort -@ {threads} -n > {output.bam}")


rule yahs:
    input:
        bam = IN_PATH + "/HiC/yahs/{sample}/{sample}_HiC.bam",
        fasta = IN_PATH + "/Assembly/{sample}/03.ctg_graph/nd.asm.fasta",
    output:
        fa = IN_PATH + "/HiC/yahs/{sample}/{sample}_scaffolds_final.fa",
        agp = IN_PATH + "/HiC/yahs/{sample}/{sample}_scaffolds_final.agp",
        b = IN_PATH + "/HiC/yahs/{sample}/{sample}.bin",
    params:
        outPrefix = IN_PATH + "/HiC/yahs/{sample}/{sample}",
    log:
        IN_PATH + "/log/yahs_{sample}.log"
    run:
        shell("samtools faidx {input.fasta}")
        shell("yahs {input.fasta} {input.bam} -o {params.outPrefix} > {log} 2>&1")


