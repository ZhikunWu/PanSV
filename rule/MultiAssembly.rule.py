"""
################################### ONT NextDenovo Assembly #############################

rule ONTFasta:
    input:
        fastq = IN_PATH + "/clean/{sample}_ONT.fastq.gz",
    output:
        fastq = IN_PATH + "/clean/{sample}_ONT.fastq",
        # fasta = IN_PATH + "/clean/{sample}_ONT.fasta",
    threads:
        THREADS
    run:
        shell("pigz -p {threads} -dc {input.fastq} > {output.fastq}")
        # shell("seqtk seq -a {input.fastq} > {output.fasta}")



rule nextDenovo:
    input:
        fastq = rules.ONTFasta.output.fastq,
    output:
        fofn = IN_PATH + "/Assembly/NextDenovo/{sample}/{sample}_WGS_LRS.fofn",
        assembly = IN_PATH + "/Assembly/NextDenovo/{sample}/03.ctg_graph/nd.asm.fasta",
    threads:
        THREADS * ThreadFold
    params:
        configfile = IN_PATH + "/config/NextDenovo_{sample}.config",
    log:
        IN_PATH + "/log/nextDenovo_{sample}.log"
    run:
        shell("realpath {input.fastq} > {output.fofn}")
        shell("nextDenovo {params.configfile} > {log} 2>&1")



rule nextStats:
    ### https://github.com/raymondkiu/sequence-stats
    input:
        assembly = IN_PATH + "/Assembly/NextDenovo/{sample}/03.ctg_graph/nd.asm.fasta",
    output:
        stat = IN_PATH + "/Assembly/NextDenovo/{sample}/03.ctg_graph/nd.asm_each_stat.txt",
    run:
        shell("sequence-stats -c  {input.assembly} > {output.stat}")


rule nextPlot:
    input:
        stat = IN_PATH + "/Assembly/NextDenovo/{sample}/03.ctg_graph/nd.asm_each_stat.txt",
    output:
        pdf = IN_PATH + "/Assembly/NextDenovo/{sample}/03.ctg_graph/nd.asm_each_stat.pdf",
    params:
        TreeMap = SCRIPT_DIR + "/TreeMap.R",
        width = 5,
        height = 4,
    log:
        IN_PATH + "/log/nextPlot_{sample}.log"
    run:
        # shell("Rscript {params.TreeMap} --input {input.stat} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")
        cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --width %s --height %s > %s 2>&1" % (params.TreeMap, input.stat, output.pdf, params.width, params.height, log)
        print(cmd)
        os.system(cmd)

#####################################################################









# ####################################### racon polish using long read ###########################
rule ontAlign:
    input:
        fastq = IN_PATH + "/clean/{sample}_ONT.fastq.gz",
        contig = rules.nextDenovo.output.assembly,
    output:
        sam = temp(IN_PATH + "/Assembly/Polish/{sample}/Racon/{sample}_ONT.sam"),
    threads:
        THREADS * ThreadFold
    log:
        IN_PATH + "/log/mm2_ont_{sample}.log"
    run:
        shell("minimap2 --MD -a -x map-ont -secondary=no -t {threads} {input.contig} {input.fastq} > {output.sam} 2>{log}")




rule racon:
    input:
        sam = rules.ontAlign.output.sam,
        contig = rules.nextDenovo.output.assembly,
        fastq = IN_PATH + "/clean/{sample}_ONT.fastq.gz",
    output:
        contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT_polish_racon.fasta",
    threads:
        THREADS * ThreadFold
    log:
        IN_PATH + "/log/racon_{sample}.log"
    run:
        shell("racon --threads {threads} {input.fastq} {input.sam} {input.contig} > {output.contig} 2>{log}")





rule ontAlign2:
    input:
        fastq = IN_PATH + "/clean/{sample}_ONT.fastq.gz",
        contig = rules.racon.output.contig,
    output:
        sam = temp(IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT2.sam"),
    threads:
        THREADS * ThreadFold
    log:
        IN_PATH + "/log/ontAlign2_{sample}.log"
    run:
        shell("minimap2 --MD -a -x map-ont -secondary=no  -t {threads} {input.contig} {input.fastq} > {output.sam} 2>{log}")


rule racon2:
    input:
        sam = rules.ontAlign2.output.sam,
        contig = rules.racon.output.contig,
        fastq = IN_PATH + "/clean/{sample}_ONT.fastq.gz",
    output:
        contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT_polish_racon2.fasta",
    threads:
        THREADS * ThreadFold
    log:
        IN_PATH + "/log/racon2_{sample}.log"
    run:
        shell("racon --threads {threads} {input.fastq} {input.sam} {input.contig} > {output.contig} 2>{log}")




rule ontAlign3:
    input:
        fastq = IN_PATH + "/clean/{sample}_ONT.fastq.gz",
        contig = rules.racon2.output.contig,
    output:
        sam = temp(IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT3.sam"),
    threads:
        THREADS * ThreadFold
    log:
        IN_PATH + "/log/ontAlign3_{sample}.log"
    run:
        shell("minimap2 --MD -a -x map-ont -secondary=no  -t {threads} {input.contig} {input.fastq} > {output.sam} 2>{log}")



rule racon3:
    input:
        sam = rules.ontAlign3.output.sam,
        contig = rules.racon2.output.contig,
        fastq = IN_PATH + "/clean/{sample}_ONT.fastq.gz",
    output:
        contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT_polish_racon3.fasta",
    threads:
        THREADS * ThreadFold
    log:
        IN_PATH + "/log/racon3_{sample}.log"
    run:
        shell("racon --threads {threads} {input.fastq} {input.sam} {input.contig} > {output.contig} 2>{log}")
#################################################
"""


########## Pilon using short reads (three times)######
rule NGSAlign:
    input:
        contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT_racon3.fasta",
        R1 = rules.fastp.output.R1,
        R2 = rules.fastp.output.R2,
    output:
        bam = temp(IN_PATH + "/Assembly/Polish/{sample}/Pilon/{sample}_SRS_sorted.bam"),
    threads:
        THREADS
    log:
        IN_PATH + "/log/NGSAlign_{sample}.log"
    run:
        shell("bwa index {input.contig}")
        shell("bwa mem -t {threads} {input.contig} {input.R1} {input.R2} | samtools view -@ {threads} -F 0x4 -b - | samtools sort - -m 5g -@ {threads} -o {output.bam} > {log} 2>&1")
        shell("samtools index -@ {threads} {output.bam}")


rule pilon:
    input:
        contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT_racon3.fasta",
        bam = rules.NGSAlign.output.bam,
    output:
        contig = IN_PATH + "/Assembly/Polish/{sample}/Pilon/{sample}.SRS_pilon.fasta",
    threads:
        THREADS
    params:
        pilon = config["pilon"],
        outdir = IN_PATH + "/Assembly/Polish/{sample}/Pilon",
        memory = "-Xmx100G",
        mindepth = 10,
    log:
        IN_PATH + "/log/pilon_{sample}.log"
    run:
        shell("java  {params.memory}  -jar {params.pilon}  --genome {input.contig}  --bam {input.bam}  --outdir {params.outdir} --output {wildcards.sample}.SRS_pilon --threads {threads} --mindepth {params.mindepth} >{log} 2>&1")





rule NGSAlign2:
    input:
        contig = rules.pilon.output.contig,
        R1 = rules.fastp.output.R1,
        R2 = rules.fastp.output.R2,
    output:
        bam = temp(IN_PATH + "/Assembly/Polish/{sample}/Pilon/{sample}_SRS_sorted2.bam"),
    threads:
        THREADS
    log:
        IN_PATH + "/log/NGSAlign2_{sample}.log"
    run:
        shell("bwa index {input.contig}")
        shell("bwa mem -t {threads} {input.contig} {input.R1} {input.R2} | samtools view -@ {threads} -F 0x4 -b - | samtools sort - -m 5g -@ {threads} -o {output.bam} > {log} 2>&1")
        shell("samtools index -@ {threads} {output.bam}")


rule pilon2:
    input:
        contig = rules.pilon.output.contig,
        bam = rules.NGSAlign2.output.bam,
    output:
        contig = IN_PATH + "/Assembly/Polish/{sample}/Pilon/{sample}.SRS_pilon2.fasta",
    threads:
        THREADS
    params:
        pilon = config["pilon"],
        outdir = IN_PATH + "/Assembly/Polish/{sample}/Pilon",
        memory = "-Xmx100G",
        mindepth = 10,
    log:
        IN_PATH + "/log/pilon2_{sample}.log"
    run:
        shell("java  {params.memory}  -jar {params.pilon}  --genome {input.contig}  --bam {input.bam}  --outdir {params.outdir} --output {wildcards.sample}.SRS_pilon2 --threads {threads} --mindepth {params.mindepth} >{log} 2>&1")





rule NGSAlign3:
    input:
        contig = rules.pilon2.output.contig,
        R1 = rules.fastp.output.R1,
        R2 = rules.fastp.output.R2,
    output:
        bam = temp(IN_PATH + "/Assembly/Polish/{sample}/Pilon/{sample}_SRS_sorted3.bam"),
    threads:
        THREADS
    log:
        IN_PATH + "/log/NGSAlign3_{sample}.log"
    run:
        shell("bwa index {input.contig}")
        shell("bwa mem -t {threads} {input.contig} {input.R1} {input.R2} | samtools view -@ {threads} -F 0x4 -b - | samtools sort - -m 5g -@ {threads} -o {output.bam} > {log} 2>&1")
        shell("samtools index -@ {threads} {output.bam}")


rule pilon3:
    input:
        contig = rules.pilon2.output.contig,
        bam = rules.NGSAlign3.output.bam,
    output:
        contig = IN_PATH + "/Assembly/Polish/{sample}/Pilon/{sample}.SRS_pilon3.fasta",
    threads:
        THREADS
    params:
        pilon = config["pilon"],
        outdir = IN_PATH + "/Assembly/Polish/{sample}/Pilon",
        memory = "-Xmx100G",
        mindepth = 10,
    log:
        IN_PATH + "/log/pilon3_{sample}.log"
    run:
        shell("java  {params.memory}  -jar {params.pilon}  --genome {input.contig}  --bam {input.bam}  --outdir {params.outdir} --output {wildcards.sample}.SRS_pilon3 --threads {threads} --mindepth {params.mindepth} >{log} 2>&1")

#########################################



########################## Referenced scaffold #########################
# rule RagTag:
#     input:
#         # contig = rules.pilon.output.contig,
#         contig = rules.racon3.output.contig,
#     output:
#         scaffold = IN_PATH + "/Assembly/RagTag/{sample}/ragtag.scaffold.fasta",
#     params:
#         RefGenome = config["RefGenome"],
#         outDir = IN_PATH + "/Assembly/RagTag/{sample}",
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/RagTag_{sample}.log",
#     run:
#         shell("ragtag.py scaffold {params.RefGenome} {input.contig}  -o {params.outDir}  -w -u -t {threads} --aligner minimap2 > {log} 2>&1")
##############################################################









################ SyRi ##############
rule alignGenome:
    input:
        scaffold = "/home/wuzhikun/database/genome/arabidopsis/At_genome_seven/{scaffold}.chr.all.v2.0.fasta",
    output:
        sam = temp(IN_PATH + "/SyRI/{scaffold}/{scaffold}.sam"),
        bam = temp(IN_PATH + "/SyRI/{scaffold}/{scaffold}_align.bam"),
        sort = IN_PATH + "/SyRI/{scaffold}/{scaffold}_sort.bam",
    params:
        RefGenome = "/home/wuzhikun/database/genome/arabidopsis/TAIR10_chr_all.fas",
    threads:
        THREADS
    log:
        IN_PATH + "/log/alignGenome_{scaffold}.log",
    run:
        shell("minimap2 -ax asm5 --eqx {input.scaffold} {params.RefGenome}  > {output.sam} 2>{log}")
        shell("samtools view -b {output.sam} > {output.bam}")
        shell("samtools sort -@ {threads}  -o {output.sort} {output.bam}")
        shell("samtools index {output.sort}")


rule SyRI:
    input:
        bam = IN_PATH + "/SyRI/{scaffold}/{scaffold}_sort.bam",
        scaffold = "/home/wuzhikun/database/genome/arabidopsis/At_genome_seven/{scaffold}.chr.all.v2.0.fasta",
    output:
        syri = IN_PATH + "/SyRI/{scaffold}/SyRI/{scaffold}_aln_out.txt",
    params:
        RefGenome = "/home/wuzhikun/database/genome/arabidopsis/TAIR10_chr_all.fas",
    threads:
        THREADS
    log:
        IN_PATH + "/log/SyRI_{scaffold}.log",
    run:
        shell("syri -c {input.bam} -r {params.RefGenome}  -q {input.scaffold} -k -F B > {output.syri} 2>{log}")

    

#####################################



