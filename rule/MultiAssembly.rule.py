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


#################### PolyPolish ################
rule bwaAlign:
    input:
        contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT_racon3.fasta",
        R1 = rules.fastp.output.R1,
        R2 = rules.fastp.output.R2,
    output:
        sam1 = temp(IN_PATH + "/Assembly/Polypolish/{sample}/{sample}_aln1.R1.sam"), 
        sam2 = temp(IN_PATH + "/Assembly/Polypolish/{sample}/{sample}_aln1.R2.sam"), 
    threads:
        THREADS
    log:
        IN_PATH + "/log/bwaAlign_{sample}.log"
    run:
        shell("bwa index {input.contig}")
        shell("bwa mem -t {threads} -a {input.contig} {input.R1} > {output.sam1} 2>{log}")
        shell("bwa mem -t {threads} -a {input.contig} {input.R2} > {output.sam2} 2>>{log}")


rule InsertFilter:
    input:
        sam1 = rules.bwaAlign.output.sam1,
        sam2 = rules.bwaAlign.output.sam2,
    output:
        sam1 = temp(IN_PATH + "/Assembly/Polypolish/{sample}/{sample}_filtered1.R1.sam"),
        sam2 = temp(IN_PATH + "/Assembly/Polypolish/{sample}/{sample}_filtered1.R2.sam"),
    threads:
        THREADS
    log:
        IN_PATH + "/log/InsertFilter_{sample}.log"
    run:
        shell("polypolish_insert_filter.py --in1 {input.sam1} --in2 {input.sam2} --out1 {output.sam1} --out2 {output.sam2} > {log} 2>&1")



rule polypolish:
    input:
        contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT_racon3.fasta",
        sam1 = rules.InsertFilter.output.sam1,
        sam2 = rules.InsertFilter.output.sam2,
    output:
        contig = IN_PATH + "/Assembly/Polypolish/{sample}/{sample}_polypolish1.fasta",
    threads:
        THREADS
    log:
        IN_PATH + "/log/polypolish_{sample}.log"
    run:
        shell("polypolish {input.contig} {input.sam1} {input.sam2} > {output.contig} >{log} 2>&1")
###############################################










##################### nextPolish #############
rule ngsAlign1:
    input:
        contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT_racon3.fasta",
        R1 = rules.fastp.output.R1,
        R2 = rules.fastp.output.R2,
    output:
        bam = temp(IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_sort1.bam"),
    threads:
        THREADS
    run:
        shell("samtools faidx {input.contig}")
        shell("bwa index {input.contig}")
        shell("bwa mem -t {threads} {input.contig} {input.R1} {input.R2} | samtools view --threads {threads} -F 0x4 -b - | samtools fixmate -m --threads {threads}  - -| samtools sort -m 2g --threads {threads} -|samtools markdup --threads {threads} -r - {output.bam}")
        shell("samtools index -@ {threads} {output.bam}")



rule nextPolish1:
    input:
        contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT_racon3.fasta",
        bam = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_sort1.bam",
    output:
        contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_nextPolish_temp.fasta",
    params:
        nextpolish1 = "/home/wuzhikun/software/NextPolish/lib/nextpolish1.py",
    threads:
        THREADS
    run:
        shell("python {params.nextpolish1} -g {input.contig} -t 1 -p {threads} -s {input.bam} > {output.contig}")





rule ngsAlign2:
    input:
        contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_nextPolish_temp.fasta",
        R1 = rules.fastp.output.R1,
        R2 = rules.fastp.output.R2,
    output:
        bam = temp(IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_sort2.bam"),
    threads:
        THREADS
    run:
        shell("samtools faidx {input.contig}")
        shell("bwa index {input.contig}")
        shell("bwa mem -t {threads} {input.contig} {input.R1} {input.R2} | samtools view --threads {threads} -F 0x4 -b - | samtools fixmate -m --threads {threads}  - -|samtools sort -m 2g --threads {threads} -|samtools markdup --threads {threads} -r - {output.bam}")
        shell("samtools index -@ {threads} {output.bam}")



rule nextPolish2:
    input:
        contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_nextPolish_temp.fasta",
        bam = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_sort2.bam",
    output:
        contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_nextPolish1.fasta",
    params:
        nextpolish1 = "/home/wuzhikun/software/NextPolish/lib/nextpolish1.py",
    threads:
        THREADS
    run:
        shell("python {params.nextpolish1} -g {input.contig} -t 2 -p {threads} -s {input.bam} > {output.contig}")


##############################################


# ######################### Referenced scaffold #########################
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
# #############################################################









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



