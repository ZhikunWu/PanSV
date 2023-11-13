
################################### ONT NextDenovo Assembly #############################

def replace_config(in_file, out_file, seq_file, outDir):
    in_h = open(in_file, "r")
    out_h = open(out_file, "w")
    for line in in_h:
        line = line.strip()
        if line.startswith("input_fofn"):
            out_h.write("input_fofn = %s\n" % seq_file)
        elif line.startswith("workdir"):
            out_h.write("workdir = %s\n" % outDir)
        else:
            out_h.write("%s\n" % line)
    in_h.close()
    out_h.close()


rule nextDenovoConfig:
    input:
        fastq = rules.SeqFilt.output.fastq,
        config = IN_PATH + "/config/NextDenovo.cfg",
    output:
        fofn = IN_PATH + "/Assembly/NextDenovo/fofn/{sample}_nextDenovo.fofn",
        config = IN_PATH + "/Assembly/NextDenovo/config/{sample}_nextDenovo.cfg",
    params:
        outDir = IN_PATH + "/Assembly/NextDenovo/{sample}",
    run:
        shell("realpath {input.fastq} > {output.fofn}")
        replace_config(input.config, output.config, output.fofn, params.outDir)



rule nextDenovo:
    input:
        config = IN_PATH + "/Assembly/NextDenovo/config/{sample}_nextDenovo.cfg",
    output:
        assembly = IN_PATH + "/Assembly/NextDenovo/{sample}/03.ctg_graph/nd.asm.fasta",
    threads:
        THREADS * ThreadFold
    log:
        IN_PATH + "/log/nextDenovo_{sample}.log"
    run:
        # shell("source activate nanovar2 && nextDenovo {input.config} > {log} 2>&1")
        cmd = "source activate Assembly3 && nextDenovo %s > %s 2>&1" % (input.config, log)
        print(cmd)
        os.system(cmd)


rule AssemblySummary:
    input:
        stats = expand(IN_PATH + "/Assembly/NextDenovo/{sample}/03.ctg_graph/nd.asm.fasta.stat", sample=SAMPLES),
    output:
        stats = IN_PATH + "/Assembly/NextDenovo/All_assembly_stats.xls",
    params:
        AssemblySummary = SCRIPT_DIR + "/AssemblySummary.py"
    run:
        Files = ",".join(input.stats)
        shell("python {params.AssemblySummary} --file {Files} --out {output.stats}")





"""
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
        sam = temp(IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT.sam"),
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
        contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT_racon.fasta",
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
        contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT_racon2.fasta",
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
        contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT_racon3.fasta",
    threads:
        THREADS * ThreadFold
    log:
        IN_PATH + "/log/racon3_{sample}.log"
    run:
        shell("racon --threads {threads} {input.fastq} {input.sam} {input.contig} > {output.contig} 2>{log}")
#################################################
"""




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
        contig = temp(IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_nextPolish_temp.fasta"),
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



##################### nextPolish 2 #############
rule ngsAlign3:
    input:
        contig =  IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_nextPolish1.fasta",
        R1 = rules.fastp.output.R1,
        R2 = rules.fastp.output.R2,
    output:
        bam = temp(IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_sort3.bam"),
    threads:
        THREADS
    run:
        shell("samtools faidx {input.contig}")
        shell("bwa index {input.contig}")
        shell("bwa mem -t {threads} {input.contig} {input.R1} {input.R2} | samtools view --threads {threads} -F 0x4 -b - | samtools fixmate -m --threads {threads}  - -| samtools sort -m 2g --threads {threads} -|samtools markdup --threads {threads} -r - {output.bam}")
        shell("samtools index -@ {threads} {output.bam}")



rule nextPolish3:
    input:
        contig =  IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_nextPolish1.fasta",
        bam = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_sort3.bam",
    output:
        contig = temp(IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_nextPolish2_temp.fasta"),
    params:
        nextpolish1 = "/home/wuzhikun/software/NextPolish/lib/nextpolish1.py",
    threads:
        THREADS
    run:
        shell("python {params.nextpolish1} -g {input.contig} -t 1 -p {threads} -s {input.bam} > {output.contig}")





rule ngsAlign4:
    input:
        contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_nextPolish2_temp.fasta",
        R1 = rules.fastp.output.R1,
        R2 = rules.fastp.output.R2,
    output:
        bam = temp(IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_sort4.bam"),
    threads:
        THREADS
    run:
        shell("samtools faidx {input.contig}")
        shell("bwa index {input.contig}")
        shell("bwa mem -t {threads} {input.contig} {input.R1} {input.R2} | samtools view --threads {threads} -F 0x4 -b - | samtools fixmate -m --threads {threads}  - -|samtools sort -m 2g --threads {threads} -|samtools markdup --threads {threads} -r - {output.bam}")
        shell("samtools index -@ {threads} {output.bam}")



rule nextPolish4:
    input:
        contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_nextPolish2_temp.fasta",
        bam = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_sort4.bam",
    output:
        contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_nextPolish2.fasta",
    params:
        nextpolish1 = "/home/wuzhikun/software/NextPolish/lib/nextpolish1.py",
    threads:
        THREADS
    run:
        shell("python {params.nextpolish1} -g {input.contig} -t 2 -p {threads} -s {input.bam} > {output.contig}")


##############################################


##################### nextPolish 3 #############
rule ngsAlign5:
    input:
        contig =  IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_nextPolish2.fasta",
        R1 = rules.fastp.output.R1,
        R2 = rules.fastp.output.R2,
    output:
        bam = temp(IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_sort5.bam"),
    threads:
        THREADS
    run:
        shell("samtools faidx {input.contig}")
        shell("bwa index {input.contig}")
        shell("bwa mem -t {threads} {input.contig} {input.R1} {input.R2} | samtools view --threads {threads} -F 0x4 -b - | samtools fixmate -m --threads {threads}  - -| samtools sort -m 2g --threads {threads} -|samtools markdup --threads {threads} -r - {output.bam}")
        shell("samtools index -@ {threads} {output.bam}")



rule nextPolish5:
    input:
        contig =  IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_nextPolish2.fasta",
        bam = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_sort5.bam",
    output:
        contig = temp(IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_nextPolish3_temp.fasta"),
    params:
        nextpolish1 = "/home/wuzhikun/software/NextPolish/lib/nextpolish1.py",
    threads:
        THREADS
    run:
        shell("python {params.nextpolish1} -g {input.contig} -t 1 -p {threads} -s {input.bam} > {output.contig}")





rule ngsAlign6:
    input:
        contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_nextPolish3_temp.fasta",
        R1 = rules.fastp.output.R1,
        R2 = rules.fastp.output.R2,
    output:
        bam = temp(IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_sort6.bam"),
    threads:
        THREADS
    run:
        shell("samtools faidx {input.contig}")
        shell("bwa index {input.contig}")
        shell("bwa mem -t {threads} {input.contig} {input.R1} {input.R2} | samtools view --threads {threads} -F 0x4 -b - | samtools fixmate -m --threads {threads}  - -|samtools sort -m 2g --threads {threads} -|samtools markdup --threads {threads} -r - {output.bam}")
        shell("samtools index -@ {threads} {output.bam}")



rule nextPolish6:
    input:
        contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_nextPolish3_temp.fasta",
        bam = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_sort6.bam",
    output:
        contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_nextPolish3.fasta",
    params:
        nextpolish1 = "/home/wuzhikun/software/NextPolish/lib/nextpolish1.py",
    threads:
        THREADS
    run:
        shell("python {params.nextpolish1} -g {input.contig} -t 2 -p {threads} -s {input.bam} > {output.contig}")


##############################################


####################### Rank assembly ##############
rule RankContigs:
    input:
        fasta = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT_racon1.fasta",
    output:
        fai = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT_racon1.fasta.fai",
        rank = IN_PATH + "/Assembly/Polish/{sample}/{sample}_contig_rank.txt",
    params:
        ScaffoldRank = SCRIPT_DIR + "/ScaffoldRank.py"
    run:
        shell("samtools faidx {input.fasta}")
        shell("python {params.ScaffoldRank} --fai {output.fai} --out {output.rank}")


def mergeRank(in_file, out_file):
    files = in_file.split(",")
    out_h = open(out_file, "w")
    for f in files:
        name = f.split("/")[-1].split('_')[0]
        in_h = open(f, "r")
        for line in in_h:
            lines = line.strip().split("\t")
            out_h.write("%s\t%s\t%s\n" % (name, lines[1], lines[2]))
        in_h.close()
    out_h.close()


rule MergeRank:
    input:
        rank = expand(IN_PATH + "/Assembly/Polish/{sample}/{sample}_contig_rank.txt", sample=SAMPLES),
    output:
        rank = IN_PATH + "/Assembly/Polish/Samples_contigs_rank.xls",
    run:
        Files = ",".join(input.rank)
        mergeRank(Files, output.rank)


rule MergeRankPlot:
    input:
        rank = IN_PATH + "/Assembly/Polish/Samples_contigs_rank.xls",
    output:
        pdf = IN_PATH + "/Assembly/Polish/Samples_contigs_rank.pdf",
    params:
        ScaffoldRank = SCRIPT_DIR + "/ScaffoldRank.R",
    run:
        cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --width 5 --height 4 " % (params.ScaffoldRank, input.rank, output.pdf)
        print(cmd)
        os.system(cmd)

######################################################




# ######################### Referenced scaffold #########################
rule RagTag:
    input:
        contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_nextPolish3.fasta",
    output:
        scaffold = IN_PATH + "/Assembly/RagTag/{sample}/ragtag.scaffold.fasta",
        stats = IN_PATH + "/Assembly/RagTag/{sample}/ragtag.scaffold.stats",
    params:
        RefGenome = config["RefGenome"],
        outDir = IN_PATH + "/Assembly/RagTag/{sample}",
    threads:
        THREADS
    log:
        IN_PATH + "/log/RagTag_{sample}.log",
    run:
        shell("ragtag.py scaffold {params.RefGenome} {input.contig}  -o {params.outDir}  -w -u -t {threads} --aligner minimap2 > {log} 2>&1")



def Merge_RagTag(in_file, out_file):
    out_h = open(out_file, "w")
    out_h.write("Sample\tplaced_sequences\tplaced_bp\tunplaced_sequences\tunplaced_bp\tgap_bp\tgap_sequences\n")
    files = in_file.split(",")
    for i in files:
        n = i.split("/")[-2]
        in_h = open(i, "r")
        header = in_h.readline().strip()
        line = in_h.readline().strip()
        out_h.write("%s\t%s\n" % (n, line))
        in_h.close()
    out_h.close()



rule MergeRagTag:
    input:
        stats = expand(IN_PATH + "/Assembly/RagTag/{sample}/ragtag.scaffold.stats", sample=SAMPLES),
    output:
        stats = IN_PATH + "/Assembly/RagTag/All_ragtag.scaffold.stats.xls",
    run:
        Files = ",".join(sorted(input.stats))
        Merge_RagTag(Files, output.stats)



rule RagTagGaps:
    input:
        assembly = IN_PATH + "/Assembly/RagTag/{sample}/ragtag.scaffold.fasta",
    output:
        gap = IN_PATH + "/Assembly/RagTagGaps/{sample}_scaffold_gap.txt",
    params:
        AssemblyGap = SCRIPT_DIR + "/AssemblyGap.py"
    threads:
        THREADS
    run:
        shell("python {params.AssemblyGap} --genome {input.assembly} --out {output.gap}")



rule RagTag2Ref:
    input:
        assembly = IN_PATH + "/Assembly/RagTag/{sample}/ragtag.scaffold.fasta",
    output:
        mashmap = IN_PATH + "/Assembly/RagTagMashmap/{sample}/mashmap.out.txt",
    params:
        RefGenome = config["RefGenome"],
        mashmap = "/home/wuzhikun/software/mashmap-Linux64-v2.0/mashmap",
        generateDotPlot = "/home/wuzhikun/software/mashmap-Linux64-v2.0/MashMap/scripts/generateDotPlot",
        outPrefix = IN_PATH + "/Assembly/RagTagMashmap/{sample}",
    threads:
        THREADS
    log:
        IN_PATH + "/log/RagTag2Ref_{sample}.log"
    run:
        shell("{params.mashmap} -r {params.RefGenome} -q {input.assembly} --perc_identity 95 --threads {threads} --segLength 50000 --filter_mode one-to-one --output {output.mashmap} > {log} 2>&1")
        shell("cd {params.outPrefix} && {params.generateDotPlot} png  large  {output.mashmap}")



rule ScaffoldStats:
    input:
        assembly = IN_PATH + "/Assembly/RagTag/{sample}/ragtag.scaffold.fasta",
    output:
        stats = IN_PATH + "/Assembly/RagTag/{sample}/ragtag.scaffold.stats.xls",
    run:
        shell("seqkit stats -aT {input.assembly} > {output.stats}")


rule ScaffoldStatsMerge:
    input:
        stats = expand(IN_PATH + "/Assembly/RagTag/{sample}/ragtag.scaffold.stats.xls", sample=SAMPLES),
    output:
        stats = IN_PATH + "/Assembly/RagTag/All_sample_ragtag.scaffold.stats.xls",
    run:
        Stats = sorted(input.stats)
        for i in range(len(Stats)):
            f = Stats[i]
            if i == 0:
                cmd = "cat %s > %s" % (f, output.stats)
            else:
                cmd = "sed '1d' %s >> %s " % (f, output.stats)
            os.system(cmd)

############################################################


##################### RagTag Revised #########
rule RevisedMap:
    input:
        assembly = IN_PATH + "/Assembly/RagTag/{sample}/ragtag.scaffold.revised.fasta",
    output:
        mashmap = IN_PATH + "/Assembly/RevisedMashmap/{sample}/mashmap.out.txt",
    params:
        RefGenome = config["RefGenome"],
        mashmap = "/home/wuzhikun/software/mashmap-Linux64-v2.0/mashmap",
        generateDotPlot = "/home/wuzhikun/software/mashmap-Linux64-v2.0/MashMap/scripts/generateDotPlot",
        outPrefix = IN_PATH + "/Assembly/RevisedMashmap/{sample}",
    threads:
        THREADS
    log:
        IN_PATH + "/log/RevisedMashmap_{sample}.log"
    run:
        shell("{params.mashmap} -r {params.RefGenome} -q {input.assembly} --perc_identity 95 --threads {threads} --segLength 50000 --filter_mode one-to-one --output {output.mashmap} > {log} 2>&1")
        shell("cd {params.outPrefix} && {params.generateDotPlot} png  large  {output.mashmap}")    

###############################################


########################## Scaffold ##############
rule Scaffold:
    input:
        # assembly = IN_PATH + "/Assembly/RagTag/{sample}/ragtag.scaffold.revised.fasta",
        assembly = IN_PATH + "/Assembly/RagTag/{sample}/ragtag.scaffold.fasta",
    output:
        assembly = IN_PATH + "/Assembly/Scaffold/{sample}.scaffold.fasta",
    threads:
        THREADS
    run:
        shell("seqkit seq -u -w 60 {input.assembly} | sed 's/_RagTag//g' > {output.assembly}")
#####################################################




##################### RagTag SV ####################
rule ont2scaffold:
    input:
        fastq = IN_PATH + "/clean/{sample}_ONT.fastq.gz",
        assembly = IN_PATH + "/Assembly/RagTag/{sample}/ragtag.scaffold.fasta",
    output:
        sam = temp(IN_PATH + "/Assembly/RagTagAlign/{sample}_ONTAlign.sam"),
        bam = IN_PATH + "/Assembly/RagTagAlign/{sample}_ONTAlign.bam",
    threads:
        THREADS
    log:
        IN_PATH + "/log/ont2scaffold_{sample}.log"
    run:
        shell("minimap2 --MD -a -x map-ont -t {threads} {input.assembly} {input.fastq} > {output.sam} 2>{log}")
        shell("samtools view -@ {threads} -b {output.sam} | samtools sort -  -@ {threads} -o {output.bam} > {log} 2>&1")
        shell("samtools index {output.bam}")


rule cuteSV:
    input:
        assembly = IN_PATH + "/Assembly/RagTag/{sample}/ragtag.scaffold.fasta",
        bam = IN_PATH + "/Assembly/RagTagAlign/{sample}_ONTAlign.bam",
    output:
        vcf = IN_PATH + "/Assembly/RagTagAlign/{sample}.cutesv.vcf",
    params:
        tempDir = IN_PATH + "/Assembly/RagTagAlign/{sample}",
    threads:
        THREADS
    log:
        IN_PATH + "/log/cuteSV_{sample}.log"
    run:
        if not os.path.exists(params.tempDir):
            os.makedirs(params.tempDir)
        cmd = "source activate nanovar && cuteSV  --threads %s --sample  %s --report_readid --min_support 3  --min_size 50 --min_siglength 50  --max_size -1 --max_cluster_bias_DEL 100  --diff_ratio_merging_DEL 0.3 --genotype %s  %s %s %s > %s  2>&1" % (threads, wildcards.sample, input.bam, input.assembly, output.vcf, params.tempDir, log)
        print(cmd)
        os.system(cmd)            


rule Sniffles:
    input:
        bam = IN_PATH + "/Assembly/RagTagAlign/{sample}_ONTAlign.bam",
    output:
        vcf = IN_PATH + "/Assembly/RagTagAlign/{sample}.sniffles.vcf",
    threads:
        THREADS 
    log:
        IN_PATH + "/log/Sniffles_{sample}.log"
    run:
        shell("sniffles --mapped_reads {input.bam} --vcf {output.vcf} --threads {threads}  --min_support 3 --min_length 50 --minmapping_qual 20 --num_reads_report -1 --min_seq_size 500  --genotype --report_BND --report-seq  >{log} 2>&1")


rule NanoSV:
    input:
        bam =IN_PATH + "/Assembly/RagTagAlign/{sample}_ONTAlign.bam",
    output:
        vcf = IN_PATH + "/Assembly/RagTagAlign/{sample}.NanoSV.vcf",
    threads:
        THREADS 
    params:
        nanosvConfig = "/home/wuzhikun/anaconda3/envs/NanoSV/lib/python3.6/site-packages/nanosv/config.ini"
    log:
        IN_PATH + "/log/NanoSV_{sample}.log"
    run:
        ### depth_support = False
        # shell("NanoSV --threads {threads} -c {params.nanosvConfig} -o {output.vcf} {input.bam} >{log} 2>&1")
        cmd = "source activate NanoSV && NanoSV --threads %s -c %s -o %s %s > %s 2>&1" % (threads, params.nanosvConfig, output.vcf, input.bam, log)
        print(cmd)
        os.system(cmd)


# rule nanovar:
#     input:
#         bam = IN_PATH + "/Assembly/RagTagAlign/{sample}_ONTAlign.bam",
#         assembly = IN_PATH + "/Assembly/RagTag/{sample}/ragtag.scaffold.fasta",
#     output:
#         total = IN_PATH + "/Assembly/RagTagAlign/{sample}/{sample}_ONT.nanovar.total.vcf",
#         vcf = IN_PATH + "/Assembly/RagTagAlign/{sample}/{sample}_ONT.nanovar.pass.vcf",
#     threads:
#         THREADS
#     params:
#         outdir = IN_PATH + "/SVCall/nanovar/{sample}",
#     log:
#         IN_PATH + "/log/nanovar_{sample}.log"
#     run:
#         ### nanovar: /home/litong/anaconda3/envs/nanovar/bin/nanovar
#         ### --data_type ont, pacbio-clr, pacbio-ccs
#         # shell("nanovar -r {params.RefGenome} -l {input.fastq} -t {threads} -o {params.outdir} >{log} 2>&1")
#         # cmd = "source activate NanoVar && nanovar -t %s --data_type ont --mincov 2 --minlen 50 %s %s %s > %s 2>&1" % (threads, input.bam, params.RefGenome, params.outdir, log)
#         cmd = "source activate nanovar && nanovar -t %s --data_type ont --mincov 3 --minlen 50 %s %s %s > %s 2>&1" % (threads, input.bam, input.assembly, params.outdir, log)
#         print(cmd)
#         os.system(cmd)

rule homoSV:
    input:
        vcf = IN_PATH + "/Assembly/RagTagAlign/{sample}.sniffles.vcf",
    output:
        homo = IN_PATH + "/Assembly/RagTagAlign/{sample}.sniffles.homoSV.txt",
    params:
        HomoSVRecord = SCRIPT_DIR + "/HomoSVRecord.py",
        AFThreshold = 0.9,
        REThreshold = 20,
    run:
        shell("python {params.HomoSVRecord} --vcf {input.vcf} --out {output.homo} --REThreshold {params.REThreshold} --AFThreshold {params.AFThreshold}")


#############################################################









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



