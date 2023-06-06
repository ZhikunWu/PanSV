
######################### Assembly ############################
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


# $ cat NextDenovo.cfg 

# [General]
# job_type = local # local, slurm, sge, pbs, lsf
# job_prefix = nextDenovo
# task = all # all, correct, assemble
# rewrite = yes # yes/no
# deltmp = yes 
# parallel_jobs = 10 # number of tasks used to run in parallel
# input_type = raw # raw, corrected
# read_type = ont # clr, ont, hifi
# input_fofn = input.fofn
# workdir = 01_rundir

# [correct_option]
# read_cutoff = 1k
# genome_size = 0.52g#1g # estimated genome size
# sort_options = -m 20g -t 15
# minimap2_options_raw = -t 6
# pa_correction = 4 # number of corrected tasks used to run in parallel, each corrected task requires ~TOTAL_INPUT_BASES/4 bytes of memory usage.
# correction_options = -p 15

# [assemble_option]
# minimap2_options_cns = -t 8 
# nextgraph_options = -a 1



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




# rule nextStats:
#     ### https://github.com/raymondkiu/sequence-stats
#     input:
#         assembly = IN_PATH + "/Assembly/NextDenovo/Vigna_unguiculata/03.ctg_graph/nd.asm.fasta",
#     output:
#         stat = IN_PATH + "/Assembly/NextDenovo/Vigna_unguiculata/03.ctg_graph/nd.asm_each_stat.txt",
#     run:
#         shell("sequence-stats -c  {input.assembly} > {output.stat}")


# rule nextPlot:
#     input:
#         stat = IN_PATH + "/Assembly/NextDenovo/Vigna_unguiculata/03.ctg_graph/nd.asm_each_stat.txt",
#     output:
#         pdf = IN_PATH + "/Assembly/NextDenovo/Vigna_unguiculata/03.ctg_graph/nd.asm_each_stat.pdf",
#     params:
#         TreeMap = SCRIPT_DIR + "/TreeMap.R",
#         width = 5,
#         height = 4,
#     log:
#         IN_PATH + "/log/nextPlot.log"
#     run:
#         # shell("Rscript {params.TreeMap} --input {input.stat} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")
#         cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --width %s --height %s > %s 2>&1" % (params.TreeMap, input.stat, output.pdf, params.width, params.height, log)
#         print(cmd)
#         os.system(cmd)

##############################################################


########################### recon Polish ############################

rule mm2_ont:
    input:
        fastq = rules.SeqFilt.output.fastq,
        contig = IN_PATH + "/Assembly/NextDenovo/{sample}/03.ctg_graph/nd.asm.fasta",
    output:
        sam = temp(IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT_align1.sam"),
    threads:
        THREADS * ThreadFold
    log:
        IN_PATH + "/log/mm2_ont_{sample}.log"
    run:
        shell("minimap2 --MD -a -x map-ont -t {threads} {input.contig} {input.fastq} > {output.sam} 2>{log}")


rule racon:
    input:
        sam = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT_align1.sam",
        contig = IN_PATH + "/Assembly/NextDenovo/{sample}/03.ctg_graph/nd.asm.fasta",
        fastq = rules.SeqFilt.output.fastq,
    output:
        contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT_racon1.fasta",
    threads:
        THREADS * ThreadFold
    log:
        IN_PATH + "/log/racon_{sample}.log"
    run:
        shell("racon --threads {threads} {input.fastq} {input.sam} {input.contig} > {output.contig} 2>{log}")





rule mm2_ont2:
    input:
        fastq = rules.SeqFilt.output.fastq,
        contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT_racon1.fasta",
    output:
        sam = temp(IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT_align2.sam"),
    threads:
        THREADS * ThreadFold
    log:
        IN_PATH + "/log/mm2_ont2_{sample}.log"
    run:
        shell("minimap2 --MD -a -x map-ont -t {threads} {input.contig} {input.fastq} > {output.sam} 2>{log}")


rule racon2:
    input:
        sam = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT_align2.sam",
        contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT_racon1.fasta",
        fastq = rules.SeqFilt.output.fastq,
    output:
        contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT_racon2.fasta",
    threads:
        THREADS * ThreadFold
    log:
        IN_PATH + "/log/racon2_{sample}.log"
    run:
        shell("racon --threads {threads} {input.fastq} {input.sam} {input.contig} > {output.contig} 2>{log}")




rule mm2_ont3:
    input:
        fastq = rules.SeqFilt.output.fastq,
        contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT_racon2.fasta",
    output:
        sam = temp(IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT_align3.sam"),
    threads:
        THREADS * ThreadFold
    log:
        IN_PATH + "/log/mm2_ont3_{sample}.log"
    run:
        shell("minimap2 --MD -a -x map-ont -t {threads} {input.contig} {input.fastq} > {output.sam} 2>{log}")



rule racon3:
    input:
        sam = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT_align3.sam",
        contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT_racon2.fasta",
        fastq = rules.SeqFilt.output.fastq,
    output:
        contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT_racon3.fasta",
    threads:
        THREADS * ThreadFold
    log:
        IN_PATH + "/log/racon3_{sample}.log"
    run:
        shell("racon --threads {threads} {input.fastq} {input.sam} {input.contig} > {output.contig} 2>{log}")




###############################################################



########## pilon ##############
# rule bwa:
#     input:
#         contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT_racon3.fasta",
#         R1 = rules.DNAQC.output.R1,
#         R2 = rules.DNAQC.output.R2,
#     output:
#         bam = temp(IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs.bam"),
#     params:
#         sa = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT_racon3.fasta.sa",
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/bwa_{sample}.log"
#     run:
#         if not os.path.exists(params.sa):
#             cmd = "bwa index %s" % (input.contig)
#             os.system(cmd)
#         # shell("bwa mem -t {threads} {input.contig} {input.R1} {input.R2} | samtools view -@ {threads} -F 0x4 -b - | samtools sort - -m 5g -@ {threads} -o {output.bam} > {log} 2>&1")
#         shell("bwa mem -t {threads} {input.contig} {input.R1} {input.R2} | samtools view -@ {threads}  -b - | samtools sort -m 5g -@ {threads} -o {output.bam} - > {log} 2>&1")
#         shell("samtools index -@ {threads} {output.bam}")




# rule pilon:
#     input:
#         contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT_racon3.fasta",
#         bam = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs.bam",
#     output:
#         contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_pilon1.fasta",
#     threads:
#         THREADS 
#     params:
#         pilon = config["pilon"],
#         outdir = IN_PATH + "/Assembly/Polish/{sample}",
#         memory = "-Xmx100G",
#         mindepth = 20,
#     log:
#         IN_PATH + "/log/pilon_{sample}.log"
#     run:
#         shell("java  {params.memory}  -jar {params.pilon}  --genome {input.contig}  --bam {input.bam}  --outdir {params.outdir} --output {wildcards.sample}_ngs_pilon1 --threads {threads} --mindepth {params.mindepth} --changes --vcf --diploid  >{log} 2>&1")




# rule bwa2:
#     input:
#         contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_pilon1.fasta",
#         R1 = rules.DNAQC.output.R1,
#         R2 = rules.DNAQC.output.R2,
#     output:
#         bam = temp(IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs2.bam"),
#     params:
#         sa = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_pilon1.fasta.sa",
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/bwa2_{sample}.log"
#     run:
#         if not os.path.exists(params.sa):
#             cmd = "bwa index %s" % (input.contig)
#             os.system(cmd)
#         # shell("bwa mem -t {threads} {input.contig} {input.R1} {input.R2} | samtools view -@ {threads} -F 0x4 -b - | samtools sort - -m 5g -@ {threads} -o {output.bam} > {log} 2>&1")
#         shell("bwa mem -t {threads} {input.contig} {input.R1} {input.R2} | samtools view -@ {threads}  -b - | samtools sort -m 5g -@ {threads} -o {output.bam} - > {log} 2>&1")
#         shell("samtools index -@ {threads} {output.bam}")




# rule pilon2:
#     input:
#         contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_pilon1.fasta",
#         bam = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs2.bam",
#     output:
#         contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_pilon2.fasta",
#     threads:
#         THREADS 
#     params:
#         pilon = config["pilon"],
#         outdir = IN_PATH + "/Assembly/Polish/{sample}",
#         memory = "-Xmx100G",
#         mindepth = 20,
#     log:
#         IN_PATH + "/log/pilon2_{sample}.log"
#     run:
#         shell("java  {params.memory}  -jar {params.pilon}  --genome {input.contig}  --bam {input.bam}  --outdir {params.outdir} --output {wildcards.sample}_ngs_pilon2 --threads {threads} --mindepth {params.mindepth} --changes --vcf --diploid  >{log} 2>&1")





# rule bwa3:
#     input:
#         contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_pilon2.fasta",
#         R1 = rules.DNAQC.output.R1,
#         R2 = rules.DNAQC.output.R2,
#     output:
#         bam = temp(IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs3.bam"),
#     params:
#         sa = IIN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_pilon2.fasta.sa",
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/bwa3_{sample}.log"
#     run:
#         if not os.path.exists(params.sa):
#             cmd = "bwa index %s" % (input.contig)
#             os.system(cmd)
#         # shell("bwa mem -t {threads} {input.contig} {input.R1} {input.R2} | samtools view -@ {threads} -F 0x4 -b - | samtools sort - -m 5g -@ {threads} -o {output.bam} > {log} 2>&1")
#         shell("bwa mem -t {threads} {input.contig} {input.R1} {input.R2} | samtools view -@ {threads}  -b - | samtools sort -m 5g -@ {threads} -o {output.bam} - > {log} 2>&1")
#         shell("samtools index -@ {threads} {output.bam}")




# rule pilon3:
#     input:
#         contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_pilon2.fasta",
#         bam = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs.bam",
#     output:
#         contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_pilon3.fasta",
#     threads:
#         THREADS 
#     params:
#         pilon = config["pilon"],
#         outdir = IN_PATH + "/Assembly/Polish/{sample}",
#         memory = "-Xmx100G",
#         mindepth = 20,
#     log:
#         IN_PATH + "/log/pilon3_{sample}.log"
#     run:
#         shell("java  {params.memory}  -jar {params.pilon}  --genome {input.contig}  --bam {input.bam}  --outdir {params.outdir} --output {wildcards.sample}_ngs_pilon3 --threads {threads} --mindepth {params.mindepth} --changes --vcf --diploid  >{log} 2>&1")
###############################



###################### nextPolish #############
# rule ngsAlign1:
#     input:
#         contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT_racon3.fasta",
#         # R1 = 
#         # R2 = 
#     output:
#         bam = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_sort1.bam",
#     threads:
#         THREADS
#     run:
#         shell("samtools faidx {input.contig}")
#         shell("bwa index {input.contig}")
#         shell("bwa mem -t {threads} {input.contig} {input.R1} {input.R2} | samtools view --threads {threads} -F 0x4 -b - | samtools fixmate -m --threads {threads}  - -|samtools sort -m 2g --threads {threads} -|samtools markdup --threads {threads} -r - {output.bam}")
#         shell("samtools index -@ {threads} {output.bam}")



# rule nextPolish1:
#     input:
#         contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT_racon3.fasta",
#         bam = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_sort1.bam",
#     output:
#         contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_nextPolish_temp.fasta",
#     params:
#         nextpolish1 = "/home/wuzhikun/software/NextPolish/lib/nextpolish1.py",
#     threads:
#         THREADS
#     run:
#         shell("python {params.nextpolish1} -g {input.contig} -t 1 -p {threads} -s {input.bam} > {output.contig}")





# rule ngsAlign2:
#     input:
#         contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_nextPolish_temp.fasta",
#         # R1 = 
#         # R2 = 
#     output:
#         bam = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_sort2.bam",
#     threads:
#         THREADS
#     run:
#         shell("samtools faidx {input.contig}")
#         shell("bwa index {input.contig}")
#         shell("bwa mem -t {threads} {input.contig} {input.R1} {input.R2} | samtools view --threads {threads} -F 0x4 -b - | samtools fixmate -m --threads {threads}  - -|samtools sort -m 2g --threads {threads} -|samtools markdup --threads {threads} -r - {output.bam}")
#         shell("samtools index -@ {threads} {output.bam}")



# rule nextPolish2:
#     input:
#         contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_nextPolish_temp.fasta",
#         bam = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_sort2.bam",
#     output:
#         contig = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_nextPolish1.fasta",
#     params:
#         nextpolish1 = "/home/wuzhikun/software/NextPolish/lib/nextpolish1.py",
#     threads:
#         THREADS
#     run:
#         shell("python {params.nextpolish1} -g {input.contig} -t 2 -p {threads} -s {input.bam} > {output.contig}")


###############################################




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
