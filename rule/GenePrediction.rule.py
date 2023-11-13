############################################################################
# rule BuildDatabase:
#     input:
#         assembly = IN_PATH + "/Assembly/Scaffold/{sample}.scaffold.fasta",
#     output:
#         nhr = IN_PATH + "/Repeat/repeatModeler/{sample}/{sample}.nhr",
#     params:
#         outPrefix = IN_PATH + "/Repeat/repeatModeler/{sample}/{sample}",
#         outDir = IN_PATH + "/Repeat/repeatModeler/{sample}",
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/BuildDatabase_{sample}.log",
#     run:
#         if not os.path.exists(params.outDir):
#             os.mkdir(params.outDir)
#         shell("cd {params.outDir} &&  BuildDatabase -name {params.outPrefix} -engine ncbi  {input.assembly} > {log} 2>&1")



# rule repeatModeler:
#     ### http://xuzhougeng.top/archives/Repeat-annotation-with-RepeatModeler
#     input:
#         nhr = IN_PATH + "/Repeat/repeatModeler/{sample}/{sample}.nhr",
#     output:
#         lib = IN_PATH + "/Repeat/repeatModeler/{sample}/{sample}-families.fa",
#     params:
#         outPrefix = IN_PATH + "/Repeat/repeatModeler/{sample}/{sample}",
#         outDir = IN_PATH + "/Repeat/repeatModeler/{sample}",
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/repeatModeler_{sample}.log",
#     run:
#         shell("cd {params.outDir} &&  RepeatModeler -pa {threads} -database {params.outPrefix}  -engine ncbi > {log} 2>&1")




rule repeatMasker0:
    input:
        assembly = IN_PATH + "/Assembly/Scaffold/{sample}.scaffold.fasta",
        lib = IN_PATH + "/BACKUP/Repeat/Vigna_unguiculata_contig-families.fa",
    output:
        mask = IN_PATH + "/Repeat/repeatMasker/{sample}/{sample}.scaffold.fasta.masked",
        out = IN_PATH + "/Repeat/repeatMasker/{sample}/{sample}.scaffold.fasta.out",
    threads:
        THREADS
    params:
        outPrefix = IN_PATH + "/Repeat/repeatMasker/{sample}",
    log:
        IN_PATH + "/log/repeatMasker_{sample}.log",
    run:
        cmd = "source activate TE2 && RepeatMasker -xsmall -nolow -norna -html -gff -pa %s -e ncbi  -poly -lib %s  -dir %s  %s > %s 2>&1" % (threads, input.lib, params.outPrefix, input.assembly, log)
        print(cmd)
        os.system(cmd)
#############################################################



####################### HISAT alignment #####################
rule HisatIndex:
    input:
        assembly = IN_PATH + "/Assembly/Scaffold/{sample}.scaffold.fasta",
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
rule fastpRNA:
    input:
        R1 = IN_PATH + "/raw/{sample}_RNA.R1.fastq.gz",
        R2 = IN_PATH + "/raw/{sample}_RNA.R2.fastq.gz",
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
        IN_PATH + "/log/trim/{sample}.log", 
    run:
        shell("fastp --in1 {input.R1} --in2 {input.R2} --out1 {output.R1} --out2 {output.R2} --thread {threads} --compression 2 --length_required {params.length_required} --qualified_quality_phred {params.qualified_quality_phred}  --unqualified_percent_limit {params.unqualified_percent_limit} --cut_front --cut_tail --cut_window_size {params.cut_window_size} --cut_mean_quality {params.cut_mean_quality} >{log} 2>&1")


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


rule stringtieFasta:
    input:
        gtf = rules.stringtie.output.gtf,
        assembly = IN_PATH + "/Assembly/Scaffold/{sample}.scaffold.fasta",
    output:
        fasta = IN_PATH + "/GenePrediction/RNA/Stringtie/{sample}/{sample}_RNA_stringtie_exon.fasta",
    threads:
        THREADS
    log:
        IN_PATH + "/log/stringtieFasta_{sample}.log",
    run:
        shell("gffread  -w {output.fasta} -g {input.assembly} {input.gtf} > {log} 2>&1")


##############################################################################


################# De novo Assembly ####################
rule TrinityAssembly:
    ### https://www.jianshu.com/p/8518a23255f8
    input:
        R1 = IN_PATH + '/clean/{sample}.RNA.R1.fq.gz',
        R2 = IN_PATH + '/clean/{sample}.RNA.R2.fq.gz',
    output:
        fasta = IN_PATH + "/GenePrediction/RNA/trinity/{sample}/trinity/Trinity.fasta",
    params:
        outDir = IN_PATH + "/GenePrediction/RNA/trinity/{sample}/trinity",
    threads:
        THREADS
    log:
        IN_PATH + "/log/TrinityAssembly_{sample}.log",
    run:
        ### Trinity --seqType fq --max_memory 100G --CPU 20  --output /home/wuzhikun/Project/Vigna/RNA/NGS/Transcript/All_assembly/trinity  --left /home/wuzhikun/data/Vigna/RNA/NGS/RNA_NGS_All.R1.fastq    --right  /home/wuzhikun/data/Vigna/RNA/NGS/RNA_NGS_All.R2.fastq  > /home/wuzhikun/Project/Vigna/log/TrinityAssemblyAll.log  2>&1
        shell("Trinity --seqType fq --max_memory 100G --CPU {threads} --output {params.outDir} --left {input.R1}  --right {input.R2} > {log} 2>&1")

#########################################################


######################## non-redundant transcripts #######
rule TansMerge:
    input:
        fa1 = rules.TrinityAssembly.output.fasta,
        fa2 = rules.stringtieFasta.output.fasta,
    output:
        merge = IN_PATH + "/GenePrediction/RNA/TransMerge/{sample}/merged_transcript.fasta",
        fasta = IN_PATH + "/GenePrediction/RNA/TransMerge/{sample}/nonreduandant_transcript.fasta",
    threads:
        THREADS
    log:
        IN_PATH + "/log/TransMerge_{sample}.log",
    run:
        shell("cat {input.fa1} {input.fa2} > {output.merge}")
        shell("cd-hit-est -i {output.merge} -o {output.fasta} -c 0.95 -aS 0.95 -d 0 >{log} 2>&1")


###########################################################





############################## prediction ###############
######### using masked fasta ?
rule snap:
    input:
        assembly = IN_PATH + "/Assembly/Scaffold/{sample}.scaffold.fasta",
    output:
        protein = IN_PATH + "/GenePrediction/snap/{sample}/{sample}_protein.fa",
        transcript = IN_PATH + "/GenePrediction/snap/{sample}/{sample}_transcript.fa",
        # gff = IN_PATH + "/GenePrediction/snap/{sample}/{sample}_transcript.gff",
    params:
        snap_hmm = "/home/wuzhikun/anaconda3/envs/Assembly3/share/snap/HMM/A.thaliana.hmm",
    log:
        IN_PATH + "/log/snap_{sample}.log",
    threads:
        THREADS
    run:
        shell("snap -gff -quiet -lcmask {params.snap_hmm}  {input.assembly} -aa {output.protein} -tx  {output.transcript} > {output.gff} 2>{log}")        


rule augustus:
    input:
        assembly = IN_PATH + "/Assembly/Scaffold/{sample}.scaffold.fasta",
        # assembly = IN_PATH + "/Repeats/repeatModeler/Vigna_unguiculata_contig/Vigna_unguiculata_assembly.fasta.masked",
    output:
        predict = IN_PATH + "/GenePrediction/augustus/{sample}/{sample}_prediction.gff",
    log:
        IN_PATH + "/log/augustus_{sample}.log",
    threads:
        THREADS
    run:
        shell("augustus --species=arabidopsis --softmasking=1 --gff3=on  {input.assembly} > {output.predict} 2>{log}")
#############################################################



################### non-redundant protein ###############
rule MergeProtein:
    input:
        fa1 = "/home/wuzhikun/database/genome/arabidopsis/TAIR10_pep_20101214.fa",
        fa2 = "/home/wuzhikun/database/genome/Vigna_unguiculata/vigun.IT97K-499-35.gnm1.ann1.zb5D.protein.faa",
    output:
        merge = IN_PATH + "/GenePrediction/Protein/Merged_protein.fasta",
        fasta = IN_PATH + "/GenePrediction/Protein/Merged_protein_nonredundant.fasta",
    log:
        IN_PATH + "/log/MergeProtein.log",
    threads:
        THREADS
    run:
        shell("cat {input.fa1} {input.fa2} > {output.merge}")
        shell("cd-hit -i {output.merge} -o {output.fasta} -c 0.9 -aS 0.9 -d 0 >{log} 2>&1")

#########################################################






##################### TE ################
# rule EDTA:
#     input:
#         assembly = IN_PATH + "/Assembly/Scaffold/{sample}.scaffold.fasta",
#     output:
#         # log = IN_PATH + "/log/EDTA_{sample}.log",
#         TE = IN_PATH + "/Repeats/EDTA/{sample}/ragtag.scaffold.fasta.mod.EDTA.anno/ragtag.scaffold.fasta.mod.EDTA.TE.fa.cln",
#     params:
#         EDTA = config["EDTA"],
#         CDSSeq = config["CDSSeq"],
#         outDir = IN_PATH + "/Repeats/EDTA/{sample}",
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/EDTA_{sample}.log",
#     run:
#         shell("cd {params.outDir} && perl {params.EDTA} --step all  --genome  {input.assembly}  --species others --cds {params.CDSSeq} --overwrite 1 --sensitive 1 --anno 1 --evaluate 1 --threads {threads} > {log} 2>&1")

########################################