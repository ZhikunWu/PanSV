############# length evaluation #######
# rule contigLen:
#     input:
#         assembly = IN_PATH + "/Assembly/NextDenovo/{sample}/03.ctg_graph/nd.asm.fasta",
#     output:
#         fai = IN_PATH + "/Assembly/NextDenovo/{sample}/03.ctg_graph/nd.asm.fasta.fai",
#         length = IN_PATH + "/Evaluation/Length/NextDenovo/{sample}.contigs.txt",
#     params:
#         ScaffoldRank = ScriptDir + "/ScaffoldRank.py",
#     run:
#         shell("samtools faidx {input.assembly}")
#         shell("python {params.ScaffoldRank} --fai {output.fai} --out {output.length}")
    
# rule AllLength:
#     input:
#         length = expand(IN_PATH + "/Evaluation/Length/NextDenovo/{sample}.contigs.txt", sample=SAMPLES),
#     output:
#         length = IN_PATH + "/Evaluation/Length/NextDenovo/All_samples.contigs.txt",
#     run:
#         Files = sorted(input.length)
#         for i in range(len(Files)):
#             f = Files[i]
#             if i == 0:
#                 cmd = """ awk '{print "%s""\t"$2"\t"$3}' %s  > %s""" % (wildcards.sample, f, output.length)
#             else:
#                 cmd = """ awk '{print "%s""\t"$2"\t"$3}' %s  >> %s""" % (wildcards.sample, f, output.length)
#             os.system(cmd)

# rule AllLengthPlot:
#     input:
#         length = IN_PATH + "/Evaluation/Length/NextDenovo/All_samples.contigs.txt",
#     output:
#         pdf = IN_PATH + "/Evaluation/Length/NextDenovo/All_samples.contigs.pdf",d
#     params:
#         ScaffoldRankR = ScriptDir + "/ScaffoldRank.R",
#     run:
#         shell("Rscript {params.ScaffoldRankR} --input {input.length} --pdf {output.pdf} --width 5 --height 5")

######################################


########### quality value #######
rule merylCount1:
    input:
        R1 = IN_PATH + "/clean/{sample}_NGS.R1.fastq.gz",
    output:
        R1 = temp(IN_PATH + "/clean/{sample}_NGS.R1.fastq"),
        kmer = directory(IN_PATH + "/QualityValue/kmer/{sample}/{sample}_R1"),
    threads:
        THREADS
    log:
        IN_PATH + "/log/merylCount_{sample}_R1.log",
    run:
        shell("pigz -p {threads} -dc {input.R1} > {output.R1}")
        shell("meryl count k=19 output {output.kmer} {output.R1} > {log} 2>&1")


rule merylCount2:
    input:
        R2 = IN_PATH + "/clean/{sample}_NGS.R2.fastq.gz",
    output:
        R2 = temp(IN_PATH + "/clean/{sample}_NGS.R2.fastq"),
        kmer = directory(IN_PATH + "/QualityValue/kmer/{sample}/{sample}_R2"),
    threads:
        THREADS
    log:
        IN_PATH + "/log/merylCount_{sample}_R2.log",
    run:
        shell("pigz -p {threads} -dc {input.R2} > {output.R2}")
        shell("meryl count k=19 output {output.kmer} {output.R2} > {log} 2>&1")



rule merylSum:
    input:
        R1 = IN_PATH + "/QualityValue/kmer/{sample}/{sample}_R1",
        R2 = IN_PATH + "/QualityValue/kmer/{sample}/{sample}_R2",
    output:
        merge = directory(IN_PATH + "/QualityValue/kmer/{sample}/{sample}_merge"),
    threads:
        THREADS
    log:
        IN_PATH + "/log/merylSum_{sample}.log",
    run:
        Directory = " ".join([input.R1, input.R2])
        shell("meryl union-sum output {output.merge} {Directory} > {log} 2>&1")



rule RawQV:
    input:
        merge = IN_PATH + "/QualityValue/kmer/{sample}/{sample}_merge",
        assembly = IN_PATH + "/Assembly/NextDenovo/{sample}/03.ctg_graph/nd.asm.fasta",
    output:
        qv = IN_PATH + "/QualityValue/NextDenovo/{sample}/{sample}.qv",
    params:
        qv = IN_PATH + "/QualityValue/NextDenovo/{sample}",
    threads:
        THREADS
    log:
        IN_PATH + "/log/RawQV_{sample}.log",
    run:
        if not os.path.exists(params.qv):
            os.makedirs(params.qv)
        shell("cd {params.qv} && merqury.sh {input.merge}  {input.assembly}  {wildcards.sample} > {log} 2>&1")


rule ONTQV:
    input:
        merge = IN_PATH + "/QualityValue/kmer/{sample}/{sample}_merge",
        assembly = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ONT_racon3.fasta",
    output:
        qv = IN_PATH + "/QualityValue/ONTPolish/{sample}/{sample}.qv",
    params:
        qv = IN_PATH + "/QualityValue/ONTPolish/{sample}",
    threads:
        THREADS
    log:
        IN_PATH + "/log/ONTQV_{sample}.log",
    run:
        if not os.path.exists(params.qv):
            os.makedirs(params.qv)
        shell("cd {params.qv} && merqury.sh {input.merge}  {input.assembly}  {wildcards.sample} > {log} 2>&1")


rule NGSQV:
    input:
        merge = IN_PATH + "/QualityValue/kmer/{sample}/{sample}_merge",
        assembly = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_nextPolish3.fasta",
    output:
        qv = IN_PATH + "/QualityValue/NGSPolish/{sample}/{sample}.qv",
    params:
        qv = IN_PATH + "/QualityValue/NGSPolish/{sample}",
    threads:
        THREADS
    log:
        IN_PATH + "/log/NGSQV_{sample}.log",
    run:
        if not os.path.exists(params.qv):
            os.makedirs(params.qv)
        shell("cd {params.qv} && merqury.sh {input.merge}  {input.assembly}  {wildcards.sample} > {log} 2>&1")
########################################






################## BUSCO ############

rule BUSCO:
    input:
        assembly = IN_PATH + "/Assembly/Polish/{sample}/{sample}_ngs_nextPolish3.fasta",
    output:
        summary = IN_PATH + "/Evaluation/BUSCO/ngsPolish/{sample}/run_embryophyta_odb10/short_summary.txt",
    threads:
        THREADS 
    params:
        buscoDB = "/home/wuzhikun/database/BUSCO/version5/embryophyta_odb10",
        outDir = IN_PATH + "/Evaluation/BUSCO/ngsPolish",
        # species = "Vigna_unguiculata", 
        # evalue = 1e-05,
        buscoConfig = "/home/wuzhikun/Project/PanVigna/config/busco.config",
    log:
        IN_PATH + '/log/BUSCO_{sample}.log' 
    run:
        if not os.path.exists(params.outDir):
            os.mkdir(params.outDir)
        cmd = "source activate busco &&  cd %s  &&  busco --cpu %s -m genome -i %s  --out %s --out_path %s  -l %s  --config %s  --force --offline > %s 2>&1" % (params.outDir, threads, input.assembly, wildcards.sample, params.outDir, params.buscoDB, params.buscoConfig, log)
        print(cmd)
        os.system(cmd)

#########################################



######################### RepeatModeler ####################
rule BuildDatabase:
    input:
        #assembly = rules.RagTag.output.scaffold,
        assembly = IN_PATH + "/Assembly/Polish/{sample}/Racon/{sample}_ONT_polish_racon3.fasta",
    output:
        nhr = IN_PATH + "/Repeat/repeatModeler/{sample}/{sample}.nhr",
    params:
        outPrefix = IN_PATH + "/Repeat/repeatModeler/{sample}/{sample}",
        outDir = IN_PATH + "/Repeat/repeatModeler/{sample}",
    threads:
        THREADS
    log:
        IN_PATH + "/log/BuildDatabase_{sample}.log",
    run:
        if not os.path.exists(params.outDir):
            os.mkdir(params.outDir)
        # shell("cd {params.outDir} &&  BuildDatabase -name {params.outPrefix} -engine ncbi  {input.assembly} > {log} 2>&1")
        cmd = "source activate TE2 && cd %s &&  BuildDatabase -name %s -engine ncbi  %s > %s 2>&1" % (params.outDir, params.outPrefix, input.assembly, log)
        print(cmd)
        os.system(cmd)



rule repeatModeler:
    ### http://xuzhougeng.top/archives/Repeat-annotation-with-RepeatModeler
    input:
        nhr = IN_PATH + "/Repeat/repeatModeler/{sample}/{sample}.nhr",
    output:
        lib = IN_PATH + "/Repeat/repeatModeler/{sample}/{sample}-families.fa",
    params:
        outPrefix = IN_PATH + "/Repeat/repeatModeler/{sample}/{sample}",
        outDir = IN_PATH + "/Repeat/repeatModeler/{sample}",
    threads:
        THREADS
    log:
        IN_PATH + "/log/repeatModeler_{sample}.log",
    run:
        # shell("cd {params.outDir} &&  RepeatModeler -pa {threads} -database {params.outPrefix}  -engine ncbi > {log} 2>&1")
        cmd = "source activate TE2 && cd %s &&  RepeatModeler -pa %s -database %s  -engine ncbi > %s 2>&1" % (params.outDir, threads, params.outPrefix, log)
        print(cmd)
        os.system(cmd)

######################################################


#################### repeatMasker ##########
rule repeatMasker:
    input:
        #assembly = rules.RagTag.output.scaffold,
        assembly = IN_PATH + "/Assembly/Polish/{sample}/Racon/{sample}_ONT_polish_racon3.fasta",
        lib = IN_PATH + "/Repeat/repeatModeler/{sample}/{sample}-families.fa",
    output:
        # mask = IN_PATH + "/repeatMasker/{species}/{sample}/{sample}.fasta.masked",
        out = IN_PATH + "/Repeat/repeatModeler/{sample}/{sample}.fasta.out",
    threads:
        THREADS
    params:
        outPrefix = IN_PATH + "/Repeat/repeatMasker/{sample}",
    log:
        IN_PATH + "/log/repeatMasker_{sample}.log",
    run:
        # -libdir <string>
        #     Path to the RepeatMasker libraries directory.
        # shell("RepeatMasker -pa {threads} -e ncbi  -gff -poly -lib {input.lib}  -dir {params.outPrefix}  {input.assembly} > {log} 2>&1")
        cmd = "source activate TE2 && RepeatMasker -pa %s -e ncbi  -gff -poly -lib %s  -dir %s  %s > %s 2>&1" % (threads, input.lib, params.outPrefix, input.assembly, log)
        print(cmd)
        os.system(cmd)


# ##########################################################
