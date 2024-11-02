'''
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
        assembly = IN_PATH + "/Scaffold/{sample}.genome.fasta",
        lib = IN_PATH + "/Repeats/repeatModeler/Merge/merge.families.cluster.fa",
    output:
        mask = IN_PATH + "/Repeats/repeatMasker/{sample}/{sample}.genome.fasta.masked",
        out = IN_PATH + "/Repeats/repeatMasker/{sample}/{sample}.genome.fasta.out",
        gff = IN_PATH + "/Repeats/repeatMasker/{sample}/{sample}.genome.fasta.out.gff",
    threads:
        THREADS
    params:
        outPrefix = IN_PATH + "/Repeats/repeatMasker/{sample}",
    log:
        IN_PATH + "/log/repeatMasker_{sample}.log",
    run:
        cmd = "source activate TE2 && RepeatMasker -xsmall -nolow -norna -html -gff -pa %s -e ncbi  -poly -lib %s  -dir %s  %s > %s 2>&1" % (threads, input.lib, params.outPrefix, input.assembly, log)
        print(cmd)
        os.system(cmd)


#############################################################



'''
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
'''


rule HisatIndex:
    input:
        #assembly = IN_PATH + "/Assembly/Scaffold/{sample}.scaffold.fasta",
        assembly = IN_PATH + "/Scaffold/{sample}.genome.fasta",
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
        bam = IN_PATH + "/GenePrediction/RNA/hisat/{sample}/ngs_RNA_hisat.bam",
    threads:
        THREADS
    log:
        IN_PATH + "/log/hisatbam_{sample}.log", 
    run:
        shell("samtools view -@ {threads} -b {input.sam} | samtools sort -  -@ {threads} -o {output.bam} > {log} 2>&1")
        shell("samtools index {output.bam}")


##############################################################################


rule BAMStats4:
    input:
        bam = IN_PATH + "/GenePrediction/RNA/hisat/{sample}/ngs_RNA_hisat.bam",
    output:
        stat = IN_PATH + '/GenePrediction/RNA/hisat/{sample}_ngs_RNA_stats.xls',
    threads:
        THREADS
    log:
        IN_PATH + "/log/BAMStats3_{sample}.log"
    run:
        shell("samtools flagstat --threads {threads}  {input.bam} > {output.stat} 2>{log}")




####################### merge gtf ################
rule TrinityGuideAssembly:
    input:
        bam = IN_PATH + "/GenePrediction/RNA/hisat/{sample}/ngs_RNA_hisat.bam",
    output:
        fa = IN_PATH + "/GenePrediction/RNA/Trinity/guideAssembly/{sample}/trinity_out_dir/Trinity-GG.fasta",
    params:
        outDir = IN_PATH + "/GenePrediction/RNA/Trinity/guideAssembly/{sample}",
    threads:
        THREADS
    log:
        IN_PATH + "/log/TrinityGuideAssembly_{sample}.log",
    run:
        ### Trinity --genome_guided_bam /home/wuzhikun/Project/Vigna/RNA/NGS/mapping/hisat/all_ngs_RNA_hisat_sort.bam   --genome_guided_max_intron 10000 --max_memory 100G --CPU 20 > /home/wuzhikun/Project/Vigna/log/TrinityGuideAssembly.log 2>&1
        if not os.path.exists(params.outDir):
            os.makedirs(params.outDir)
        shell("cd {params.outDir} && Trinity --genome_guided_bam {input.bam} --genome_guided_max_intron 10000 --max_memory 50G --CPU {threads} > {log} 2>&1")



rule stringtie:
    input:
        bam = IN_PATH + "/GenePrediction/RNA/hisat/{sample}/ngs_RNA_hisat.bam",
    output:
        gtf = IN_PATH + "/GenePrediction/RNA/Stringtie/{sample}/{sample}_RNA_stringtie.gtf",
    threads:
        THREADS
    log:
        IN_PATH + "/log/stringtie_{sample}.log",
    run:
        shell("stringtie {input.bam} -p {threads} -o {output.gtf} > {log} 2>&1")

rule gtf2gff:
    input:
        gtf = IN_PATH + "/GenePrediction/RNA/Stringtie/{sample}/{sample}_RNA_stringtie.gtf",
    output:
        gff = IN_PATH + "/GenePrediction/RNA/Stringtie/{sample}/{sample}_RNA_stringtie.gff",
    run:
        shell("gffread -o {output.gff} {input.gtf}")





rule stringtieFasta:
    input:
        gtf = rules.stringtie.output.gtf,
        assembly = IN_PATH + "/Scaffold/{sample}.genome.fasta",
    output:
        fasta = IN_PATH + "/GenePrediction/RNA/Stringtie/{sample}/{sample}_RNA_stringtie_exon.fasta",
    threads:
        THREADS
    log:
        IN_PATH + "/log/stringtieFasta_{sample}.log",
    run:
        shell("gffread  -w {output.fasta} -g {input.assembly} {input.gtf} > {log} 2>&1")


rule cufflinks:
    input:
        bam = IN_PATH + "/GenePrediction/RNA/hisat/{sample}/ngs_RNA_hisat.bam",
    output:
        gtf = IN_PATH + "/GenePrediction/RNA/cufflinks/{sample}/transcripts.gtf",
    params:
        outDir = IN_PATH + "/GenePrediction/RNA/cufflinks/{sample}",
    threads:
        THREADS
    log:
        IN_PATH + "/log/cufflinks_{sample}.log",
    run:
        # shell("cufflinks -p {threads} --library-type fr-unstranded  -o {params.outDir} {input.bam} > {log} 2>&1")
        cmd = "source activate Assembly && cufflinks -p %s --library-type fr-unstranded  -o %s %s > %s 2>&1" % (threads, params.outDir, input.bam, log)
        print(cmd)
        os.system(cmd)


rule gtf2fasta1:
    input:
        assembly = IN_PATH + "/Scaffold/{sample}.genome.fasta",
        gtf = IN_PATH + "/GenePrediction/RNA/cufflinks/{sample}/transcripts.gtf",
    output:
        fasta = IN_PATH + "/GenePrediction/RNA/cufflinks/{sample}/transcripts.fasta",
    run:
        shell("gffread {input.gtf} -g {input.assembly} -w {output.fasta}")




rule scallop:
    input:
        bam = IN_PATH + "/GenePrediction/RNA/hisat/{sample}/ngs_RNA_hisat.bam",
    output:
        gtf = IN_PATH + "/GenePrediction/RNA/scallop/{sample}_transcripts.gtf",
    threads:
        THREADS
    log:
        IN_PATH + "/log/scallop_{sample}.log",
    run:
        shell("scallop -i {input.bam} -o {output.gtf} > {log} 2>&1")



rule gtf2fasta2:
    input:
        assembly = IN_PATH + "/Scaffold/{sample}.genome.fasta",
        gtf = IN_PATH + "/GenePrediction/RNA/scallop/{sample}_transcripts.gtf",
    output:
        fasta = IN_PATH + "/GenePrediction/RNA/scallop/{sample}_transcripts.fasta",
    run:
        shell("gffread {input.gtf} -g {input.assembly} -w {output.fasta}")





rule cuffmerge:
    input:
        stringtie = IN_PATH + "/GenePrediction/RNA/Stringtie/{sample}/{sample}_RNA_stringtie.gtf",
        cufflinks = IN_PATH + "/GenePrediction/RNA/cufflinks/{sample}/transcripts.gtf",
        scallop = IN_PATH + "/GenePrediction/RNA/scallop/{sample}_transcripts.gtf",
        assembly = IN_PATH + "/Scaffold/{sample}.genome.fasta",
    output:
        transList = IN_PATH + "/GenePrediction/RNA/cuffmerge/{sample}/assembly_merge_list.txt",
        gtf = IN_PATH + "/GenePrediction/RNA/cuffmerge/{sample}/merged.gtf",
    params:
        outDir = IN_PATH + "/GenePrediction/RNA/cuffmerge/{sample}",
    log:
        IN_PATH + "/log/cuffmerge_{sample}.log",
    run:
        if not os.path.exists(params.outDir):
            os.makedirs(params.outDir)
        shell("echo {input.stringtie} > {output.transList}")
        shell("echo {input.cufflinks} >> {output.transList}")
        shell("echo {input.scallop} >> {output.transList}")
        cmd = "source activate Assembly &&  cuffmerge -s %s -p %s -o %s  %s > %s 2>&1" % (input.assembly, threads, params.outDir, output.transList, log)
        print(cmd)
        os.system(cmd)


rule gtf2fasta:
    input:
        assembly = IN_PATH + "/Scaffold/{sample}.genome.fasta",
        gtf = IN_PATH + "/GenePrediction/RNA/cuffmerge/{sample}/merged.gtf",
    output:
        fasta = IN_PATH + "/GenePrediction/RNA/cuffmerge/{sample}/merged.fasta",
    run:
        shell("gffread {input.gtf} -g {input.assembly} -w {output.fasta}")


#########################################################



############## TransDecoder ###################

rule extractSeq:
    input:
        gtf = IN_PATH + "/GenePrediction/RNA/Stringtie/{sample}/{sample}_RNA_stringtie.gtf",
        assembly = IN_PATH + "/Scaffold/{sample}.genome.fasta",
    output:
        fa = IN_PATH + "/GenePrediction/RNA/Stringtie/{sample}/{sample}_RNA_stringtie.fasta",
    run:
        shell("/home/wuzhikun/anaconda3/envs/Anno/bin/gtf_genome_to_cdna_fasta.pl {input.gtf}  {input.assembly} > {output.fa}")


rule gtf2gff3:
    input:
        gtf = IN_PATH + "/GenePrediction/RNA/Stringtie/{sample}/{sample}_RNA_stringtie.gtf",
    output:
        gff3 = IN_PATH + "/GenePrediction/RNA/Stringtie/{sample}/{sample}_RNA_stringtie.gff3",
    run:
        shell("/home/wuzhikun/anaconda3/envs/Anno/bin/gtf_to_alignment_gff3.pl {input.gtf}  > {output.gff3}")

rule extractORF:
    input:
        fa = IN_PATH + "/GenePrediction/RNA/Stringtie/{sample}/{sample}_RNA_stringtie.fasta",
    output:
        gff3 = IN_PATH + "/GenePrediction/RNA/Stringtie/{sample}/{sample}_RNA_stringtie.fasta.transdecoder_dir/longest_orfs.gff3",
    params:
        outDir = IN_PATH + "/GenePrediction/RNA/Stringtie/{sample}",
    log:
        IN_PATH + "/log/extractORF_{sample}.log"
    run:
        shell("/home/wuzhikun/anaconda3/envs/Anno/bin/TransDecoder.LongOrfs -t {input.fa} --output_dir {params.outDir} > {log} 2>&1")

rule ORF:
    input:
        fa = IN_PATH + "/GenePrediction/RNA/Stringtie/{sample}/{sample}_RNA_stringtie.fasta",
        transgff = IN_PATH + "/GenePrediction/RNA/Stringtie/{sample}/{sample}_RNA_stringtie.gff3",
        orfgff = IN_PATH + "/GenePrediction/RNA/Stringtie/{sample}/{sample}_RNA_stringtie.fasta.transdecoder_dir/longest_orfs.gff3",
    output:
        gff3 = IN_PATH + "/GenePrediction/RNA/Stringtie/{sample}/{sample}_RNA_stringtie.fasta.transdecoder.genome.gff3",
    run:
        shell("/home/wuzhikun/anaconda3/envs/Anno/bin/cdna_alignment_orf_to_genome_orf.pl  {input.orfgff}  {input.transgff} {input.fa} > {output.gff3}")

#############################################



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
        shell("Trinity --seqType fq --max_memory 50G --CPU {threads} --output {params.outDir} --left {input.R1}  --right {input.R2} > {log} 2>&1")

#########################################################


######################## non-redundant transcripts #######
rule TansMerge:
    input:
        fa1 = IN_PATH + "/GenePrediction/RNA/Stringtie/{sample}/{sample}_RNA_stringtie_exon.fasta",
        fa2 = IN_PATH + "/GenePrediction/RNA/cufflinks/{sample}/transcripts.fasta",
        fa3 = IN_PATH + "/GenePrediction/RNA/scallop/{sample}_transcripts.fasta",
    output:
        merge = IN_PATH + "/GenePrediction/RNA/TransMerge/{sample}/merged_transcript.fasta",
        fasta = IN_PATH + "/GenePrediction/RNA/TransMerge/{sample}/nonreduandant_transcript.fasta",
    threads:
        THREADS
    log:
        IN_PATH + "/log/TransMerge_{sample}.log",
    run:
        shell("cat {input.fa1} {input.fa2} {input.fa3} > {output.merge}")
        shell("cd-hit-est -i {output.merge} -o {output.fasta} -c 0.95 -aS 0.95 -d 0 >{log} 2>&1")


###########################################################


rule homo:
    input:
        assembly = IN_PATH + "/Scaffold/{sample}.genome.fasta",
    output:
        homo = IN_PATH + "/GenePrediction/homo/miniprot/{sample}_align.gff",
    params:
        protein = "/home/wuzhikun/data/Vigna/Protein/Faba_10genome_protein_cdhit_0.8",
    threads:
        THREADS
    run:
        shell("miniprot -t {threads}  --gff {input.assembly}  {params.protein} > {output.homo}")



rule gffFilt:
    input:
        gff = IN_PATH + "/GenePrediction/homo/miniprot/{sample}_align.gff",
    output:
        gff = IN_PATH + "/GenePrediction/homo/miniprot/{sample}_align_filt.gff",
    params:
        miniProtIdenFilter = SCRIPT_DIR + "/miniProtIdenFilter.py"
    run:
        shell("python {params.miniProtIdenFilter} --gff {input.gff} --out {output.gff} --idenThreshold 0.8") 



rule gffAddID:
    input:
        gff = IN_PATH + "/GenePrediction/homo/miniprot/{sample}_align.gff",
    output:
        gff = IN_PATH + "/GenePrediction/homo/miniprot/{sample}_align_addID.gff",
    params:
        gffAddID = SCRIPT_DIR + "/gffAddID.py",
    run:
        shell("python {params.gffAddID} --gff {input.gff} --out {output.gff}")





############################## prediction ###############
######### using masked fasta ?
rule snap:
    input:
        assembly = IN_PATH + "/Scaffold/{sample}.genome.fasta",
    output:
        protein = IN_PATH + "/GenePrediction/snap/{sample}/{sample}_protein.fa",
        transcript = IN_PATH + "/GenePrediction/snap/{sample}/{sample}_transcript.fa",
        gff = IN_PATH + "/GenePrediction/snap/{sample}/{sample}_transcript.gff",
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
        assembly = IN_PATH + "/Scaffold/{sample}.genome.fasta",
        # assembly = IN_PATH + "/Repeats/repeatModeler/Vigna_unguiculata_contig/Vigna_unguiculata_assembly.fasta.masked",
    output:
        predict = IN_PATH + "/GenePrediction/augustus/{sample}_augustus.gff",
    log:
        IN_PATH + "/log/augustus_{sample}.log",
    threads:
        THREADS
    run:
        #shell("augustus --species=arabidopsis --softmasking=1 --gff3=on  {input.assembly} > {output.predict} 2>{log}")
        cmd = "source activate Braker && augustus  --species=arabidopsis --gff3=on %s > %s" % (input.assembly, output.predict)
        print(cmd)
        os.system(cmd)





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



'''

rule EVidenceModeler:
    input:
        genome = IN_PATH + "/Scaffold/{sample}.genome.fasta",
        augustus = IN_PATH + "/GenePrediction/augustus/{sample}_augustus.gff",
        stringtie = IN_PATH + "/GenePrediction/RNA/Stringtie/{sample}/{sample}_RNA_stringtie.fasta.transdecoder.genome.gff3",
        miniprot = IN_PATH + "/GenePrediction/homo/miniprot/{sample}_align_addID.gff",
        #repeat = IN_PATH + "/Repeats/repeatMasker/{sample}/{sample}.genome.fasta.out.gff"
    output:
        evidence = IN_PATH + "/GenePrediction/evidencemodeler/{sample}/{sample}.evm.out",
        EVM = IN_PATH + "/GenePrediction/evidencemodeler/{sample}/{sample}.EVM.gff3",
    threads:
        THREADS
    params:
        weight = IN_PATH + "/config/Trans_weight.txt",
        outDir = IN_PATH + "/GenePrediction/evidencemodeler/{sample}",
    log:
        IN_PATH + "/log/EVidenceModeler_{sample}.log"
    run:
        ### EVidenceModeler --sample_id CPG11 --segmentSize 100000 --overlapSize 10000  --genome /home/wuzhikun/Project/PanSV/Scaffold/CPG11.genome.fasta  --weights /home/wuzhikun/Project/BAssembly/pipeline/evidencemodeler/Trans_weight.txt --gene_predictions  /home/wuzhikun/Project/PanSV/GenePrediction/augustus/CPG11_augustus.gff --protein_alignments  /home/wuzhikun/Project/PanSV/GenePrediction/homo/miniprot/CPG11_align.gff --transcript_alignments  /home/wuzhikun/Project/PanSV/GenePrediction/RNA/Stringtie/CPG11/CPG11_RNA_stringtie.gtf  --repeats /home/wuzhikun/Project/PanSV/Repeats/repeatMasker/CPG11/CPG11.genome.fasta.out.gff --CPU 16  > /home/wuzhikun/Project/BAssembly/pipeline/evidencemodeler/CPG11.evm.out 2>CPG11.evm.out.log  
        if not os.path.exists(params.outDir):
            os.makedirs(params.outDir)
        ### --repeats %s
        cmd = "source activate Anno &&  cd %s &&  EVidenceModeler --sample_id %s --segmentSize 100000 --overlapSize 10000  --genome %s  --weights %s --gene_predictions  %s --protein_alignments  %s --transcript_alignments  %s   --CPU %s  > %s 2>%s " % (params.outDir, wildcards.sample, input.genome, params.weight, input.augustus, input.miniprot, input.stringtie, threads, output.evidence, log)
        print(cmd)
        os.system(cmd)






################## rename gene ###################
rule gffRename:
    input:
        gff = IN_PATH + "/GenePrediction/evidencemodeler/{sample}/{sample}.EVM.gff3",
    output:
        gff = IN_PATH + "/GenePrediction/Gene/{sample}/{sample}.gene.gff3",
        gene = IN_PATH + "/GenePrediction/Gene/{sample}/{sample}.rename.txt",
    params:
        EvmGffRename = SRC_DIR + '/EvmGffRename.py',
    run:
        shell("python {params.EvmGffRename} --gff {input.gff} --out {output.gff} --out2 {output.gene}")




rule rename2:
    input:
        gff = IN_PATH + "/GenePrediction/Gene/{sample}/{sample}.gene.gff3",
        pep = IN_PATH + "/GenePrediction/Gene/{sample}/{sample}.pep.fasta",
        cds = IN_PATH + "/GenePrediction/Gene/{sample}/{sample}.cds.fasta",
    output:
        gff = IN_PATH + "/GenePrediction/Gene2/{sample}/{sample}.gene.gff3",
        pep = IN_PATH + "/GenePrediction/Gene2/{sample}/{sample}.pep.fasta",
        cds = IN_PATH + "/GenePrediction/Gene2/{sample}/{sample}.cds.fasta",
    params:
        GeneAddName = SCRIPT_DIR + "/GeneAddName.py",
    run:
        shell("python {params.GeneAddName} --input {input.gff} --out {output.gff} --format gff --name {wildcards.sample}")
        shell("python {params.GeneAddName} --input {input.pep} --out {output.pep} --format fasta --name {wildcards.sample}")
        shell("python {params.GeneAddName} --input {input.cds} --out {output.cds} --format fasta --name {wildcards.sample}")

##########################################################






'''



rule barrnap:
    input:
        assembly = IN_PATH + "/Scaffold/{sample}.genome.fasta",
    output:
        fa = IN_PATH + "/GenePrediction/barrnap/{sample}.barrnap_seq.fasta",
        predict = IN_PATH + "/GenePrediction/barrnap/{sample}.barrnap_predict.txt",
    threads:
        THREADS
    run:
        ### Env: Anno
        shell("barrnap --kingdom euk --threads {threads} --outseq {output.fa} {input.assembly} > {output.predict}")


rule rfamInfernalPan:
    input:
        assembly = IN_PATH + "/Scaffold/{sample}.genome.fasta",
    output:
        table = IN_PATH + "/GenePrediction/Rfam/{sample}_Rfam_table.txt",
        hit = IN_PATH + "/GenePrediction/Rfam/{sample}_Rfam_predict.txt",
    params:
        Rfam = "/home/wuzhikun/database/Rfam/Rfam.cm"
    threads:
        THREADS 
    log:
        IN_PATH + '/log/rfam_{sample}.log' 
    run:
        ### infernal
        # shell("cmpress {params.Rfam}")
        shell("cmscan --cpu {threads}  --tblout {output.table}  -o {output.hit}   {params.Rfam} {input.assembly} > {log} 2>&1")


rule tRNAscanPan:
    input:
        assembly = IN_PATH + "/Scaffold/{sample}.genome.fasta",
    output:
        table = IN_PATH + "/GenePrediction/tRNAscan/{sample}_tRNAscan_table.txt",
        stats = IN_PATH + "/GenePrediction/tRNAscan/{sample}_tRNAscan_stats.txt",
        structure = IN_PATH + "/GenePrediction/tRNAscan/{sample}_tRNAscan_structure.txt",
        bed = IN_PATH + "/GenePrediction/tRNAscan/{sample}_tRNAscan_bed.txt",
    threads:
        THREADS 
    log:
        IN_PATH + '/log/tRNAscan_{sample}.log' 
    run:
        ### -g ~/soft/miniconda3/pkgs/trnascan-se-2.0.6-pl526h516909a_0/lib/tRNAscan-SE/gcode/gcode.vertmito
        shell("tRNAscan-SE -qQ --detail --thread {threads} -o {output.table} -m {output.stats} -f {output.structure} -b {output.bed}   {input.assembly} > {log} 2>&1")

'''


rule RNAIndStats:
    input:
        tRNAscan = IN_PATH + "/GenePrediction/tRNAscan/{sample}_tRNAscan_bed.txt",
        barrnap = IN_PATH + "/GenePrediction/barrnap/{sample}.barrnap_predict.txt",
        Rfam = IN_PATH + "/GenePrediction/Rfam/{sample}_Rfam_table.txt",
    output:
        RNA = IN_PATH + "/GenePrediction/RNA/{sample}_RNA_stats.txt",
    run:
        shell("echo {wildcards.sample} > {output.RNA}") 
        shell("wc -l {input.tRNAscan} | cut -f 1 -d ' ' >> {output.RNA}") ## tRNA
        shell("grep -c rRNA {input.barrnap} >> {output.RNA}") ## rRNA
        shell("grep -c microRNA {input.Rfam} >> {output.RNA}") ## miRNA
        shell("grep -c 'Small nucleolar' {input.Rfam} >> {output.RNA}") ##snoRNA


rule MergeRNA:
    input:
        RNA = expand(IN_PATH + "/GenePrediction/RNA/{sample}_RNA_stats.txt", sample=SAMPLES),
    output:
        head = IN_PATH + "/GenePrediction/RNA/All_samples_RNA_header.xls", 
        RNA = IN_PATH + "/GenePrediction/RNA/All_samples_RNA_stats.xls",
    run:
        Files = " ".join(sorted(input.RNA))
        shell("echo Sample > {output.head}")
        shell("echo tRNA >> {output.head}")
        shell("echo rRNA >> {output.head}")
        shell("echo miRNA >> {output.head}")
        shell("echo snoRNA >> {output.head}")
        shell("paste {output.head} {Files} > {output.RNA}")















'''
rule hisatEachCPG01:
    input:
        ht = IN_PATH + "/GenePrediction/RNA/HisatIndex/{sample}_hisat.1.ht2",
        R1 = IN_PATH + '/clean/CPG01.RNA.R1.fq.gz',
        R2 = IN_PATH + '/clean/CPG01.RNA.R2.fq.gz',
    output:
        sam = temp(IN_PATH + "/GenePrediction/RNA/hisatCPG01/{sample}/ngs_RNA_hisat.sam"),
    params:
        prefix = IN_PATH + "/GenePrediction/RNA/HisatIndex/{sample}_hisat",
    threads:
        THREADS
    log:
        IN_PATH + "/log/hisatEach_{sample}.log", 
    run:
        shell("hisat2 --dta -q --sensitive --threads {threads} -x {params.prefix} -1 {input.R1} -2 {input.R2} -S {output.sam} > {log} 2>&1")


rule hisatBAMCPG01:
    input:
        sam = IN_PATH + "/GenePrediction/RNA/hisatCPG01/{sample}/ngs_RNA_hisat.sam",
    output:
        bam = IN_PATH + "/GenePrediction/RNA/hisatCPG01/{sample}/ngs_RNA_hisat.bam",
    threads:
        THREADS
    log:
        IN_PATH + "/log/hisatbam_{sample}.log", 
    run:
        shell("samtools view -@ {threads} -b {input.sam} | samtools sort -  -@ {threads} -o {output.bam} > {log} 2>&1")
        shell("samtools index {output.bam}")


rule BAMStatsCPG01:
    input:
        bam = IN_PATH + "/GenePrediction/RNA/hisatCPG01/{sample}/ngs_RNA_hisat.bam",
    output:
        stat = IN_PATH + '/GenePrediction/RNA/hisatCPG01/{sample}_ngs_RNA_stats.xls',
    threads:
        THREADS
    log:
        IN_PATH + "/log/BAMStats3_{sample}.log"
    run:
        shell("samtools flagstat --threads {threads}  {input.bam} > {output.stat} 2>{log}")


'''



rule gtf:
    input:
        gff = IN_PATH + "/GenePrediction/evidencemodeler/{sample}/{sample}.EVM.gff3",
    output:
        gtf = IN_PATH + "/GenePrediction/evidencemodeler/{sample}/{sample}.EVM.gtf",
    run:
        shell("gffread  {input.gff} -T -o {output.gtf}")


rule STARIndex:
    input:
        scaffold = IN_PATH + "/Scaffold/{sample}.genome.fasta",
        gtf = IN_PATH + "/GenePrediction/evidencemodeler/{sample}/{sample}.EVM.gtf",
    output:
        SA = IN_PATH + "/GenePrediction/STARIndex/{sample}/SA",
    params:
        outDir = IN_PATH + "/GenePrediction/STARIndex/{sample}",
        tempDir = IN_PATH + "/GenePrediction/TempDir/{sample}",
    threads:
        THREADS
    log:
        IN_PATH + "/log/STARIndex_{sample}.log"
    run:
        if not os.path.exists(params.outDir):
            os.makedirs(params.outDir)
        shell("STAR --runThreadN {threads} --genomeSAindexNbases 13 --runMode genomeGenerate --genomeDir {params.outDir} --genomeFastaFiles {input.scaffold} --sjdbGTFfile {input.gtf} --sjdbOverhang 149 --outTmpDir {params.tempDir}  > {log} 2>&1")





rule STAR:
    input:
        R1 = IN_PATH + '/clean/{sample}.RNA.R1.fq.gz',
        R2 = IN_PATH + '/clean/{sample}.RNA.R2.fq.gz',
        SA = IN_PATH + "/GenePrediction/STARIndex/{sample}/SA",
    output:
        bam = IN_PATH + '/RNA/NGS/mapping/{sample}Aligned.sortedByCoord.out.bam',
    params:
        IndexDir = IN_PATH + "/GenePrediction/STARIndex/{sample}",
        outDir = IN_PATH + '/RNA/NGS/mapping/{sample}',
    threads:
        THREADS
    log:
        IN_PATH + "/log/STAR/{sample}.log"
    run:
        if not os.path.exists(params.outDir):
            os.makedirs(params.outDir)
        shell(' STAR --runThreadN {threads} --genomeDir {params.IndexDir}  --readFilesIn {input.R1}  {input.R2}  --readFilesCommand gunzip -c --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {params.outDir}   > {log} 2>&1')
        shell("samtools index {output.bam}")




rule StarIndex:
    input:
        bam = IN_PATH + '/RNA/NGS/mapping/{sample}Aligned.sortedByCoord.out.bam',
    output:
        bai = IN_PATH + '/RNA/NGS/mapping/{sample}Aligned.sortedByCoord.out.bam.bai',
    run:
        shell("samtools index {input.bam}")        



rule BAMStats2:
    input:
        bam = IN_PATH + '/RNA/NGS/mapping/{sample}Aligned.sortedByCoord.out.bam',
        bai = IN_PATH + '/RNA/NGS/mapping/{sample}Aligned.sortedByCoord.out.bam.bai',
    output:
        stat = IN_PATH + '/RNA/NGS/mapping/{sample}.bam_stats.xls',
    threads:
        THREADS
    log:
        IN_PATH + "/log/BAMStats2_{sample}.log"
    run:
        shell("samtools stats --threads {threads} {input.bam} > {output.stat} 2>{log}")

rule BAMStats3:
    input:
        bam = IN_PATH + "/RNA/NGS/mapping/{sample}Aligned.sortedByCoord.out.bam",
        bai = IN_PATH + '/RNA/NGS/mapping/{sample}Aligned.sortedByCoord.out.bam.bai',
    output:
        stat = IN_PATH + '/RNA/NGS/mapping/{sample}_bam_flag_stats.xls',
    threads:
        THREADS
    log:
        IN_PATH + "/log/BAMStats3_{sample}.log"
    run:
        shell("samtools flagstat --threads {threads} -O tsv {input.bam} > {output.stat} 2>{log}")





rule featureCount:
    input:
        IN_PATH + "/RNA/NGS/mapping/{sample}Aligned.sortedByCoord.out.bam",
    output:
        IN_PATH + '/RNA/NGS/Counts/{sample}_Reads_count.xls'
    threads:
        THREADS
    log:
        IN_PATH + '/log/{sample}_featurecount.log'
    params:
        feature_type = "exon", # "gene",  #config['feature_type'],
        feature_attribute = "Parent", # "ID", # config['feature_attribute'],
        GFF = IN_PATH + "/GenePrediction/Gene/{sample}/{sample}.gene.gff3",
    run:
        shell('featureCounts -T {threads} -d 50 -D 1000 -C -s 0 -t {params.feature_type} -g {params.feature_attribute} --primary -O -M --fraction -a {params.GFF} -o {output} {input} > {log} 2>&1' )



rule merge_featurecount:
    input:
        expand(IN_PATH + '/RNA/NGS/Counts/{sample}_Reads_count.xls', sample=SAMPLES)
    output:
        IN_PATH + '/RNA/NGS/Counts/All_samples_reads_counts.xls',
    params:
        merge_featurecount = SCRIPT_DIR + '/merge_featurecount.py'
    log:
        IN_PATH + '/log/All_reads_counts.log',
    run:
        INPUT = ','.join(input)
        shell('python {params.merge_featurecount} -i {INPUT} -o {output}')


