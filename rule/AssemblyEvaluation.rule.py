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

######################## scaffold stats ##############
rule scaffoldStats:
    input:
        assembly = IN_PATH + "/Assembly/Scaffold/{sample}.scaffold.fasta",
    output:
        stat = IN_PATH + "/Assembly/Scaffold/Stats/{sample}.scaffold.stats.txt",
    run:
        shell("seqkit stats -aT {input.assembly} > {output.stat}")

rule anchorRatio:
    input:
        assembly = IN_PATH + "/Assembly/Scaffold/{sample}.scaffold.fasta",
    output:
        anchor = IN_PATH + "/Assembly/Scaffold/Stats/{sample}.anchor.ratio.txt",
    params:
        AnchorRatio = SCRIPT_DIR + "/AnchorRatio.py",
    run:
        shell("python {params.AnchorRatio} --fasta {input.assembly} --out {output.anchor}")

"""
rule AnchorRatioSum:
    input:
        anchor = expand(IN_PATH + "/Summary/anchor/{species}/{sample}.anchor_ratio.txt", species=SPECIES, sample=SAMPLES),
    output:
        summary = IN_PATH + "/Summary/anchor/{species}/scaffold.anchor_ratio.summary.xls",
    run:
        Files = sorted(input.anchor)
        fileNum = len(Files)
        for i in range(fileNum):
            f = Files[i]
            if i == 0:
                cmd = "cat %s > %s" % (f, output.summary)
            else:
                cmd = "sed '1d' %s >> %s" % (f, output.summary)
            os.system(cmd)
"""
#####################################################


###################### Coverage #########################

rule ONTCov:
    input:
        bam = IN_PATH + "/Assembly/RagTagAlign/{sample}_ONTAlign.bam",
    output:
        bw = IN_PATH + "/Evaluation/Coverage/{sample}/{sample}_ONT_WGS.bw",
    log:
        IN_PATH + "/log/ONTCov_{sample}.log"
    run:
        shell("bamCoverage -b {input.bam} -o {output.bw} > {log} 2>&1")



rule ONTDepth:
    input:
        bam = IN_PATH + "/Assembly/RagTagAlign/{sample}_ONTAlign.bam",
    output:
        depth = IN_PATH +  "/Evaluation/Coverage/{sample}/{sample}_ONT.per-base.bed.gz",
    threads:
        THREADS
    params:
        outPrefix = IN_PATH + '/Evaluation/Coverage/{sample}/{sample}_ONT',
    log:
        IN_PATH + "/log/ONTDepth_{sample}.log"
    run:
        # shell('mosdepth --threads {threads} --fast-mode --flag 256  {params.outPrefix} {input.bam} > {log} 2>&1')
        cmd = "source activate NanoSV && mosdepth --threads %s --fast-mode --flag 256  %s %s > %s 2>&1" % (threads, params.outPrefix, input.bam, log)
        print(cmd)
        os.system(cmd)
##################################################################


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



rule MergeQV:
    input:
        nextDenovo = expand(IN_PATH + "/QualityValue/NextDenovo/{sample}/{sample}.qv", sample=SAMPLES),
        ONT = expand(IN_PATH + "/QualityValue/ONTPolish/{sample}/{sample}.qv", sample=SAMPLES),
        NGS = expand(IN_PATH + "/QualityValue/NGSPolish/{sample}/{sample}.qv", sample=SAMPLES),
    output:
        nextDenovo = IN_PATH + "/QualityValue/NextDenovo_samples.qv.xls",
        ONT = IN_PATH + "/QualityValue/ONTPolish_samples.qv.xls",
        NGS = IN_PATH + "/QualityValue/NGSPolish_samples.qv.xls",
    run:
        NextDenovos = ",".join(sorted(input.nextDenovo))
        ONTS = ",".join(sorted(input.ONT))
        NGSS = ",".join(sorted(input.NGS))
        merge_qv(NextDenovos, output.nextDenovo)
        merge_qv(ONTS, output.ONT)
        merge_qv(NGSS, output.NGS)


def merge_qv(in_file, out_file):
    out_h = open(out_file, "w")
    files = in_file.split(",")
    for i in files:
        n = i.split("/")[-2]
        in_h = open(i, "r")
        line = in_h.readline().strip()
        out_h.write("%s\t%s\n" % (n, line))
        in_h.close()
    out_h.close()

########################################


################## centromeres ##########
rule SatSeqAlign:
    input:
        fa =  "/home/wuzhikun/Project/Vigna/BACKUP_230424/Vigna_SRF/Vigna_satellite.fa",
        assembly = IN_PATH + "/Assembly/Scaffold/{sample}.scaffold.fasta",
    output:
        paf = IN_PATH + "/Centromere/SRF/{sample}_srf_aln.paf",
    threads:
        THREADS
    params:
        srfutils = "/home/wuzhikun/software/srf/srfutils.js",
    log:
        IN_PATH + "/log/SatSeqAlign_{sample}.log"
    run:
        shell("minimap2 -c -t {threads}  -N1000000 -f1000 -r100,100 <({params.srfutils} enlong {input.fa}) {input.assembly} > {output.paf}")


rule paf2bed:
    input:
        paf =  IN_PATH + "/Centromere/SRF/{sample}_srf_aln.paf",
    output:
        bed = IN_PATH + "/Centromere/SRF/{sample}_srf_aln.bed",
        abun = IN_PATH + "/Centromere/SRF/{sample}_srf_aln.len",
    params:
        srfutils = "/home/wuzhikun/software/srf/srfutils.js",
    threads:
        THREADS
    run:
        shell("{params.srfutils} paf2bed {input.paf} > {output.bed}")
        shell("{params.srfutils} bed2abun {output.bed} > {output.abun}")



rule satellite:
    input:
        bed = IN_PATH + "/Centromere/SRF/{sample}_srf_aln.bed",
    output:
        satellite = IN_PATH + "/Centromere/SRF/{sample}_srf_aln.satellite.txt",
    threads:
        THREADS
    params:
        CenRegionRepeat = SCRIPT_DIR + "/CenRegionRepeat.py",
    run:
        shell("python {params.CenRegionRepeat} --record {input.bed} --satellite {output.satellite}")



rule cluster:
    input:
        bed = IN_PATH + "/Centromere/SRF/{sample}_srf_aln.bed",
    output:
        bed = temp(IN_PATH + "/Centromere/SRF/{sample}_srf_aln_temp.bed"),
        cluster = IN_PATH + "/Centromere/SRF/{sample}_srf_aln_cluster.bed",
    threads:
        THREADS
    run:
        shell("cut -f 1-3 {input.bed} > {output.bed}")
        shell("bedtools cluster -i {output.bed} -d 100000 > {output.cluster}")


rule targetRegion:
    input:
        cluster = IN_PATH + "/Centromere/SRF/{sample}_srf_aln_cluster.bed",
    output:
        target = IN_PATH + "/Centromere/SRF/{sample}_srf_aln_cluster_target.txt",
    params:
        ClusterLongest = SCRIPT_DIR + "/ClusterLongest.py",
    threads:
        THREADS
    run:
        shell("python {params.ClusterLongest} --cluster {input.cluster} --out {output.target}")


#########################################


################# centromere2 ##########
rule SatSeqAlign2:
    input:
        fa =  "/home/wuzhikun/Project/PanVigna/BACKUP/centromere/cent_satellite_seq.fasta",
        assembly = IN_PATH + "/Assembly/Scaffold/{sample}.scaffold.fasta",
    output:
        paf = IN_PATH + "/Centromere/CEN/{sample}_satellite_aln.paf",
    threads:
        THREADS
    params:
        srfutils = "/home/wuzhikun/software/srf/srfutils.js",
    log:
        IN_PATH + "/log/SatSeqAlign2_{sample}.log"
    run:
        shell("minimap2 -c -t {threads}  -N1000000 -f1000 -r100,100 <({params.srfutils} enlong {input.fa}) {input.assembly} > {output.paf}")


rule paf2bed2:
    input:
        paf =  IN_PATH + "/Centromere/CEN/{sample}_satellite_aln.paf",
    output:
        bed = IN_PATH + "/Centromere/CEN/{sample}_satellite_aln.bed",
        abun = IN_PATH + "/Centromere/CEN/{sample}_satellite_aln.len",
    params:
        srfutils = "/home/wuzhikun/software/srf/srfutils.js",
    threads:
        THREADS
    run:
        shell("{params.srfutils} paf2bed {input.paf} > {output.bed}")
        shell("{params.srfutils} bed2abun {output.bed} > {output.abun}")


rule cluster2:
    input:
        bed = IN_PATH + "/Centromere/CEN/{sample}_satellite_aln.bed",
    output:
        bed = temp(IN_PATH + "/Centromere/CEN/{sample}_srf_aln_temp.bed"),
        cluster = IN_PATH + "/Centromere/CEN/{sample}_srf_aln_cluster.bed",
    threads:
        THREADS
    run:
        shell("cut -f 1-3 {input.bed} > {output.bed}")
        shell("bedtools cluster -i {output.bed} -d 100000 > {output.cluster}")

########################################



################## BUSCO ############

rule BUSCO:
    input:
        assembly = IN_PATH + "/Assembly/RagTag/{sample}/ragtag.scaffold.fasta",
    output:
        summary = IN_PATH + "/Evaluation/BUSCO/ragtag/{sample}/run_embryophyta_odb10/short_summary.txt",
    threads:
        THREADS 
    params:
        buscoDB = "/home/wuzhikun/database/BUSCO/version5/embryophyta_odb10",
        outDir = IN_PATH + "/Evaluation/BUSCO/ragtag",
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


rule buscoSummary:
    input:
        stats = expand(IN_PATH + "/Evaluation/BUSCO/ragtag/{sample}/run_embryophyta_odb10/short_summary.txt", sample=SAMPLES),
    output:
        summary = IN_PATH + "/Evaluation/BUSCO/ragtag/All_busco_summary.xls",
    params:
        buscoSummary = SCRIPT_DIR + "/buscoSummary.py"
    run:
        Files = ",".join(input.stats)
        shell("python {params.buscoSummary} --file {Files} --out {output.summary}")



#########################################



########### inspector ###############
# rule inspector:
#     input:
#         assembly = IN_PATH + "/Assembly/Scaffold/{sample}.scaffold.fasta",
#         read = IN_PATH + "/clean/{sample}_ONT.fastq.gz",
#     output:
#         summary = IN_PATH + "/Evaluation/Inspector/{sample}/summary_statistics",
#     params:
#         inspector = "/home/wuzhikun/software/Inspector/inspector.py",
#         outDir = IN_PATH + "/Evaluation/Inspector/{sample}",
#     threads:
#         THREADS
#     log:
#         IN_PATH + '/log/inspector_{sample}.log' 
#     run:
#         shell("python {params.inspector} --contig {input.assembly} --read {input.read} --datatype nanopore  --thread {threads} --outpath {params.outDir} > {log} 2>&1")
#####################################




# ######################### RepeatModeler ####################
# rule BuildDatabase:
#     input:
#         #assembly = rules.RagTag.output.scaffold,
#         assembly = IN_PATH + "/Assembly/Polish/{sample}/Racon/{sample}_ONT_polish_racon3.fasta",
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
#         # shell("cd {params.outDir} &&  BuildDatabase -name {params.outPrefix} -engine ncbi  {input.assembly} > {log} 2>&1")
#         cmd = "source activate TE2 && cd %s &&  BuildDatabase -name %s -engine ncbi  %s > %s 2>&1" % (params.outDir, params.outPrefix, input.assembly, log)
#         print(cmd)
#         os.system(cmd)



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
#         # shell("cd {params.outDir} &&  RepeatModeler -pa {threads} -database {params.outPrefix}  -engine ncbi > {log} 2>&1")
#         cmd = "source activate TE2 && cd %s &&  RepeatModeler -pa %s -database %s  -engine ncbi > %s 2>&1" % (params.outDir, threads, params.outPrefix, log)
#         print(cmd)
#         os.system(cmd)

# ######################################################


# #################### repeatMasker ##########
rule repeatMasker:
    input:
        #assembly = rules.RagTag.output.scaffold,
        assembly = IN_PATH + "/Assembly/Scaffold/{sample}.scaffold.fasta",
        lib = IN_PATH + "/BACKUP/Repeat/Vigna_unguiculata_contig-families.fa",
    output:
        mask = IN_PATH + "/Repeat/RepeatMasker/{sample}/{sample}.scaffold.fasta.masked",
        tbl = IN_PATH + "/Repeat/RepeatMasker/{sample}/{sample}.scaffold.fasta.tbl",
    threads:
        THREADS
    params:
        outPrefix = IN_PATH + "/Repeat/RepeatMasker/{sample}",
    log:
        IN_PATH + "/log/repeatMasker_{sample}.log",
    run:
        # -libdir <string>
        #     Path to the RepeatMasker libraries directory.
        # shell("RepeatMasker -pa {threads} -e ncbi  -gff -poly -lib {input.lib}  -dir {params.outPrefix}  {input.assembly} > {log} 2>&1")
        cmd = "source activate TE2 && RepeatMasker -xsmall -nolow -norna -html -gff -pa %s -e ncbi  -poly -lib %s  -dir %s  %s > %s 2>&1" % (threads, input.lib, params.outPrefix, input.assembly, log)
        print(cmd)
        os.system(cmd)


###########################################################
















###################### ONT ##########################
rule ONTAlignSelf:
    input:
        fastq = IN_PATH + "/clean/{sample}_ONT.fastq.gz",
        assembly = IN_PATH + "/Scaffold/{sample}.genome.fasta"
    output:
        sam = temp(IN_PATH + "/Evaluation/mapping/ONT/{sample}_ONT.sam"),
        bam = IN_PATH + "/Evaluation/mapping/ONT/{sample}_ONT.bam",
    threads:
        THREADS
    log:
        IN_PATH + "/log/ONTAlignselt_{sample}.log"
    run:
        shell("minimap2 --MD -a -x map-ont -t {threads} {input.assembly} {input.fastq} > {output.sam} 2>{log}")
        shell("samtools view -@ {threads} -b {output.sam} | samtools sort -  -@ {threads} -o {output.bam} > {log} 2>&1")
        shell("samtools index {output.bam}")



rule ONTCovSelf:
    input:
        bam = IN_PATH + "/Evaluation/mapping/ONT/{sample}_ONT.bam",
    output:
        bw = IN_PATH + "/Evaluation/mapping/ONT/{sample}_ONT_WGS.bw",
    log:
        IN_PATH + "/log/ONTCov_{sample}.log"
    run:
        shell("/home/wuzhikun/anaconda3/envs/RNA2/bin/bamCoverage -b {input.bam} -o {output.bw} > {log} 2>&1")



rule ONTDepthSelf:
    input:
        bam = IN_PATH + "/Evaluation/mapping/ONT/{sample}_ONT.bam",
    output:
        depth = IN_PATH +  "/Evaluation/mapping/ONT/{sample}/{sample}_ONT.per-base.bed.gz",
    threads:
        THREADS
    params:
        outPrefix = IN_PATH + '/Evaluation/mapping/ONT/{sample}/{sample}_ONT',
    log:
        IN_PATH + "/log/ONTDepthSelf_{sample}.log"
    run:
        # shell('mosdepth --threads {threads} --fast-mode --flag 256  {params.outPrefix} {input.bam} > {log} 2>&1')
        cmd = "source activate NanoSV && mosdepth --threads %s --fast-mode --flag 256  %s %s > %s 2>&1" % (threads, params.outPrefix, input.bam, log)
        print(cmd)
        os.system(cmd)



rule SnifflesSelf:
    input:
        bam = IN_PATH + "/Evaluation/mapping/ONT/{sample}_ONT.bam",
    output:
        vcf = IN_PATH + "/Evaluation/mapping/ONT/{sample}.sniffles.vcf",
    threads:
        THREADS 
    log:
        IN_PATH + "/log/SnifflesSelf_{sample}.log"
    run:
        shell("sniffles --mapped_reads {input.bam} --vcf {output.vcf} --threads {threads}  --min_support 3 --min_length 50 --minmapping_qual 20 --num_reads_report -1 --min_seq_size 500  --genotype --report_BND --report-seq  >{log} 2>&1")



rule cuteSVSelf:
    input:
        bam = IN_PATH + "/Evaluation/mapping/ONT/{sample}_ONT.bam",
    output:
        vcf = IN_PATH + "/Evaluation/mapping/ONT/{sample}.cutesv.vcf",
    threads:
        THREADS 
    params:
        RefGenome = "/Scaffold/{sample}.genome.fasta",
        tempDir = IN_PATH + "/Evaluation/mapping/ONT/cuteSVTemp/{sample}",
    log:
        IN_PATH + "/log/cuteSVSelf_{sample}.log"
    run:
        if not os.path.exists(params.tempDir):
            os.makedirs(params.tempDir)
        cmd = "source activate nanovar && cuteSV --max_cluster_bias_INS 100   --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --threads %s --sample  %s --report_readid --min_support 5  --min_size 50 --min_siglength 50  --max_size -1 --genotype %s  %s %s %s > %s  2>&1" % (threads, wildcards.sample, input.bam, params.RefGenome, output.vcf, params.tempDir, log)
        print(cmd)
        os.system(cmd)





############################################




###################### NGS ###################

rule SelfBWAIndex:
    input:
        ref = IN_PATH + "/Scaffold/{sample}.genome.fasta",
    output:
        bwt = IN_PATH + "/Scaffold/{sample}.genome.fasta.bwt",
    run:
        if not os.path.exists(output.bwt):
            cmd = 'bwa index %s' % input.ref
            os.system(cmd)


rule SelfBWAAlign:
    input:
        ref = IN_PATH + "/Scaffold/{sample}.genome.fasta",
        bwt = IN_PATH + "/Scaffold/{sample}.genome.fasta.bwt",
        R1 = IN_PATH + '/clean/{sample}_NGS.R1.fastq.gz',
        R2 = IN_PATH + '/clean/{sample}_NGS.R2.fastq.gz',
    output:
        sam = temp(IN_PATH + '/Evaluation/mapping/NGS/{sample}.sam'),
    threads:
        THREADS
        # rg = "@RG\tID:01\tPL:ILLUMINA\tPU:{sample}\tSM:{sample}",
    log:
        IN_PATH + "/log/{sample}.SelfBWAAlign.log",
    run:
        shell('bwa mem -M -t {threads} {input.ref} {input.R1} {input.R2} > {output.sam} 2>{log}')



rule SelfSAM2BAM:
    input:
        sam = IN_PATH + '/Evaluation/mapping/NGS/{sample}.sam',
    output:
        tempbam = temp(IN_PATH + '/Evaluation/mapping/NGS/{sample}_temp.bam'),
        bam = IN_PATH + '/Evaluation/mapping/NGS/{sample}_align.bam',
    threads:
        THREADS
    log:
        IN_PATH + "/log/{sample}.SelfSAM2BAM.log",
    run:
        shell('sambamba view --nthreads {threads} --sam-input --format bam  -o {output.tempbam} {input.sam} >{log} 2>&1' )
        shell('samtools sort --threads {threads} -o {output.bam} {output.tempbam} 2>>{log}')
        shell("samtools index {output.bam}")




rule SelfBAMStats:
    input:
        bam = IN_PATH + '/Evaluation/mapping/NGS/{sample}_align.bam',
    output:
        stat = IN_PATH + '/Evaluation/mapping/NGS/{sample}_bam_stats.xls',
    threads:
        THREADS
    log:
        IN_PATH + "/log/SelfBAMStats_{sample}.log"
    run:
        shell("samtools flagstat --threads {threads} {input.bam} > {output.stat} 2>{log}")



rule NGSCovSelf:
    input:
        bam = IN_PATH + "/Evaluation/mapping/NGS/{sample}_align.bam",
    output:
        bw = IN_PATH + "/Evaluation/mapping/NGS/{sample}_NGS_WGS.bw",
    log:
        IN_PATH + "/log/NGSCov_{sample}.log"
    run:
        shell("/home/wuzhikun/anaconda3/envs/RNA2/bin/bamCoverage -b {input.bam} -o {output.bw} > {log} 2>&1")



rule NGSDepthSelf:
    input:
        bam = IN_PATH + "/Evaluation/mapping/NGS/{sample}_align.bam",
    output:
        depth = IN_PATH +  "/Evaluation/mapping/NGS/{sample}/{sample}_NGS.per-base.bed.gz",
    threads:
        THREADS
    params:
        outPrefix = IN_PATH + '/Evaluation/mapping/NGS/{sample}/{sample}_NGS',
    log:
        IN_PATH + "/log/NGSDepthSelf_{sample}.log"
    run:
        # shell('mosdepth --threads {threads} --fast-mode --flag 256  {params.outPrefix} {input.bam} > {log} 2>&1')
        cmd = "source activate NanoSV && mosdepth --threads %s --fast-mode --flag 256  %s %s > %s 2>&1" % (threads, params.outPrefix, input.bam, log)
        print(cmd)
        os.system(cmd)


###################################################################




rule AssemblyMashmap:
    input:
        assembly = IN_PATH + "/Scaffold/{sample}.genome.fasta",
    output:
        mashmap = IN_PATH + "/Evaluation/mashmap/{sample}/mashmap.out.txt",
    params:
        RefGenome = "/home/wuzhikun/Project/PanSV/BACKUP/Assembly/Fengchan6/Vigna_unguiculata_assembly.fasta",
        mashmap = "/home/wuzhikun/software/mashmap-Linux64-v2.0/mashmap",
        generateDotPlot = "/home/wuzhikun/software/mashmap-Linux64-v2.0/MashMap/scripts/generateDotPlot",
        outPrefix = IN_PATH + "/Evaluation/mashmap/{sample}",
    threads:
        THREADS
    log:
        IN_PATH + "/log/AssemblyMashmap_{sample}.log"
    run:
        shell("{params.mashmap} -r {params.RefGenome} -q {input.assembly} --perc_identity 95 --threads {threads} --segLength 50000 --filter_mode one-to-one --output {output.mashmap} > {log} 2>&1")
        shell("cd {params.outPrefix} && {params.generateDotPlot} png  large  {output.mashmap}")



