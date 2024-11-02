



# ############### Evaluation #############

# rule ONTAlign:
#     input:
#         fastq = IN_PATH + "/clean/{sample}_ONT.fastq.gz",
#         assembly = rules.RagTag.output.scaffold,
#     output:
#         sam = temp(IN_PATH + "/SVCall/mapping/{sample}_ONT.sam"),
#         bam = IN_PATH + "/SVCall/mapping/{sample}_ONT.bam",
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/ONTAlign_{sample}.log"
#     run:
#         shell("minimap2 --MD -a -x map-ont -t {threads} {input.assembly} {input.fastq} > {output.sam} 2>{log}")
#         shell("samtools view -@ {threads} -b {output.sam} | samtools sort -  -@ {threads} -o {output.bam} > {log} 2>&1")
#         shell("samtools index {output.bam}")



# # rule ONTCov:
# #     input:
# #         bam = rules.ONTAlign.output.bam,
# #     output:
# #         bw = IN_PATH + "/Evaluation/Coverage/{sample}/{sample}_ONT_WGS.bw",
# #     log:
# #         IN_PATH + "/log/ONTCov_{sample}.log"
# #     run:
# #         shell("bamCoverage -b {input.bam} -o {output.bw} > {log} 2>&1")


# # rule ONTDepth:
# #     input:
# #         bam = rules.ONTAlign.output.bam,
# #     output:
# #         depth = IN_PATH +  "/Evaluation/Coverage/{sample}/{sample}_ONT.per-base.bed.gz",
# #     threads:
# #         THREADS
# #     params:
# #         outPrefix = IN_PATH + '/Evaluation/Coverage/{sample}/{sample}_ONT',
# #     log:
# #         IN_PATH + "/log/ONTDepth_{sample}.log"
# #     run:
# #         # shell('mosdepth --threads {threads} --fast-mode --flag 256  {params.outPrefix} {input.bam} > {log} 2>&1')
# #         cmd = "source activate NanoSV && mosdepth --threads %s --fast-mode --flag 256  %s %s > %s 2>&1" % (threads, params.outPrefix, input.bam, log)
# #         print(cmd)
# #         os.system(cmd)

# ########################################################


'''

# ######################### SV ###########################
rule ONTAlignRef:
    input:
        fastq = IN_PATH + "/clean/{sample}_ONT.fastq.gz",
        assembly = "/home/wuzhikun/Project/Vigna/Final/Final/Vigna_unguiculata_assembly.fasta",
    output:
        sam = temp(IN_PATH + "/SVCall/mapping/{sample}_ONT.sam"),
        bam = IN_PATH + "/SVCall/mapping/{sample}_ONT.bam",
    threads:
        THREADS
    log:
        IN_PATH + "/log/ONTAlign_{sample}.log"
    run:
        shell("minimap2 --MD -a -x map-ont -t {threads} {input.assembly} {input.fastq} > {output.sam} 2>{log}")
        shell("samtools view -@ {threads} -b {output.sam} | samtools sort -  -@ {threads} -o {output.bam} > {log} 2>&1")
        shell("samtools index {output.bam}")
'''
        

rule BAMStats2:
    input:
        bam = IN_PATH + "/SVCall/mapping/{sample}_ONT.bam",
    output:
        stat = IN_PATH + '/SVCall/mapping/{sample}_ONT_bam_stats.xls',
    threads:
        THREADS
    log:
        IN_PATH + "/log/BAMStats2_{sample}.log"
    run:
        shell("samtools stats --threads {threads} {input.bam} > {output.stat} 2>{log}")

rule BAMStats3:
    input:
        bam = IN_PATH + "/SVCall/mapping/{sample}_ONT.bam",
    output:
        stat = IN_PATH + '/SVCall/mapping/{sample}_ONT_flag_stats.xls',
    threads:
        THREADS
    log:
        IN_PATH + "/log/BAMStats3_{sample}.log"
    run:
        shell("samtools flagstat --threads {threads} -O tsv {input.bam} > {output.stat} 2>{log}")



rule Sniffles1:
    input:
        bam = IN_PATH + "/SVCall/mapping/{sample}_ONT.bam",
    output:
        vcf = IN_PATH + "/SVCall/SV/{sample}.sniffles.vcf",
    threads:
        THREADS 
    log:
        IN_PATH + "/log/Sniffles1_{sample}.log"
    run:
        shell("sniffles --mapped_reads {input.bam} --vcf {output.vcf} --threads {threads}  --min_support 5 --min_length 50 --minmapping_qual 20 --num_reads_report -1 --min_seq_size 500  --genotype --report_BND --report-seq  >{log} 2>&1")


rule NanoSV1:
    input:
        bam = IN_PATH + "/SVCall/mapping/{sample}_ONT.bam",
    output:
        vcf = IN_PATH + "/SVCall/SV/{sample}.NanoSV.vcf",
    threads:
        THREADS 
    params:
        nanosvConfig = "/home/wuzhikun/anaconda3/envs/NanoSV/lib/python3.6/site-packages/nanosv/config.ini"
    log:
        IN_PATH + "/log/NanoSV_{sample}.log"
    run:
        ### depth_support = False
        cmd = "source activate NanoSV && NanoSV --threads %s -c %s -o %s %s > %s 2>&1" % (threads, params.nanosvConfig, output.vcf, input.bam, log)
        print(cmd)
        os.system(cmd)



rule cuteSV1:
    input:
        bam = IN_PATH + "/SVCall/mapping/{sample}_ONT.bam",
    output:
        vcf = IN_PATH + "/SVCall/SV/{sample}.cutesv.vcf",
    threads:
        THREADS 
    params:
        RefGenome = "/home/wuzhikun/Project/Vigna/Final/Final/Vigna_unguiculata_assembly.fasta",
        tempDir = IN_PATH + "/SVCall/SV/cuteSVTemp/{sample}",
    log:
        IN_PATH + "/log/cuteSV_{sample}.log"
    run:
	    # For ONT data:
		# --max_cluster_bias_INS		100
		# --diff_ratio_merging_INS	0.3
		# --max_cluster_bias_DEL	100
		# --diff_ratio_merging_DEL	0.3
        if not os.path.exists(params.tempDir):
            os.makedirs(params.tempDir)
        cmd = "source activate nanovar && cuteSV --max_cluster_bias_INS 100   --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --threads %s --sample  %s --report_readid --min_support 5  --min_size 50 --min_siglength 50  --max_size -1 --genotype %s  %s %s %s > %s  2>&1" % (threads, wildcards.sample, input.bam, params.RefGenome, output.vcf, params.tempDir, log)
        print(cmd)
        os.system(cmd)



rule nanovar:
    input:
        bam = IN_PATH + "/SVCall/mapping/{sample}_ONT.bam",
    output:
        total = IN_PATH + "/SVCall/nanovar/{sample}/{sample}.nanovar.total.vcf",
        vcf = IN_PATH + "/SVCall/nanovar/{sample}/{sample}.nanovar.pass.vcf",
    threads:
        THREADS
    params:
        RefGenome = "/home/wuzhikun/Project/Vigna/Final/Final/Vigna_unguiculata_assembly.fasta",
        outdir = IN_PATH + "/SVCall/nanovar/{sample}",
    log:
        IN_PATH + "/log/nanovar_{sample}.log"
    run:
        ### nanovar: /home/litong/anaconda3/envs/nanovar/bin/nanovar
        ### --data_type ont, pacbio-clr, pacbio-ccs
        # shell("nanovar -r {params.RefGenome} -l {input.fastq} -t {threads} -o {params.outdir} >{log} 2>&1")
        # cmd = "source activate NanoVar && nanovar -t %s --data_type ont --mincov 2 --minlen 50 %s %s %s > %s 2>&1" % (threads, input.bam, params.RefGenome, params.outdir, log)
        cmd = "source activate nanovar && nanovar -t %s --data_type ont --mincov 5 --minlen 50 %s %s %s > %s 2>&1" % (threads, input.bam, params.RefGenome, params.outdir, log)
        print(cmd)
        os.system(cmd)



rule SVision:
    input:
        bam = IN_PATH + "/SVCall/mapping/{sample}_ONT.bam",
    output:
        vcf = IN_PATH + "/SVCall/SVision/{sample}/{sample}.svision.s3.graph.vcf",
    log:
        IN_PATH + "/log/SVision_{sample}.log", 
    params:
        model = "/home/zhengjingjing/software/SVision/svision_model/svision-cnn-model.ckpt",
        RefGenome = "/home/wuzhikun/Project/Vigna/Final/Final/Vigna_unguiculata_assembly.fasta",
        outDir = IN_PATH + "/SVCall/SVision/{sample}",
    threads:
        THREADS
    run:
        # shell("source activate /home/zhengjingjing/anaconda3/envs/svisionenv && SVision  -o {params.outDir} -b {input.bam} -m {params.model} -g {params.RefGenome} -n {wildcards.sample} -s 3 --graph --qname -t {threads} > {log} 2>&1")
        cmd = "source activate /home/zhengjingjing/anaconda3/envs/svisionenv && SVision  -o %s -b %s -m %s -g %s -n %s -s 3 --graph --qname -t %s > %s 2>&1" % (params.outDir, input.bam, params.model, params.RefGenome, wildcards.sample, threads, log)
        print(cmd)
        os.system(cmd)



# ##########################################################


rule SVIM:
    input:
        bam = IN_PATH + "/SVCall/mapping/{sample}_ONT.bam",
    output:
        vcf = IN_PATH + "/SVCall/SVIM/{sample}/variants.vcf",
    params:
        RefGenome =  "/home/wuzhikun/Project/Vigna/Final/Final/Vigna_unguiculata_assembly.fasta",
        min_length = 50,
        minmapping_qual = 20,
        outdir = IN_PATH + "/SVCall/SVIM/{sample}",
    threads:
        THREADS
    log:
        IN_PATH + "/log/SVIM_{sample}.log"
    run:
        # shell("svim alignment --min_sv_size {params.min_length} --min_mapq {params.minmapping_qual} --sample {wildcards.sample} --insertion_sequences  {params.outdir} {input.bam} {params.RefGenome} >{log} 2>&1")
        cmd = "source activate NanoSV && svim alignment --min_sv_size 50 --min_mapq 20 --sample %s --insertion_sequences  %s %s  %s > %s 2>&1" % (wildcards.sample, params.outdir, input.bam, params.RefGenome, log)
        print(cmd)
        os.system(cmd)


rule SVIMFilt:
    input:
        vcf = IN_PATH + "/SVCall/SVIM/{sample}/variants.vcf",
    output:
        vcf = IN_PATH + "/SVCall/SV/{sample}.SVIM.vcf",
    run:
        shell("grep '^#' {input.vcf} > {output.vcf}")
        shell("grep '0/1' {input.vcf} >> {output.vcf}")
        shell("grep '1/1' {input.vcf} >> {output.vcf}")

# ##########################################################
rule snf:
    input:
        bam = IN_PATH + "/SVCall/mapping/{sample}_ONT.bam",
    output:
        snf = IN_PATH + "/SVCall/Sniffles2/snf/{sample}.snf",
        vcf = IN_PATH + "/SVCall/Sniffles2/snf/{sample}.vcf",
    threads:
        THREADS
    params:
        RefGenome = "/home/wuzhikun/Project/Vigna/Final/Final/Vigna_unguiculata_assembly.fasta",
        annotation = "/home/wuzhikun/Project/PanSV/Vigna_SRF/Vigna_satellite_assembly_aln_3column.bed",
    log:
        IN_PATH + "/log/snf_{sample}.log"
    run:
        ### --tandem-repeats {params.annotation}
        shell("~/anaconda3/envs/sniffles2/bin/sniffles --threads {threads} --minsupport 5 --minsvlen 50 --sample-id {wildcards.sample} --reference {params.RefGenome} --tandem-repeats {params.annotation} --input {input.bam} --vcf {output.vcf} --snf {output.snf}  --long-ins-length 50000  > {log} 2>&1")



# rule MergeVCF:
#     input:
#         snf = expand(IN_PATH + "/SVCall/Sniffles2/snf/{sample}.snf", sample=SAMPLES),
#     output:
#         vcf = IN_PATH + "/SVCall/Sniffles2/Samples_SV_merge.vcf",
#     threads:
#         THREADS
#     params:
#         RefGenome = "/home/wuzhikun/Project/Vigna/Final/Final/Vigna_unguiculata_assembly.fasta",
#         annotation = "/home/wuzhikun/Project/PanVigna/Vigna_SRF/Vigna_satellite_assembly_aln_3column.bed",
#     log:
#         IN_PATH + "/log/MergeVCF.log"
#     run:
#         ### --tandem-repeats {params.annotation}
#         Samples = " ".join(sorted(input.snf))
#         shell("~/anaconda3/envs/sniffles2/bin/sniffles --input {Samples} --vcf {output.vcf} --threads {threads}  --reference {params.RefGenome} --tandem-repeats {params.annotation} --minsupport 5 --minsvlen 50  --long-ins-length 50000 > {log} 2>&1")


# rule VCFFilt:
#     input:
#         vcf = IN_PATH + "/SVCall/Sniffles2/Samples_SV_merge.vcf",
#     output:
#         vcf = IN_PATH + "/SVCall/Sniffles2/Samples_SV_merge_filt.vcf",
#     params:
#         MergeVCFFilt = SCRIPT_DIR + "/MergeVCFFilt.py",
#         vcf = IN_PATH + "/SVCall/Sniffles2/Samples_SV_merge_filt_out.vcf",
#     run:
#         shell("python {params.MergeVCFFilt} --vcf {input.vcf} --out {output.vcf} > {params.vcf}")



# rule VCF2vgFormat:
#     input:
#         vcf = IN_PATH + "/SVCall/Sniffles2/Samples_SV_merge_filt.vcf",
#     output:
#         vcf = IN_PATH + "/SVCall/PanGraph/Samples_SV_merge_vgformat.vcf.gz",
#     params:
#         RefGenome = "/home/wuzhikun/Project/Vigna/Final/Final/Vigna_unguiculata_assembly.fasta",
#         vcf = IN_PATH + "/SVCall/PanGraph/Samples_SV_merge_vgformat.vcf",
#         VCF2vgFormat = SCRIPT_DIR + "/VCF2vgFormat.py",
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/VCF2vgFormat.log"
#     run:
#         shell("python {params.VCF2vgFormat} --vcf {input.vcf} --out {params.vcf} --fasta {params.RefGenome} > {log} 2>&1")
#         shell("bgzip {params.vcf}")
#         shell("tabix -p vcf {output.vcf}")
# #################################################################


# ############## Extract #################
# rule extract:
#     input:
#         vcf = IN_PATH + "/SVCall/PanGraph/Samples_SV_merge_vgformat.vcf",
#     output:
#         vcf = IN_PATH + "/SVCall/PanGraph/Extract/{sample}_SV_extract.vcf",
#     params:
#         PanSVExtract = SCRIPT_DIR + "/PanSVExtract.py",
#     run:
#         shell("python {params.PanSVExtract} --vcf {input.vcf} --out {output.vcf} --sample {wildcards.sample}")
        
# ######################################




# # rule simpleVG:
# #     input:
# #         vcf = IN_PATH + "/SVCall/PanGraph/Samples_SV_merge_vgformat.vcf.gz",
# #     output:
# #         vg = IN_PATH + "/SVCall/PanGraph/Samples_SV_merge_construct.vg",
# #     threads:
# #         THREADS
# #     params:
# #         RefGenome = "/home/wuzhikun/Project/Vigna/Final/Final/Vigna_unguiculata_assembly.fasta",
# #     log:
# #         IN_PATH + "/log/simpleVG.log"
# #     run:
# #         shell("vg construct -p  --threads {threads}  -r  {params.RefGenome} -v {input.vcf} > {output.vg} 2>> {log}")



# ############################### giraffe ########################
# rule GiraffeIndex:
#     input:
#         vcf = IN_PATH + "/SVCall/PanGraph/Samples_SV_merge_vgformat.vcf.gz",
#     output:
#         dist = IN_PATH + "/SVCall/PanGraph/Samples_SV_merge_giraffe.dist",
#         gbz = IN_PATH + "/SVCall/PanGraph/Samples_SV_merge_giraffe.giraffe.gbz",
#         m = IN_PATH + "/SVCall/PanGraph/Samples_SV_merge_giraffe.min",
#     params:
#         RefGenome = "/home/wuzhikun/Project/Vigna/Final/Final/Vigna_unguiculata_assembly.fasta",
#         outPrefix = IN_PATH + "/SVCall/PanGraph/Samples_SV_merge_giraffe",
#     log:
#         IN_PATH + "/log/GiraffeIndex.log"
#     run:
#         shell("/home/wuzhikun/software/vg autoindex --workflow giraffe -r {params.RefGenome} -v {input.vcf} -p {params.outPrefix} > {log} 2>&1")


# rule vgConvert:
#     input:
#         gbz = IN_PATH + "/SVCall/PanGraph/Samples_SV_merge_giraffe.giraffe.gbz",
#     output:
#         xg = IN_PATH + "/SVCall/PanGraph/Samples_SV_merge_giraffe.giraffe.xg",
#     run:
#         shell("/home/wuzhikun/software/vg convert -x {input.gbz}  > {output.xg}")


# rule snarls:
#     input:
#         gbz = IN_PATH + "/SVCall/PanGraph/Samples_SV_merge_giraffe.giraffe.gbz",
#     output:
#         snarls = IN_PATH + "/SVCall/PanGraph/Samples_SV_merge_giraffe.giraffe.snarls.pb",
#     threads:
#         THREADS
#     run:
#         shell("/home/wuzhikun/software/vg snarls -t {threads} {input.gbz} > {output.snarls}")

#############################################################



###################### call SV ###############################

# rule Giraffe:
#     input:
#         dist = IN_PATH + "/SVCall/PanGraph/Samples_SV_merge_giraffe.dist",
#         gbz = IN_PATH + "/SVCall/PanGraph/Samples_SV_merge_giraffe.giraffe.gbz",
#         m = IN_PATH + "/SVCall/PanGraph/Samples_SV_merge_giraffe.min",
#         xg = IN_PATH + "/SVCall/PanGraph/Samples_SV_merge_giraffe.giraffe.xg",
#         R1 = IN_PATH + "/clean/{sample}_NGS.R1.fastq.gz",
#         R2 = IN_PATH + "/clean/{sample}_NGS.R2.fastq.gz",
#     output:
#         gam = IN_PATH + "/SVCall/NGS/{sample}.giraffe.gam",
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/Giraffe_{sample}.log"
#     run:
#         shell("vg giraffe --gbz-name {input.gbz} --dist-name {input.dist} --minimizer-name {input.m} --xg-name {input.xg} --fastq-in {input.R1} --fastq-in {input.R2}   --threads {threads} --sample {wildcards.sample}  -o gam  >  {output.gam} 2>{log}")



rule pack:
    input:
        gbz = IN_PATH + "/SVCall/PanGraph/Samples_SV_merge_giraffe.giraffe.gbz",
        gam = IN_PATH + "/SVCall/NGS/{sample}.giraffe.gam",
    output:
        pack = IN_PATH + "/SVCall/NGS/{sample}.giraffe.pack",
    threads:
        THREADS
    log:
        IN_PATH + "/log/pack_{sample}.log"
    run:
        shell("/home/wuzhikun/software/vg pack -t {threads} -x {input.gbz} -g {input.gam} -Q5 -o {output.pack} > {log} 2>&1")




rule gamstat:
    input:
        gam = IN_PATH + "/SVCall/NGS/{sample}.giraffe.gam",
    output:
        stats = IN_PATH + "/SVCall/NGS/{sample}.giraffe.stats.txt",
    run:
        shell("/home/wuzhikun/software/vg stats -a {input.gam} > {output.stats}")


rule call:
    input:
        gbz = IN_PATH + "/SVCall/PanGraph/Samples_SV_merge_giraffe.giraffe.gbz",
        snarls = IN_PATH + "/SVCall/PanGraph/Samples_SV_merge_giraffe.giraffe.snarls.pb",
        pack = IN_PATH + "/SVCall/NGS/{sample}.giraffe.pack",
    output:
        vcf = IN_PATH + "/SVCall/NGS/{sample}.giraffe.vcf.gz",
    params:
        vcf = IN_PATH + "/SVCall/NGS/{sample}.giraffe.vcf",
    threads:
        THREADS
    run:
        shell("/home/wuzhikun/software/vg call -t {threads} -k {input.pack} --min-support 2,2 -z -a -r {input.snarls} -z -s {wildcards.sample} {input.gbz} > {params.vcf}")
        shell("bgzip {params.vcf}")
        shell("tabix -p vcf {output.vcf}")



# # rule CallFilt:
# #     input:
# #         vcf = IN_PATH + "/SVCall/NGS/{sample}.giraffe.vcf.gz",
# #     output:
# #         vcf0 = temp(IN_PATH + "/SVCall/NGS/{sample}.giraffe.temp.vcf"),
# #         vcf = IN_PATH + "/SVCall/NGS/{sample}.giraffe.filt.vcf",
# #     params:
# #         NGSgiraffeFilt = SCRIPT_DIR + "/NGSgiraffeFilt.py",
# #     threads:
# #         THREADS
# #     log:
# #         IN_PATH + "/log/CallFilt_{sample}.log"
# #     run:
# #         shell("pigz -p {threads} -dc {input.vcf} > {output.vcf0}")
# #         shell("python {params.NGSgiraffeFilt} --vcf {output.vcf0} --out {output.vcf} --altThreshold 10 > {log} 2>&1")



# rule PanSV:
#     input:
#         vcf = expand(IN_PATH + "/SVCall/NGS/{sample}.giraffe.vcf.gz", sample=SAMPLES),
#     output:
#         vcf = IN_PATH + "/SVCall/PanSV/Dsamples.giraffe.sv.vcf.gz",
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/PanSV.log"
#     run:
#         Files = " ".join(sorted(input.vcf))
#         shell("bcftools merge --threads {threads} -Oz -o {output.vcf} {Files} > {log} 2>&1")


# rule TwoAllele:
#     input:
#         vcf = IN_PATH + "/SVCall/PanSV/Dsamples.giraffe.sv.vcf.gz",   
#     output:
#         vcf = IN_PATH + "/SVCall/PanSV/Dsamples.giraffe.sv.vcf",
#         recode = IN_PATH + "/SVCall/PanSV/Dsamples.twoallel.recode.vcf",
#         gt = IN_PATH + "/SVCall/PanSV/Dsamples.twoallel.GT.vcf",
#     params:
#         outPrefix = IN_PATH + "/SVCall/PanSV/Dsamples.twoallel",
#     log:
#         IN_PATH + "/log/TwoAllele.log",
#     run:
#         shell("pigz -dc {input.vcf} > {output.vcf}")
#         shell("vcftools --vcf {output.vcf} --min-alleles 2 --max-alleles 2   --max-missing 0.2 --out {params.outPrefix} --recode > {log} 2>&1")
#         shell("bcftools annotate -x ^FORMAT/GT {output.recode} > {output.gt}")
        



# rule SVMAFAllele:
#     input:
#         vcf = IN_PATH + "/SVCall/PanSV/Dsamples.twoallel.GT.vcf",  
#     output:
#         vcf = IN_PATH + "/SVCall/PanSV/Dsamples.twoallel.maf005.recode.vcf",
#     params:
#         outPrefix = IN_PATH + "/SVCall/PanSV/Dsamples.twoallel.maf005",
#     log:
#         IN_PATH + "/log/SVMAFAllele.log",
#     run:
#         shell("vcftools --vcf {input.vcf}  --maf 0.05 --out {params.outPrefix} --recode > {log} 2>&1")



# rule SVConvertPed:
#     input:
#         vcf = IN_PATH + "/SVCall/PanSV/Dsamples.twoallel.maf005.recode.vcf",
#     output:
#         Ped = IN_PATH + "/SVCall/PanSV/Dsamples.twoallel.maf005.ped",
#         Map = IN_PATH + "/SVCall/PanSV/Dsamples.twoallel.maf005.map",
#     params:
#         outPrefix = IN_PATH + "/SVCall/PanSV/Dsamples.twoallel.maf005",
#     log:
#         IN_PATH + "/log/SVConvertPed.log",
#     run:
#         shell("vcftools --vcf {input.vcf} --plink --out {params.outPrefix} > {log} 2>&1")


# def map_id(vcf_file, map_file, out_file):
#     IDRegion = {}
#     in_h = open(vcf_file, "r")
#     for line in in_h:
#         line = line.strip()
#         if line.startswith("#"):
#             continue
#         else:
#             lines = line.split("\t")
#             Chr, Pos, Id = lines[:3]
#             IDRegion[Id] = (Chr, Pos)
#     in_h.close()
#     map_h = open(map_file, "r")
#     out_h = open(out_file, "w")
#     for line in map_h:
#         lines = line.strip().split("\t")
#         mapid = lines[1]
#         region = IDRegion[mapid]
#         c, p = region
#         newline = "%s\t%s:%s\t0\t%s" % (c, c, p, p) 
#         out_h.write("%s\n" % newline)
#     map_h.close()
#     out_h.close()


# rule updateMap:
#     input:
#         vcf = IN_PATH + "/SVCall/PanSV/Dsamples.twoallel.maf005.recode.vcf",
#         Map = IN_PATH + "/SVCall/PanSV/Dsamples.twoallel.maf005.map",
#     output:
#         Map = IN_PATH + "/SVCall/PanSV/Dsamples.twoallel.maf005-1.map",
#     run:
#         map_id(input.vcf, input.Map, output.Map)
#         # shell("cp {output.Map} {input.Map}")


# rule SVmakebed:
#     input:
#         Map0 = IN_PATH + "/SVCall/PanSV/Dsamples.twoallel.maf005-1.map",
#         Ped = IN_PATH + "/SVCall/PanSV/Dsamples.twoallel.maf005.ped",
#         Map = IN_PATH + "/SVCall/PanSV/Dsamples.twoallel.maf005.map",
#     output:
#         bed = IN_PATH + "/SVCall/PanSV/Dsamples.twoallel.maf005.bed",
#         fam = IN_PATH + "/SVCall/PanSV/Dsamples.twoallel.maf005.fam",
#     params:
#         outPrefix = IN_PATH + "/SVCall/PanSV/Dsamples.twoallel.maf005",
#     log:
#         IN_PATH + "/log/SVmakebed.log",
#     run:
#         shell("plink --file {params.outPrefix} --make-bed --out {params.outPrefix} --allow-extra-chr  > {log} 2>&1")

rule modifySVfam:
    input:
        fam = IN_PATH + "/SVCall/PanSV/Dsamples.twoallel.maf005.fam",
    output:
        fam = IN_PATH + "/SVCall/PanSV/Dsamples.twoallel.maf005-1.fam",
    run:
        cmd = """awk '{print $1" "$2" "$3" "$4" "$5" "1}' %s > %s """ % (input.fam, output.fam)
        os.system(cmd)
        shell("cp {output.fam} {input.fam}")



rule kinship:
    input:
        bed = IN_PATH + "/SVCall/PanSV/Dsamples.twoallel.maf005.bed",
        fam = IN_PATH + "/SVCall/PanSV/Dsamples.twoallel.maf005-1.fam",
    output:
        ks = IN_PATH + "/SVCall/PanSV/GWASGemma/output/geno_kS.cXX.txt",
    params:
        inPrefix = IN_PATH + "/SVCall/PanSV/Dsamples.twoallel.maf005",
        outDir = IN_PATH + "/SVCall/PanSV/GWASGemma",
        outPrefix = "geno_kS",
    log:
        IN_PATH + "/log/kinship.log", 
    run:
        if not os.path.exists(params.outDir):
            os.makedirs(params.outDir)
        # cmd = "source activate WGS && cd %s && gemma -bfile %s -gk 1 -o %s > %s 2>&1" % (params.outDir, params.inPrefix, params.outPrefix, log)
        # print(cmd)
        # os.system(cmd)
        shell("cd {params.outDir} && gemma -bfile {params.inPrefix} -gk 1 -o {params.outPrefix} > {log} 2>&1")




rule SVassociation:
    input:
        bed = IN_PATH + "/SVCall/PanSV/Dsamples.twoallel.maf005.bed",
        ks = IN_PATH + "/SVCall/PanSV/GWASGemma/output/geno_kS.cXX.txt",
        pheno1 = IN_PATH + "/SVCall/PanSV/phenotype/{trait}.txt",
        # pheno1 = IN_PATH + "/NGS/phenotype/Glucose_pheno.txt",
        # pheno2 = IN_PATH + "/NGS/phenotype/Levulose_pheno.txt",
        # pheno3 = IN_PATH + "/NGS/phenotype/Sucrose_pheno.txt",
    output:
        asso1 = IN_PATH + "/SVCall/PanSV/GWASGemma/output/{trait}.assoc.txt",
    threads:
        THREADS
    params:
        inPrefix = IN_PATH + "/SVCall/PanSV/Dsamples.twoallel.maf005",
        outDir = IN_PATH + "/SVCall/PanSV/GWASGemma",
    log:
        IN_PATH + "/log/SVassociation_{trait}.log", 
    run:
        ### gemma -bfile /home/wuzhikun/Project/Population/population/plink/Sample_SV_geno -k /home/wuzhikun/Project/Population/population/GWAS/output/SV_geno_kS.cXX.txt -lmm 1 -o pheno1 -p /home/wuzhikun/Project/Population/phenotype/CN329_phenotype_value-1.txt
        # shell("cut -f {wildcards.trait} {input.pheno} > {output.pheno}")
        shell("cd {params.outDir} && gemma -bfile {params.inPrefix} -k {input.ks} -lmm 1 -o {wildcards.trait} -p {input.pheno1} > {log} 2>&1")




def changeChrString(in_file, out_file):
    in_h = open(in_file, "r")
    out_h = open(out_file, "w")
    header = in_h.readline().strip()
    out_h.write("%s\n" % header)
    for line in in_h:
        lines = line.strip().split("\t")
        Chr = lines[1].strip("Vu").lstrip("0")
        out_h.write("%s\t%s\t%s\t%s\n" % (lines[0], Chr, lines[2], lines[3]))
    in_h.close()
    out_h.close()


rule SVManhattanPlot0:
    input:
        assoc = IN_PATH + "/SVCall/PanSV/GWASGemma/output/{trait}.assoc.txt",
    output:
        assoc = IN_PATH + "/SVCall/PanSV/GWASGemma/output/plot/{trait}.assoc.txt",
        assoc1 = IN_PATH + "/SVCall/PanSV/GWASGemma/output/plot/{trait}.assoc-1.txt",
    run:
        shell("echo 'SNP CHR BP P' | sed 's/ /\t/g'   > {output.assoc}")
        cmd = """sed '1d' %s | grep -v "_un" | awk '{print $2"\t"$1"\t"$3"\t"$12}' >> %s"""  % (input.assoc, output.assoc)
        print(cmd)
        os.system(cmd)
        changeChrString(output.assoc, output.assoc1)



rule ManhattanPlot1:
    input:
        assoc = IN_PATH + "/SVCall/PanSV/GWASGemma/output/plot/{trait}.assoc-1.txt",
    output:
        manhattan = IN_PATH + "/SVCall/PanSV/GWASGemma/output/plot/{trait}.assoc_mahattan_SV.jpeg",
        qq = IN_PATH + "/SVCall/PanSV/GWASGemma/output/plot/{trait}.assoc_qq_SV.jpeg",
    params:
        GWASPlot = SCRIPT_DIR + "/GWASPlot.R",
    log:
        IN_PATH + "/log/SVManhattanPlot_{trait}.log", 
    run:
        shell("Rscript {params.GWASPlot} --input {input.assoc} --manhattan {output.manhattan} --qq {output.qq} > {log} 2>&1")


rule SVSigLoci:
    input:
        asso = IN_PATH + "/SVCall/PanSV/GWASGemma/output/plot/{trait}.assoc-1.txt",
    output:
        sig = IN_PATH + "/SVCall/PanSV/GWASGemma/output/plot/{trait}.assoc_sig.txt",
    run:
        cmd = "awk '{if ($4 < 2.9e-5){print $0}}'  %s > %s" % (input.asso, output.sig)
        print(cmd)
        os.system(cmd)


##############################################################


# rule MergeSV:
#     input:
#         vcf = expand(IN_PATH + "/SVCall/cuteSV/{sample}.cutesv.vcf", sample=SAMPLES),
#     output:
#         fileList = IN_PATH + "/SVCall/Merge/samples_file_list_DEL.txt",
#         vcf = IN_PATH + "/SVCall/Merge/Sample_SVs_jasmine.vcf",
#     threads:
#         THREADS * ThreadFold
#     log:
#         IN_PATH + "/log/MergeSV.log",
#     run:
#         Files = "\n".join(sorted(input.vcf))
#         out_h = open(output.fileList, "w")
#         out_h.write("%s\n" % Files)
#         out_h.close()
#         # shell("jasmine --output_genotypes --ignore_strand --keep_var_ids threads={threads} file_list={output.fileList} out_file={output.vcf} > {log} 2>&1")
#         cmd = "source activate nanovar && jasmine --output_genotypes --ignore_strand --keep_var_ids threads=%s file_list=%s out_file=%s > %s 2>&1" % (threads, output.fileList, output.vcf, log)
#         os.system(cmd)







rule straglr:
    input:
        bam = IN_PATH + "/SVCall/mapping/{sample}_ONT.bam",
    output:
        TR = IN_PATH + "/TR/straglr/{sample}/{chr}_straglr_TR.tsv",
        bed = IN_PATH + "/TR/straglr/{sample}/{chr}_straglr_TR.bed",
    threads:
        THREADS
    params:
        RefGenome = "/home/wuzhikun/Project/Vigna/Final/Final/Vigna_unguiculata_assembly.fasta",
        outdir = IN_PATH + "/TR/straglr/{sample}",
        tmpdir = IN_PATH + "/TR/Temp/{sample}",
        TR = IN_PATH + "/TR/straglr/{sample}/{chr}_straglr_TR",
    log:
        IN_PATH + "/log/straglr_{sample}_{chr}.log",
    run:
        ### --exclude hg38.exclude.bed
        ### source actiuvate HiC-Pro
        if not os.path.exists(params.outdir):
            os.makedirs(params.outdir)
        if not os.path.exists(params.tmpdir):
            os.makedirs(params.tmpdir)
        # shell("source activate HiC-Pro && cd {params.outdir} && straglr.py {input.bam}  {params.RefGenome} {params.TR}  --min_str_len 2 --max_str_len 100 --min_ins_size 50 --min_support 5  --chroms {wildcards.chr} --nprocs {threads} --tmpdir {params.tmpdir}  --genotype_in_size > {log} 2>&1")
        cmd = "source activate HiC-Pro && cd %s && straglr.py %s  %s %s  --min_str_len 2 --max_str_len 100 --min_ins_size 50 --min_support 5  --chroms %s --nprocs %s --tmpdir %s  --genotype_in_size > %s 2>&1" % (params.outdir, input.bam, params.RefGenome, params.TR, wildcards.chr, threads, params.tmpdir, log)
        print(cmd)
        os.system(cmd)



rule mergeRT:
    input:
        TR = expand(IN_PATH + "/TR/straglr/{sample}/{chr}_straglr_TR.tsv", sample=SAMPLES, chr=CHRS),
    output:
        TR = IN_PATH + "/TR/straglr/{sample}_all_chr_straglr_TR.txt",
    run:
        Target = []
        files = sorted(input.TR)
        for i in files:
            s = i.split("/")[-2]
            if s == wildcards.sample:
                Target.append(i)
        print(Target)
        for j in range(len(Target)):
            if j == 0:
                cmd = "cat %s > %s" % (Target[j], output.TR)
            else:
                cmd = "grep -v  '^#' %s >> %s" % (Target[j], output.TR)
            os.system(cmd)



rule mergeRT2:
    input:
        bed = expand(IN_PATH + "/TR/straglr/{sample}/{chr}_straglr_TR.bed", sample=SAMPLES, chr=CHRS),
    output:
        bed = IN_PATH + "/TR/straglr/{sample}_all_chr_straglr_TR.bed",
    run:
        Target = []
        files = sorted(input.bed)
        for i in files:
            s = i.split("/")[-2]
            if s == wildcards.sample:
                Target.append(i)
        print(Target)
        for j in range(len(Target)):
            if j == 0:
                cmd = "cat %s > %s" % (Target[j], output.bed)
            else:
                cmd = "grep -v  '^#' %s >> %s" % (Target[j], output.bed)
            os.system(cmd)

################################
