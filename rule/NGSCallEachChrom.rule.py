
# ######################################################################
# rule HaplotypeCaller:
#     input:
#         bam = IN_PATH + '/NGS/mapping/{sample}/{sample}.sorted.bam',
#     output:
#         vcf = temp(IN_PATH + '/NGS/VCF/{sample}/{sample}_{Chr}.g.vcf'),
#         # idx = IN_PATH + '/VCF/{sample}_{Chr}.gvcf.idx',
#     params:
#         GATK4 = config['GATK4'],
#         REF = config['RefGenome'],
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/{sample}_{Chr}.HaplotypeCaller.log",
#     run:
#         ### -L 20:10,000,000-10,200,000
#         ### -L hg38.exome.regions.bed
#         ### --dbsnp dbSNP.vcf
#         ### --variant_index_type LINEAR
#         ### --variant_index_parameter 128000
#         ##### -nt / --num_threads controls the number of data threads sent to the processor (acting at the machine level)
#         ##### -nct / --num_cpu_threads_per_data_thread controls the number of CPU threads allocated to each data thread (acting at the core level)
#         ###  -variant_index_type LINEAR -variant_index_parameter 128000 
#         shell('java -Xmx10g -jar {params.GATK4} HaplotypeCaller -R {params.REF} -L {wildcards.Chr} -I {input.bam} -O {output.vcf} -ERC GVCF --native-pair-hmm-threads {threads}  >{log} 2>&1') 
#         ### -gt_mode DISCOVERY -stand_call_conf 30  -stand_emit_conf 30



# rule bgzip:
#     input:
#         vcf = IN_PATH + '/NGS/VCF/{sample}/{sample}_{Chr}.g.vcf',
#     output:
#         vcf = IN_PATH + '/NGS/VCF/{sample}/{sample}_{Chr}.g.vcf.gz',
#     run:
#         shell("bgzip {input.vcf}")
#         # shell("tabix -p vcf {output.vcf}")


# rule bgzip2:
#     input:
#         vcf = IN_PATH + '/NGS/VCF/{sample}/{sample}_{Chr}.g.vcf.gz',
#     output:
#         vcf = IN_PATH + '/NGS/VCF/{sample}/{sample}_{Chr}.g.vcf.gz.tbi',
#     run:
#         shell("tabix -p vcf {input.vcf}")





# rule HaplotypeCallerContig:
#     input:
#         bam = IN_PATH + '/NGS/mapping/{sample}/{sample}.sorted.bam',
#     output:
#         vcf = temp(IN_PATH + '/NGS/VCF/{sample}/{sample}.contig.g.vcf'),
#     params:
#         GATK4 = config['GATK4'],
#         REF = config['RefGenome'],
#         contig = "/home/wuzhikun/Project/PanSV/BACKUP/Assembly/Fengchan6/Vigna_unguiculata_assembly.contig.bed",
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/{sample}_HaplotypeCallerContig.log",
#     run:
#         shell('java -Xmx10g -jar {params.GATK4} HaplotypeCaller -R {params.REF} -L {params.contig} -I {input.bam} -O {output.vcf} -ERC GVCF --native-pair-hmm-threads {threads}  >{log} 2>&1') 
#         ### -gt_mode DISCOVERY -stand_call_conf 30  -stand_emit_conf 30


# rule bgzipContig:
#     input:
#         vcf = IN_PATH + '/NGS/VCF/{sample}/{sample}.contig.g.vcf',
#     output:
#         vcf = IN_PATH + '/NGS/VCF/{sample}/{sample}.contig.g.vcf.gz',
#     run:
#         shell("bgzip {input.vcf}")
#         shell("tabix -p vcf {output.vcf}")

##################################################


######################### concat chromosomes  ##############
rule concat:
    input:
        vcf = expand(IN_PATH + '/NGS/VCF/{sample}/{sample}_{Chr}.g.vcf.gz', sample=SAMPLES, Chr=CHRS),
    output:
        vcf = IN_PATH + '/NGS/Joint/gVCF/{sample}_NGS.g.vcf.gz',
    run:
        target = []
        files = input.vcf
        t = wildcards.sample + "_"
        for f in files:
            if t in f:
                target.append(f)
        sortTarget = " ".join(sorted(target))
        shell("bcftools concat -Oz -o {output.vcf} {sortTarget} ")



#############################################################

