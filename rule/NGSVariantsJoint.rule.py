# #################### joint call ###########
# rule tabix:
#     input:
#         vcf = IN_PATH + '/NGS/VCF/{sample}/{sample}_{Chr}.g.vcf.gz',
#     output:
#         tbi = IN_PATH + '/NGS/VCF/{sample}/{sample}_{Chr}.g.vcf.gz.tbi',
#     run:
#         shell("tabix -p vcf {input.vcf}")


# rule GenomicsDBImport:
#     input:
#         gvcf = expand(IN_PATH + '/NGS/VCF/{sample}/{sample}_{Chr}.g.vcf.gz', sample=SAMPLES, Chr=CHRS),
#     output:
#         genomicsdb = directory(IN_PATH + '/NGS/Joint/genomicsdb/{Chr}_genomicsdb'),
#     params:
#         GATK4 = config['GATK4'],
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/{Chr}.GenomicsDBImport.log",
#     run:
#         Files = input.gvcf
#         Targets = []
#         for f in Files:
#             if wildcards.Chr in f:
#                 Targets.append(f)
#         GVCFs = " -V ".join(sorted(Targets))
#         print(sorted(Targets))
#         shell("gatk GenomicsDBImport   -V {GVCFs} --genomicsdb-workspace-path {output.genomicsdb} -L {wildcards.Chr} > {log} 2>&1")
    

# rule jointCall:
#     input:
#         genomicsdb = IN_PATH + '/NGS/Joint/genomicsdb/{Chr}_genomicsdb',
#     output:
#         vcf = IN_PATH + '/NGS/Joint/VCF/{Chr}_Fengchan6_joint_variant.vcf',
#     params:
#         GATK4 = config['GATK4'],
#         REF = config['RefGenome'],
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/{Chr}.jointCall.log",
#     run:
#         shell("gatk GenotypeGVCFs  -R {params.REF} -V gendb://{input.genomicsdb} -O {output.vcf} --allow-old-rms-mapping-quality-annotation-data > {log} 2>&1")

# ##########################################################



#################### joint call2 ###########
rule tabix2:
    input:
        vcf = IN_PATH + '/NGS/Joint/gVCF/{sample}_NGS.g.vcf.gz',
    output:
        tbi = IN_PATH + '/NGS/Joint/gVCF/{sample}_NGS.g.vcf.gz.tbi',
    run:
        shell("tabix -p vcf {input.vcf}")


rule GenomicsDBImport2:
    input:
        gvcf = expand(IN_PATH + '/NGS/Joint/gVCF/{sample}_NGS.g.vcf.gz', sample=SAMPLES),
    output:
        genomicsdb = directory(IN_PATH + '/NGS/Joint/genomicsdb/{Chr}_genomicsdb'),
    params:
        GATK4 = config['GATK4'],
    threads:
        THREADS
    log:
        IN_PATH + "/log/{Chr}.GenomicsDBImport2.log",
    run:
        GVCFs = " -V ".join(sorted(input.gvcf))
        print(GVCFs)
        shell("gatk GenomicsDBImport   -V {GVCFs} --genomicsdb-workspace-path {output.genomicsdb} -L {wildcards.Chr} > {log} 2>&1")
    

rule jointCall2:
    input:
        genomicsdb = IN_PATH + '/NGS/Joint/genomicsdb/{Chr}_genomicsdb',
    output:
        vcf = IN_PATH + '/NGS/Joint/VCF/G98_{Chr}_joint_variant.vcf',
    params:
        GATK4 = config['GATK4'],
        REF = config['RefGenome'],
    threads:
        THREADS
    log:
        IN_PATH + "/log/{Chr}.jointCall2.log",
    run:
        shell("gatk GenotypeGVCFs  -R {params.REF} -V gendb://{input.genomicsdb} -O {output.vcf} --allow-old-rms-mapping-quality-annotation-data > {log} 2>&1")




rule GenomicsDBImportContig2:
    input:
        gvcf = expand(IN_PATH + '/NGS/Joint/gVCF/{sample}_NGS.g.vcf.gz', sample=SAMPLES),
    output:
        genomicsdb = directory(IN_PATH + '/NGS/Joint/genomicsdb/Contigs_genomicsdb'),
    params:
        GATK4 = config['GATK4'],
        # Chrs = "Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11",
        contig = "/home/wuzhikun/database/genome/Vigna_unguiculata/G98/Lachesis_assembly_changed.contig.bed",
    threads:
        THREADS
    log:
        IN_PATH + "/log/contig.GenomicsDBImport2.log",
    run:
        GVCFs = " -V ".join(sorted(input.gvcf))
        print(GVCFs)
        shell("gatk GenomicsDBImport   -V {GVCFs} --genomicsdb-workspace-path {output.genomicsdb} -L {params.contig} > {log} 2>&1")


rule jointCallContig:
    input:
        genomicsdb = IN_PATH + '/NGS/Joint/genomicsdb/Contigs_genomicsdb',
    output:
        vcf = IN_PATH + '/NGS/Joint/VCF/G98_Contigs_joint_variant.vcf',
    params:
        GATK4 = config['GATK4'],
        REF = config['RefGenome'],
    threads:
        THREADS
    log:
        IN_PATH + "/log/contigs.jointCall.log",
    run:
        shell("gatk GenotypeGVCFs  -R {params.REF} -V gendb://{input.genomicsdb} -O {output.vcf} --allow-old-rms-mapping-quality-annotation-data  > {log} 2>&1")

##########################################################





# rule MergeChr:
#     input:
#         vcf = expand(IN_PATH + '/NGS/Joint/VCF/{Chr}_joint_variant.vcf', Chr=CHRS),
#         contig = IN_PATH + '/NGS/Joint/VCF/Contigs_joint_variant.vcf',
#     output:
#         vcf = IN_PATH + '/NGS/Joint/VCF/Samples_Chrs_joint_variant.vcf',
#     run:
#         Chrs = " ".join(sorted(input.vcf))
#         shell("bcftools concat -o {output.vcf} {Chrs} {input.contig} ")
        




################## test #############
# rule JointCallContig:
#     input:
#         gvcf = expand(IN_PATH + '/NGS/Joint/gVCF/{sample}_NGS.g.vcf.gz', sample=SAMPLES),
#     output:
#         vcf = IN_PATH + '/NGS/GATK3/Samples_joint_call.Contig.vcf',
#     params:
#         GATK3 = config['GATK3'],
#         REF = config['RefGenome'],
#         Memory = "-Xmx50g",
#         contig = "/home/wuzhikun/database/genome/Vigna_unguiculata/G98/Lachesis_assembly_changed.contig.bed",
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/JointCallContig.log",
#     run:
#         Files = sorted(input.gvcf)
#         IN = "  --variant ".join(Files)
#         shell("java {params.Memory} -jar {params.GATK3} -T GenotypeGVCFs -L {params.contig} -nt {threads} -R {params.REF} --variant {IN} -o {output.vcf} > {log} 2>&1")
####################################







######################## joint call ######################
### https://nbisweden.github.io/workshop-ngsintro/2105/lab_vc.html#2_Variant_calling_in_cohort

# rule JointChrs:
#     input:
#         tbi = expand(IN_PATH + '/NGS/Joint/gVCF/{sample}_NGS.g.vcf.gz.tbi', sample=SAMPLES),
#         gvcf = expand(IN_PATH + '/NGS/Joint/gVCF/{sample}_NGS.g.vcf.gz', sample=SAMPLES),
#     output:
#         gvcf = IN_PATH + '/NGS/DSamples/Samples_jointcall_{Chr}.g.vcf',    
#     params:
#         REF = config['RefGenome'],
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/JointChrs_{Chr}.log",  
#     run:
#         Files = sorted(input.gvcf)
#         IN = "  -V ".join(Files)
#         shell("gatk --java-options -Xmx20g CombineGVCFs -R {params.REF} -V {IN} -O {output.gvcf} -L {wildcards.Chr} > {log} 2>&1")  




# rule JointContig:
#     input:
#         gvcf = expand(IN_PATH + '/NGS/Joint/gVCF/{sample}_NGS.g.vcf.gz', sample=SAMPLES),
#     output:
#         gvcf = IN_PATH + '/NGS/GATK3/Samples_joint_call.Contig.g.vcf',    
#     params:
#         # GATK3 = config['GATK3'],
#         REF = config['RefGenome'],
#         # Memory = "-Xmx50g",
#         contig = "/home/wuzhikun/database/genome/Vigna_unguiculata/G98/Lachesis_assembly_changed.contig.bed",    
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/JointContig.log",  
#     run:
#         Files = sorted(input.gvcf)
#         IN = "  -V ".join(Files)
#         shell("gatk --java-options -Xmx20g CombineGVCFs -R {params.REF} -V {IN} -O {output.gvcf} -L {params.contig} > {log} 2>&1")  




# rule JointCallContig:
#     input:
#         gvcf = IN_PATH + '/NGS/GATK3/Samples_joint_call.Contig.g.vcf', 
#     output:
#         vcf = temp(IN_PATH + '/NGS/GATK3/Samples_joint_call.variant.Contig.vcf'),
#     params:
#         REF = config['RefGenome'], 
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/JointCallContig.log",  
#     run:
#         shell("gatk --java-options -Xmx20g  GenotypeGVCFs -R {params.REF} -V {input.gvcf} -O {output.vcf} --allow-old-rms-mapping-quality-annotation-data > {log} 2>&1")


# rule JointCallChrs:
#     input:
#         gvcf = IN_PATH + '/NGS/DSamples/Samples_jointcall_{Chr}.g.vcf', 
#     output:
#         vcf = temp(IN_PATH + '/NGS/DSamples/Sample_jointcall.variant.{Chr}.vcf'),
#     params:
#         REF = config['RefGenome'], 
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/JointCallChr_{Chr}.log",  
#     run:
#         shell("gatk --java-options -Xmx20g  GenotypeGVCFs -R {params.REF} -V {input.gvcf} -L {wildcards.Chr}  -O {output.vcf} --allow-old-rms-mapping-quality-annotation-data  > {log} 2>&1")


# rule vcfbgip:
#     input:
#         vcf = IN_PATH + '/NGS/DSamples/Sample_jointcall.variant.{Chr}.vcf',
#     output:
#         vcf = IN_PATH + '/NGS/DSamples/Sample_jointcall.variant.{Chr}.vcf.gz',
#     run:
#         shell("bgzip {input.vcf}")
#         shell("tabix -p vcf {output.vcf}")


# rule contigbgip:
#     input:
#         vcf = IN_PATH + '/NGS/GATK3/Samples_joint_call.variant.Contig.vcf'
#     output:
#         vcf = IN_PATH + '/NGS/GATK3/Samples_joint_call.variant.Contig.vcf.gz',
#     run:
#         shell("bgzip {input.vcf}")
#         shell("tabix -p vcf {output.vcf}")


rule MergeChrVCF:
    input:
        vcf = expand(IN_PATH + '/NGS/DSamples/Sample_jointcall.variant.{Chr}.vcf.gz', Chr=CHRS),
        # contig = IN_PATH + '/NGS/GATK3/Samples_joint_call.variant.Contig.vcf.gz',
    output:
        vcf = IN_PATH + '/NGS/DSamples/Sample_all_chrom_variant.vcf.gz',
    run:
        # Files = " ".join(sorted(input.vcf) + [input.contig])
        Files = " ".join(sorted(input.vcf))
        shell("bcftools concat -Oz -o {output.vcf} {Files} ")
        shell("tabix -p vcf {output.vcf}")

##############################################################################



#################### joint call SNV ##############
rule JointChrs3:
    input:
        tbi = expand(IN_PATH + '/NGS/Joint/gVCF/{sample}_NGS.g.vcf.gz.tbi', sample=SAMPLES),
        gvcf = expand(IN_PATH + '/NGS/Joint/gVCF/{sample}_NGS.g.vcf.gz', sample=SAMPLES),
    output:
        gvcf = IN_PATH + '/NGS/AllSamples/Samples_jointcall_{Chr}.g.vcf',    
    params:
        REF = "/home/wuzhikun/Project/Vigna/Final/Final/Vigna_unguiculata_assembly.fasta",
    threads:
        THREADS
    log:
        IN_PATH + "/log/JointChrs_{Chr}.log",  
    run:
        Files = sorted(input.gvcf)
        IN = "  -V ".join(Files)
        shell("gatk --java-options -Xmx20g CombineGVCFs -R {params.REF} -V {IN} -O {output.gvcf} -L {wildcards.Chr} > {log} 2>&1") 



rule JointCallChrs3:
    input:
        gvcf = IN_PATH + '/NGS/AllSamples/Samples_jointcall_{Chr}.g.vcf', 
    output:
        vcf = temp(IN_PATH + '/NGS/AllSamples/Sample_jointcall.variant.{Chr}.vcf'),
    params:
        REF = "/home/wuzhikun/Project/Vigna/Final/Final/Vigna_unguiculata_assembly.fasta",
    threads:
        THREADS
    log:
        IN_PATH + "/log/JointCallChr_{Chr}.log",  
    run:
        shell("gatk --java-options -Xmx20g  GenotypeGVCFs -R {params.REF} -V {input.gvcf} -L {wildcards.Chr}  -O {output.vcf} --allow-old-rms-mapping-quality-annotation-data  > {log} 2>&1")


rule vcfbgip3:
    input:
        vcf = IN_PATH + '/NGS/AllSamples/Sample_jointcall.variant.{Chr}.vcf',
    output:
        vcf = IN_PATH + '/NGS/AllSamples/Sample_jointcall.variant.{Chr}.vcf.gz',
    run:
        shell("bgzip {input.vcf}")
        shell("tabix -p vcf {output.vcf}")



rule MergeChrVCF3:
    input:
        vcf = expand(IN_PATH + '/NGS/AllSamples/Sample_jointcall.variant.{Chr}.vcf.gz', Chr=CHRS),
    output:
        vcf = IN_PATH + '/NGS/AllSamples/Sample_all_chrom_variant.vcf.gz',
    run:
        Files = " ".join(sorted(input.vcf))
        shell("bcftools concat -Oz -o {output.vcf} {Files} ")
        shell("tabix -p vcf {output.vcf}")


#################################################################





