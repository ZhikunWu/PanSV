

# # ################## Variant #########################

# # rule sekectVariant:
# #     input:
# #         # vcf = IN_PATH + "/NGS/GATK3/Sample_all_chrom_variant.vcf.gz",
# #         vcf = IN_PATH + "/NGS/DSamples/Sample_all_chrom_variant.vcf.gz",
# #     output:
# #         snp = temp(IN_PATH + "/NGS/Variant/Sample_jointcall.snp.vcf"),
# #         indel = temp(IN_PATH + "/NGS/Variant/Sample_jointcall.indel.vcf"),
# #     params:
# #         REF = config['RefGenome'], 
# #     log:
# #         IN_PATH + "/log/sekectVariant.log",
# #     run:
# #         shell("gatk --java-options -Xmx30g SelectVariants -R {params.REF} -V {input.vcf} --select-type-to-include SNP  -O {output.snp} > {log} 2>&1")
# #         shell("gatk --java-options -Xmx30g SelectVariants -R {params.REF} -V {input.vcf} --select-type-to-include INDEL  -O {output.indel} >> {log} 2>&1")


# # rule VariantFilt:
# #     input:
# #         snp = IN_PATH + "/NGS/Variant/Sample_jointcall.snp.vcf",
# #         indel = IN_PATH + "/NGS/Variant/Sample_jointcall.indel.vcf",
# #     output:
# #         snp = temp(IN_PATH + "/NGS/Variant/Sample_jointcall.snp.filt.vcf"),
# #         indel = temp(IN_PATH + "/NGS/Variant/Sample_jointcall.indel.filt.vcf"),
# #     params:
# #         REF = config['RefGenome'], 
# #     log:
# #         IN_PATH + "/log/VariantFilt.log",
# #     run:
# #         shell('gatk --java-options -Xmx30g VariantFiltration -R {params.REF} -V {input.snp} -O {output.snp} --filter-name QDfilter --filter-expression "QD < 2.0" --filter-name MQfilter --filter-expression "MQ < 40.0" --filter-name FSfilter --filter-expression "FS > 60.0"  > {log} 2>&1')
# #         shell('gatk --java-options -Xmx30g VariantFiltration -R {params.REF} -V {input.indel} -O {output.indel} --filter-name QDfilter --filter-expression "QD < 2.0" --filter-name FSfilter --filter-expression "FS > 200.0"  >> {log} 2>&1')

# # rule MergeVariant:
# #     input:
# #         snp = IN_PATH + "/NGS/Variant/Sample_jointcall.snp.filt.vcf",
# #         indel = IN_PATH + "/NGS/Variant/Sample_jointcall.indel.filt.vcf",
# #     output:
# #         variant = IN_PATH + "/NGS/Variant/Sample_jointcall.variant.filt.vcf",
# #     log:
# #         IN_PATH + "/log/MergeVariant.log",
# #     run:
# #         shell("gatk --java-options -Xmx30g MergeVcfs -I {input.snp} -I {input.indel} -O {output.variant} > {log} 2>&1")


# # ####### the unmber is identical, 过滤无效？


# # # rule RemoveIndividual:
# # #     input:
# # #         vcf = IN_PATH + "/NGS/Variant/Sample_jointcall.variant.filt.vcf",
# # #     output:
# # #         vcf = IN_PATH + "/NGS/Variant/Sample_jointcall.variant.removeIndiv.recode.vcf",    
# # #     params:
# # #         outPrefix = IN_PATH + "/NGS/Variant//Sample_jointcall.variant.removeIndiv",
# # #         removeIndiv = IN_PATH + "/NGS/remove_individual.txt",
# # #     log:
# # #         IN_PATH + "/log/RemoveIndividual.log",
# # #     run:
# # #         shell("vcftools --vcf {input.vcf}  --remove {params.removeIndiv} --out {params.outPrefix} --recode > {log} 2>&1")


# # rule Recode:
# #     input:
# #         vcf = IN_PATH + "/NGS/Variant/Sample_jointcall.variant.filt.vcf",
# #     output:
# #         vcf = IN_PATH + "/NGS/Variant/Sample_jointcall.variant.recode.vcf",    
# #     params:
# #         outPrefix = IN_PATH + "/NGS/Variant//Sample_jointcall.variant",
# #     log:
# #         IN_PATH + "/log/RemoveIndividual.log",
# #     run:
# #         shell("vcftools --vcf {input.vcf}  --out {params.outPrefix} --recode > {log} 2>&1")



# rule MAFAllele01GT:
#     input:
#         vcf = IN_PATH + "/NGS/Variant/Sample_jointcall.variant.removeIndiv.recode.vcf",
#     output:
#         vcf = IN_PATH + "/NGS/Variant/Sample_jointcall.variant.removeIndiv.GT.vcf",
#     log:
#         IN_PATH + "/log/MAFAllele01GT.log",
#     run:
#         shell("bcftools annotate -x ^FORMAT/GT {input.vcf} > {output.vcf}")



# rule VariantStats0:
#     input:
#         vcf = IN_PATH + "/NGS/Variant/Sample_jointcall.variant.removeIndiv.GT.vcf",
#     output:
#         stats = IN_PATH + "/NGS/Variant/Sample_jointcall.variant.stats.xls",
#         heter = IN_PATH + "/NGS/Variant/Sample_jointcall.variant.heterzygoisty.xls",
#     params:
#         VariantStats = SCRIPT_DIR + "/VariantStats.py",
#     run:
#         shell("python {params.VariantStats} --vcf {input.vcf} --stats {output.stats} --heter {output.heter}")

# #######################################################################


# ###################### MAF 0.01  LD decay #############################
# rule MAFAllele01:
#     input:
#         vcf = IN_PATH + "/NGS/Variant/Sample_jointcall.variant.removeIndiv.recode.vcf",
#     output:
#         vcf = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf001.recode.vcf",
#         gt = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf001.GT.vcf",
#     params:
#         outPrefix = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf001",
#     log:
#         IN_PATH + "/log/MAFAllele01.log",
#     run:
#         shell("vcftools --vcf {input.vcf} --min-alleles 2 --max-alleles 2  --maf 0.01 --max-missing 0.2 --out {params.outPrefix} --recode > {log} 2>&1")
#         shell("bcftools annotate -x ^FORMAT/GT {output.vcf} > {output.gt}")



rule LDAll:
    input:
        vcf1 = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf001.GT.vcf",
    output:
        ld1 = IN_PATH + "/NGS/Variant/MAF001/all/LD_{Chr}.ld.gz", 
    params:
        q1 = IN_PATH + "/NGS/Variant/MAF001/all/LD_{Chr}",
    run:
        shell("plink --vcf {input.vcf1} --double-id --allow-extra-chr --set-missing-var-ids @:# --maf 0.05 --geno 0.1 --mind 0.5 --chr {wildcards.Chr} --thin 0.1 -r2 gz --ld-window 100 --ld-window-kb 1000 --ld-window-r2 0 --make-bed --out {params.q1}")




rule MergeAllChr:
    input:
        ld1 = expand(IN_PATH + "/NGS/Variant/MAF001/all/LD_{Chr}.ld.gz", Chr = CHRS),
    output:
        ld1 = IN_PATH + "/NGS/Variant/MAF001/all/All_Chr_merge.ld.gz",    
    params:
        ld1 = IN_PATH + "/NGS/Variant/MAF001/all/All_Chr_merge.ld",
    run:
        LD1 = sorted(input.ld1)
        FLD1 = " ".join(LD1)
        first1 = LD1[0]
        shell("pigz -dc {first1} | sed -n '1p' > {params.ld1} ") 
        shell("pigz -dc {FLD1} | grep -v 'CHR_A' >> {params.ld1} ") 
        shell("pigz {params.ld1}")


rule LDStatsAll:
    input:
        ld1 = IN_PATH + "/NGS/Variant/MAF001/all/All_Chr_merge.ld.gz",   
    output:
        q1 = IN_PATH + "/NGS/Variant/MAF001/all/All_Chr_merge.ld_decay_bins",    
    params:
        ld_decay_calc = "/home/wuzhikun/github/speciationgenomics/ld_decay_calc.py",
        q1 = IN_PATH + "/NGS/Variant/MAF001/all/All_Chr_merge",
    run:
        shell("python {params.ld_decay_calc} -i {input.ld1} -o {params.q1}")










rule groupVCF001:
    input:
        vcf = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf001.GT.vcf",
    output:
        vcf1 = IN_PATH + "/NGS/Variant/MAF001/Q1_samples.vcf",
        vcf2 = IN_PATH + "/NGS/Variant/MAF001/Q2_samples.vcf",
        vcf3 = IN_PATH + "/NGS/Variant/MAF001/Q3_samples.vcf",
        vcf4 = IN_PATH + "/NGS/Variant/MAF001/Q4_samples.vcf",
    threads:
        THREADS
    params:
        q1 = IN_PATH + "/NGS/Variant/MAF001/Q1_sample.txt",
        q2 = IN_PATH + "/NGS/Variant/MAF001/Q2_sample.txt",
        q3 = IN_PATH + "/NGS/Variant/MAF001/Q3_sample.txt",
        q4 = IN_PATH + "/NGS/Variant/MAF001/Q4_sample.txt",
    run:
        shell("bcftools view --threads {threads} --samples-file {params.q1} -o {output.vcf1} {input.vcf} ")
        shell("bcftools view --threads {threads} --samples-file {params.q2} -o {output.vcf2} {input.vcf} ")
        shell("bcftools view --threads {threads} --samples-file {params.q3} -o {output.vcf3} {input.vcf} ")
        shell("bcftools view --threads {threads} --samples-file {params.q4} -o {output.vcf4} {input.vcf} ")

rule LD:
    input:
        vcf1 = IN_PATH + "/NGS/Variant/MAF001/Q1_samples.vcf",
        vcf2 = IN_PATH + "/NGS/Variant/MAF001/Q2_samples.vcf",
        vcf3 = IN_PATH + "/NGS/Variant/MAF001/Q3_samples.vcf",
        vcf4 = IN_PATH + "/NGS/Variant/MAF001/Q4_samples.vcf",
    output:
        ld1 = IN_PATH + "/NGS/Variant/MAF001/Q1/LD_{Chr}.ld.gz", 
        ld2 = IN_PATH + "/NGS/Variant/MAF001/Q2/LD_{Chr}.ld.gz", 
        ld3 = IN_PATH + "/NGS/Variant/MAF001/Q3/LD_{Chr}.ld.gz", 
        ld4 = IN_PATH + "/NGS/Variant/MAF001/Q4/LD_{Chr}.ld.gz", 
    params:
        q1 = IN_PATH + "/NGS/Variant/MAF001/Q1/LD_{Chr}",
        q2 = IN_PATH + "/NGS/Variant/MAF001/Q2/LD_{Chr}",
        q3 = IN_PATH + "/NGS/Variant/MAF001/Q3/LD_{Chr}",
        q4 = IN_PATH + "/NGS/Variant/MAF001/Q4/LD_{Chr}",
    run:
        shell("plink --vcf {input.vcf1} --double-id --allow-extra-chr --set-missing-var-ids @:# --maf 0.05 --geno 0.1 --mind 0.5 --chr {wildcards.Chr} --thin 0.1 -r2 gz --ld-window 100 --ld-window-kb 1000 --ld-window-r2 0 --make-bed --out {params.q1}")
        shell("plink --vcf {input.vcf2} --double-id --allow-extra-chr --set-missing-var-ids @:# --maf 0.05 --geno 0.1 --mind 0.5 --chr {wildcards.Chr} --thin 0.1 -r2 gz --ld-window 100 --ld-window-kb 1000 --ld-window-r2 0 --make-bed --out {params.q2}")
        shell("plink --vcf {input.vcf3} --double-id --allow-extra-chr --set-missing-var-ids @:# --maf 0.05 --geno 0.1 --mind 0.5 --chr {wildcards.Chr} --thin 0.1 -r2 gz --ld-window 100 --ld-window-kb 1000 --ld-window-r2 0 --make-bed --out {params.q3}")
        shell("plink --vcf {input.vcf4} --double-id --allow-extra-chr --set-missing-var-ids @:# --maf 0.05 --geno 0.1 --mind 0.5 --chr {wildcards.Chr} --thin 0.1 -r2 gz --ld-window 100 --ld-window-kb 1000 --ld-window-r2 0 --make-bed --out {params.q4}")


rule MergeChr:
    input:
        ld1 = expand(IN_PATH + "/NGS/Variant/MAF001/Q1/LD_{Chr}.ld.gz", Chr = CHRS),
        ld2 = expand(IN_PATH + "/NGS/Variant/MAF001/Q2/LD_{Chr}.ld.gz",  Chr = CHRS),
        ld3 = expand(IN_PATH + "/NGS/Variant/MAF001/Q3/LD_{Chr}.ld.gz",  Chr = CHRS),
        ld4 = expand(IN_PATH + "/NGS/Variant/MAF001/Q4/LD_{Chr}.ld.gz",  Chr = CHRS),
    output:
        ld1 = IN_PATH + "/NGS/Variant/MAF001/Q1/All_Chr_merge.ld.gz",    
        ld2 = IN_PATH + "/NGS/Variant/MAF001/Q2/All_Chr_merge.ld.gz",    
        ld3 = IN_PATH + "/NGS/Variant/MAF001/Q3/All_Chr_merge.ld.gz",    
        ld4 = IN_PATH + "/NGS/Variant/MAF001/Q4/All_Chr_merge.ld.gz",  
    params:
        ld1 = IN_PATH + "/NGS/Variant/MAF001/Q1/All_Chr_merge.ld",
        ld2 = IN_PATH + "/NGS/Variant/MAF001/Q2/All_Chr_merge.ld",
        ld3 = IN_PATH + "/NGS/Variant/MAF001/Q3/All_Chr_merge.ld",
        ld4 = IN_PATH + "/NGS/Variant/MAF001/Q4/All_Chr_merge.ld",
    run:
        LD1 = sorted(input.ld1)
        FLD1 = " ".join(LD1)
        first1 = LD1[0]
        shell("pigz -dc {first1} | sed -n '1p' > {params.ld1} ") 
        shell("pigz -dc {FLD1} | grep -v 'CHR_A' >> {params.ld1} ") 
        shell("pigz {params.ld1}")
        LD2 = sorted(input.ld2)
        FLD2 = " ".join(LD2)
        first2 = LD2[0]
        shell("pigz -dc {first2} | sed -n '1p' > {params.ld2} ") 
        shell("pigz -dc {FLD2} | grep -v 'CHR_A' >> {params.ld2} ") 
        shell("pigz {params.ld2}")
        LD3 = sorted(input.ld3)
        FLD3 = " ".join(LD3)
        first3 = LD3[0]
        shell("pigz -dc {first3} | sed -n '1p' > {params.ld3} ") 
        shell("pigz -dc {FLD3} | grep -v 'CHR_A' >> {params.ld3} ") 
        shell("pigz {params.ld3}")
        LD4 = sorted(input.ld4)
        FLD4 = " ".join(LD4)
        first4 = LD4[0]
        shell("pigz -dc {first4} | sed -n '1p' > {params.ld4} ") 
        shell("pigz -dc {FLD4} | grep -v 'CHR_A' >> {params.ld4} ") 
        shell("pigz {params.ld4}")


### https://github.com/speciationgenomics/scripts/tree/master
### python2 (Assembly)
rule LDStats1:
    input:
        ld1 = IN_PATH + "/NGS/Variant/MAF001/Q1/All_Chr_merge.ld.gz",   
        ld2 = IN_PATH + "/NGS/Variant/MAF001/Q2/All_Chr_merge.ld.gz",    
        ld3 = IN_PATH + "/NGS/Variant/MAF001/Q3/All_Chr_merge.ld.gz",    
        ld4 = IN_PATH + "/NGS/Variant/MAF001/Q4/All_Chr_merge.ld.gz",  
    output:
        q1 = IN_PATH + "/NGS/Variant/MAF001/Q1/All_Chr_merge.ld_decay_bins",  
        q2 = IN_PATH + "/NGS/Variant/MAF001/Q2/All_Chr_merge.ld_decay_bins",   
        q3 = IN_PATH + "/NGS/Variant/MAF001/Q3/All_Chr_merge.ld_decay_bins",     
        q4 = IN_PATH + "/NGS/Variant/MAF001/Q4/All_Chr_merge.ld_decay_bins",   
    params:
        ld_decay_calc = "/home/wuzhikun/github/speciationgenomics/ld_decay_calc.py",
        q1 = IN_PATH + "/NGS/Variant/MAF001/Q1/All_Chr_merge",
        q2 = IN_PATH + "/NGS/Variant/MAF001/Q2/All_Chr_merge",
        q3 = IN_PATH + "/NGS/Variant/MAF001/Q3/All_Chr_merge",
        q4 = IN_PATH + "/NGS/Variant/MAF001/Q4/All_Chr_merge",
    run:
        shell("python {params.ld_decay_calc} -i {input.ld1} -o {params.q1}")
        shell("python {params.ld_decay_calc} -i {input.ld2} -o {params.q2}")
        shell("python {params.ld_decay_calc} -i {input.ld3} -o {params.q3}")
        shell("python {params.ld_decay_calc} -i {input.ld4} -o {params.q4}")


def Chr_LD_average(ld_files, out_file):
    import collections
    import statistics
    PosLD = collections.defaultdict(lambda: collections.defaultdict(list))
    out_h = open(out_file, "w")
    out_h.write("Group\tDistance\tR2\n")
    files = ld_files.split(",")
    for f in files:
        f = f.strip()
        sample = f.split("/")[-2]
        in_h = open(f, "r")
        header = in_h.readline().strip()
        for line in in_h:
            lines = line.strip().split("\t")
            Chr, Pos, ld = lines[:3]
            ld = float(ld)
            PosLD[sample][Pos].append(ld)
        in_h.close()
    ### merge groups
    for s in PosLD:
        POS = PosLD[s]
        for p in POS:
            lds = POS[p]
            ldAve = statistics.mean(lds)
            ldAve = "%.6f" % ldAve
            out_h.write("%s\t%s\t%s\n" % (s, p, ldAve)) 
    out_h.close()


rule ChrAve:
    input:
        q1 = IN_PATH + "/NGS/Variant/MAF001/Q1/All_Chr_merge.ld_decay_bins",  
        q2 = IN_PATH + "/NGS/Variant/MAF001/Q2/All_Chr_merge.ld_decay_bins",   
        q3 = IN_PATH + "/NGS/Variant/MAF001/Q3/All_Chr_merge.ld_decay_bins",  
        q4 = IN_PATH + "/NGS/Variant/MAF001/Q4/All_Chr_merge.ld_decay_bins",  
        a = IN_PATH + "/NGS/Variant/MAF001/all/All_Chr_merge.ld_decay_bins",
    output:
        q1 = IN_PATH + "/NGS/Variant/MAF001/All_group_all_Chr_merge.ld_decay_bins_all.txt",  
    run:
        files = ",".join([input.a, input.q1, input.q2, input.q3, input.q4])
        Chr_LD_average(files, output.q1)


rule LDPlot:
    input:
        q1 = IN_PATH + "/NGS/Variant/MAF001/All_group_all_Chr_merge.ld_decay_bins_all.txt",  
    output:
        q1 = IN_PATH + "/NGS/Variant/MAF001/All_group_all_Chr_merge.ld_decay_bins_all.pdf",
    params:
        LDdecay = SCRIPT_DIR + "/LDdecay.R",
    run:
        shell("Rscript {params.LDdecay} --input {input.q1}  --pdf {output.q1} --width 5 --height 4")



##############################################################



##################### MAF 0.05 ################################



rule MAFAllele:
    input:
        vcf = IN_PATH + "/NGS/Variant/Sample_jointcall.variant.removeIndiv.recode.vcf",
        # vcf = IN_PATH + "/NGS/Variant/Sample_jointcall.variant.recode.vcf",   
    output:
        vcf = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.recode.vcf",
        gt = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.GT.vcf",
    params:
        outPrefix = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005",
    log:
        IN_PATH + "/log/MAFAllele.log",
    run:
        shell("vcftools --vcf {input.vcf} --min-alleles 2 --max-alleles 2  --maf 0.05 --max-missing 0.2 --out {params.outPrefix} --recode > {log} 2>&1")
        shell("bcftools annotate -x ^FORMAT/GT {output.vcf} > {output.gt}")



rule FiltMiss:
    input:
        vcf = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.GT.vcf",
    output:
        vcf = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.filtmiss.recode.vcf",
    params:
        outPrefix = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.filtmiss",
    log:
        IN_PATH + "/log/FiltMiss.log",
    run:
        shell("vcftools --vcf {input.vcf}  --max-missing 0.2 --out {params.outPrefix} --recode > {log} 2>&1")



rule VariantStats:
    input:
        vcf = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.filtmiss.recode.vcf",
    output:
        stats = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.filtmiss.stats.xls",
        heter = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.filtmiss.heterzygoisty.xls",
    params:
        VariantStats = SCRIPT_DIR + "/VariantStats.py",
    run:
        shell("python {params.VariantStats} --vcf {input.vcf} --stats {output.stats} --heter {output.heter}")

##############################################################




rule snpEff:
    input:
        vcf = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.filtmiss.recode.vcf",
    output:
        vcf = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.filtmiss.recode.annotation.vcf",
    run:
        shell("/home/wuzhikun/anaconda3/envs/PanSV/bin/snpEff  ann G98 {input.vcf} > {output.vcf}")


rule snpEffStats:
    input:
        vcf = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.filtmiss.recode.annotation.vcf",
    output:
        vcf = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.filtmiss.recode.annotation.stats.txt",
        cut = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.filtmiss.recode.annotation.cut.txt",
    params:
        snpEffStats = SCRIPT_DIR + "/snpEffStats.py",
    run:
        shell("python {params.snpEffStats} --annotation {input.vcf} --out {output.vcf}")
        shell("grep -v '^##' {input.vcf} | cut -f 1-8 > {output.cut}")


############# Pi ##############

rule PopPi:
    input:
        vcf = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.filtmiss.recode.vcf",
    output:
        q1 = IN_PATH + "/NGS/Variant/Group/Pi/population.windowed.pi",
    threads:
        THREADS
    params:
        Q1 = IN_PATH + "/NGS/Variant/Group/Pi/population",
    log:
        q1 = IN_PATH + "/log/PopPi.log",

    run:
        shell("vcftools --vcf {input.vcf}  --window-pi  100000 --window-pi-step 0 --out {params.Q1} > {log.q1} 2>&1")

###################################




################# Distribution  ##########
rule dist:
    input:
        vcf = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.filtmiss.recode.vcf",
    output:
        pos = IN_PATH + "/NGS/Variant/Sample_jointcall.maf005_pos.txt",
        pdf = IN_PATH + "/NGS/Variant/Sample_jointcall.maf005_pos_distribution.pdf",
    params:
        SNVDistribution = SCRIPT_DIR + "/SNVDistribution.R"
    run:
        shell("grep -v '^#' {input.vcf}| cut -f 1-2 > {output.pos} ")
        shell("Rscript {params.SNVDistribution} --input {output.pos} --pdf {output.pdf}  --maxNum 600")

#########################################




########################## Tree #################
rule distance:
    input:
        vcf = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.filtmiss.recode.vcf",
    output:
        mat = IN_PATH + "/NGS/Variant/tree/Sample_distance.mat",
    run:
        shell("/home/wuzhikun/software/VCF2Dis-1.50/bin/VCF2Dis -InPut {input.vcf} -OutPut {output.mat}")

# rule phylip:
#     input:
#         mat = IN_PATH + "/NGS/Variant/tree/Sample_distance.mat",
#     output:
#         tre = IN_PATH + "/NGS/Variant/tree/Sample_distance.tre",
#     run:
#         shell("/home/wuzhikun/software/phylip-3.697/exe/neighbor {input.mat}")
#################################################



############################ population stats ################

rule groupVCF:
    input:
        vcf = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.filtmiss.recode.vcf",
    output:
        vcf1 = IN_PATH + "/NGS/Variant/Group/Q1_samples.vcf",
        vcf2 = IN_PATH + "/NGS/Variant/Group/Q2_samples.vcf",
        vcf3 = IN_PATH + "/NGS/Variant/Group/Q3_samples.vcf",
    threads:
        THREADS
    params:
        q1 = IN_PATH + "/Q1_sample.txt",
        q2 = IN_PATH + "/Q2_sample.txt",
        q3 = IN_PATH + "/Q3_sample.txt",
    run:
        shell("bcftools view --threads {threads} --samples-file {params.q1} -o {output.vcf1} {input.vcf} ")
        shell("bcftools view --threads {threads} --samples-file {params.q2} -o {output.vcf2} {input.vcf} ")
        shell("bcftools view --threads {threads} --samples-file {params.q3} -o {output.vcf3} {input.vcf} ")


rule pi:
    input:
        vcf1 = IN_PATH + "/NGS/Variant/Group/Q1_samples.vcf",
        vcf2 = IN_PATH + "/NGS/Variant/Group/Q2_samples.vcf",
        vcf3 = IN_PATH + "/NGS/Variant/Group/Q3_samples.vcf",
    output:
        q1 = IN_PATH + "/NGS/Variant/Group/Pi/Q1.windowed.pi",
        q2 = IN_PATH + "/NGS/Variant/Group/Pi/Q2.windowed.pi",
        q3 = IN_PATH + "/NGS/Variant/Group/Pi/Q3.windowed.pi",
    params:
        q1 = IN_PATH + "/NGS/Variant/Group/Pi/Q1",
        q2 = IN_PATH + "/NGS/Variant/Group/Pi/Q2",
        q3 = IN_PATH + "/NGS/Variant/Group/Pi/Q3",
    log:
        q1 = IN_PATH + "/log/Q1_pi.log",
        q2 = IN_PATH + "/log/Q2_pi.log",
        q3 = IN_PATH + "/log/Q3_pi.log",
    run:
        shell("vcftools --vcf {input.vcf1} --out {params.q1}  --window-pi 100000 --window-pi-step 20000 > {log.q1} 2>&1")
        shell("vcftools --vcf {input.vcf2} --out {params.q2}  --window-pi 100000 --window-pi-step 20000 > {log.q2} 2>&1")
        shell("vcftools --vcf {input.vcf3} --out {params.q3}  --window-pi 100000 --window-pi-step 20000 > {log.q3} 2>&1")


rule piSite:
    input:
        vcf1 = IN_PATH + "/NGS/Variant/Group/Q1_samples.vcf",
        vcf2 = IN_PATH + "/NGS/Variant/Group/Q2_samples.vcf",
        vcf3 = IN_PATH + "/NGS/Variant/Group/Q3_samples.vcf",
    output:
        q1 = IN_PATH + "/NGS/Variant/Group/Pi/Q1.sites.pi",
        q2 = IN_PATH + "/NGS/Variant/Group/Pi/Q2.sites.pi",
        q3 = IN_PATH + "/NGS/Variant/Group/Pi/Q3.sites.pi",
    params:
        q1 = IN_PATH + "/NGS/Variant/Group/Pi/Q1",
        q2 = IN_PATH + "/NGS/Variant/Group/Pi/Q2",
        q3 = IN_PATH + "/NGS/Variant/Group/Pi/Q3",
    log:
        q1 = IN_PATH + "/log/Q1_pi_site.log",
        q2 = IN_PATH + "/log/Q2_pi_site.log",
        q3 = IN_PATH + "/log/Q3_pi_site.log",
    run:
        shell("vcftools --vcf {input.vcf1} --out {params.q1}  --site-pi > {log.q1} 2>&1")
        shell("vcftools --vcf {input.vcf2} --out {params.q2}  --site-pi > {log.q2} 2>&1")
        shell("vcftools --vcf {input.vcf3} --out {params.q3}  --site-pi > {log.q3} 2>&1")


rule pistats:
    input:
        q1 = IN_PATH + "/NGS/Variant/Group/Pi/Q1.windowed.pi",
        q2 = IN_PATH + "/NGS/Variant/Group/Pi/Q2.windowed.pi",
        q3 = IN_PATH + "/NGS/Variant/Group/Pi/Q3.windowed.pi",
    output:
        stat = IN_PATH + "/NGS/Variant/Group/Pi/groups.sites.pi.stats.txt",
    run:
        cmd1 = """ awk '{ sum += $5; n++ } END {  print sum / n; }' %s > %s """ % (input.q1, output.stat)
        os.system(cmd1)
        cmd2 = """ awk '{ sum += $5; n++ } END {  print sum / n; }' %s >> %s """ % (input.q2, output.stat)
        os.system(cmd2)
        cmd3 = """ awk '{ sum += $5; n++ } END {  print sum / n; }' %s >> %s """ % (input.q3, output.stat)
        os.system(cmd3)


rule Fst:
    input:
        vcf = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.filtmiss.recode.vcf",
    output:
        q1 = IN_PATH + "/NGS/Variant/Group/Fst/Q1_Q2.windowed.weir.fst",
        q2 = IN_PATH + "/NGS/Variant/Group/Fst/Q1_Q3.windowed.weir.fst",
        q3 = IN_PATH + "/NGS/Variant/Group/Fst/Q2_Q3.windowed.weir.fst",
    threads:
        THREADS
    params:
        q1 = IN_PATH + "/Q1_sample.txt",
        q2 = IN_PATH + "/Q2_sample.txt",
        q3 = IN_PATH + "/Q3_sample.txt",
        Q1 = IN_PATH + "/NGS/Variant/Group/Fst/Q1_Q2",
        Q2 = IN_PATH + "/NGS/Variant/Group/Fst/Q1_Q3",
        Q3 = IN_PATH + "/NGS/Variant/Group/Fst/Q2_Q3",
    log:
        q1 = IN_PATH + "/log/Q1_Q2_fst.log",
        q2 = IN_PATH + "/log/Q1_Q3_fst.log",
        q3 = IN_PATH + "/log/Q2_Q3_fst.log",
    run:
        shell("vcftools --vcf {input.vcf} --weir-fst-pop {params.q1}  --weir-fst-pop {params.q2} --out {params.Q1} --fst-window-size 100000 --fst-window-step 20000 > {log.q1} 2>&1")
        shell("vcftools --vcf {input.vcf} --weir-fst-pop {params.q1}  --weir-fst-pop {params.q3} --out {params.Q2} --fst-window-size 100000 --fst-window-step 20000 > {log.q2} 2>&1")
        shell("vcftools --vcf {input.vcf} --weir-fst-pop {params.q2}  --weir-fst-pop {params.q3} --out {params.Q3} --fst-window-size 100000 --fst-window-step 20000 > {log.q3} 2>&1")

############################################################################



# rule AlleleIdenty:
#     input:
#         gt = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.GT.vcf",
#     output:
#         ratio = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.GT.identity.xls",
#     params:
#         AllelIdentityRatio = SCRIPT_DIR + "/AllelIdentityRatio.py",
#     run:
#         shell("python {params.AllelIdentityRatio} --vcf {input.gt} --out {output.ratio}")



rule ConvertPed:
    input:
        # vcf = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.GT.vcf",
        vcf = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.filtmiss.recode.vcf",
    output:
        Ped = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.ped",
        Map = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.map",
    params:
        outPrefix = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005",
        Map = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005-1.map",
    log:
        IN_PATH + "/log/ConvertPed.log",
    run:
        shell("vcftools --vcf {input.vcf} --plink --out {params.outPrefix} > {log} 2>&1")
        cmd = """ cut -f 2 %s | sed 's/:/\t/' | awk '{print $1"\t"$1":"$2"\t"0"\t"$2}' > %s """ % (output.Map, params.Map)
        os.system(cmd)
        shell("mv {params.Map} {output.Map}")


# rule PruneIN:
#     input:
#         Ped = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.ped",
#         Map = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.map",
#     output:
#         prune = IN_PATH + "/NGS/Variant/Sample_variant_PruneLD.prune.in",
#         bed = IN_PATH + "/NGS/Variant/Sample_variant_PruneLD.bed",
#     params:
#         outPrefix = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005",
#         prunePrefix = IN_PATH + "/NGS/Variant/Sample_variant_PruneLD",
#     log:
#         IN_PATH + "/log/PruneIN.log",
#     run:
#         shell("plink --file {params.outPrefix} --indep-pairwise 50 10 0.5 --out {params.prunePrefix} --allow-extra-chr  > {log} 2>&1")
#         shell("plink --file {params.outPrefix} --extract {output.prune} -make-bed  --out {params.prunePrefix} --allow-extra-chr >> {log} 2>&1")
# ####################################################################

# ####################### PCA ##########################
# rule PCA:
#     input:
#         bed = IN_PATH + "/NGS/Variant/Sample_variant_PruneLD.bed",
#     output:
#         eigenval = IN_PATH + "/NGS/Variant/PCA/Samples_PCA.eigenval",
#         eigenvec = IN_PATH + "/NGS/Variant/PCA/Samples_PCA.eigenvec",
#     params:
#         prunePrefix = IN_PATH + "/NGS/Variant/Sample_variant_PruneLD",
#         outPrefix = IN_PATH + "/NGS/Variant/PCA/Samples_PCA",
#     log:
#         IN_PATH + "/log/PCA.log",
#     run:
#         shell("plink --bfile {params.prunePrefix} --pca --out {params.outPrefix} --allow-extra-chr > {log} 2>&1")


# rule PCA1:
#     input:
#         pca12 = IN_PATH + "/NGS/Variant/PCA/Sample_PCA12.txt",
#         pca34 = IN_PATH + "/NGS/Variant/PCA/Sample_PCA34.txt",
#     output:
#         pca12Group = IN_PATH + "/NGS/Variant/PCA/Sample_group_PCA12.txt",
#         pca34Group = IN_PATH + "/NGS/Variant/PCA/Sample_group_PCA34.txt",
#     params:
#         PCAaddGroup = SCRIPT_DIR + "/PCAaddGroup.py",
#         metafile = IN_PATH + "/meta_infor.txt",
#     log:
#         IN_PATH + "/log/PCA1.log", 
#     run:
#         shell("python {params.PCAaddGroup} --structure {params.metafile} --pca {input.pca12} --out {output.pca12Group} --column Type > {log} 2>&1")
#         shell("python {params.PCAaddGroup} --structure {params.metafile} --pca {input.pca34} --out {output.pca34Group} --column Type >> {log} 2>&1")



# rule PCAPlot:
#     input:
#         pca12Group = IN_PATH + "/NGS/Variant/PCA/Sample_group_PCA12.txt",
#         pca34Group = IN_PATH + "/NGS/Variant/PCA/Sample_group_PCA34.txt",
#     output:
#         pdf12 = IN_PATH + "/NGS/Variant/PCA/Sample_group_PCA12.pdf",
#         pdf34 = IN_PATH + "/NGS/Variant/PCA/Sample_group_PCA34.pdf",
#     params:
#         genoPCA = SCRIPT_DIR + "/genoPCA.R",
#         width = 5,
#         height = 4,
#     log:
#         IN_PATH + "/log/PCAPlot.log", 
#     run:
#         cmd2 = "source activate Rmeta && Rscript %s --input %s --pdf %s --width %s --height %s > %s 2>&1" % (params.genoPCA, input.pca12Group, output.pdf12, params.width, params.height, log)
#         print(cmd2)
#         os.system(cmd2)
#         cmd2 = "source activate Rmeta && Rscript %s --input %s --pdf %s --width %s --height %s > %s 2>&1" % (params.genoPCA, input.pca34Group, output.pdf34, params.width, params.height, log)
#         print(cmd2)
#         os.system(cmd2)

# ###################################################

# ##################### structure ##########################
# rule Structure:
#     input:
#         bed = IN_PATH + "/NGS/Variant/Sample_variant_PruneLD.bed",
#     output:
#         p = IN_PATH + "/NGS/Variant/Structure/struct.{popu}.meanP",
#     params:
#         prunePrefix = IN_PATH + "/NGS/Variant/Sample_variant_PruneLD",
#         outPrefix = IN_PATH + "/NGS/Variant/Structure/struct",
#     log:
#         IN_PATH + "/log/Structure_{popu}.log",
#     run:
#         cmd = "source activate Assembly && python /home/wuzhikun/anaconda3/envs/Assembly/bin/structure.py -K %s --input %s --out %s --full --seed 100 " % (wildcards.popu, params.prunePrefix, params.outPrefix)
#         os.system(cmd)
    

# rule chooseK:
#     input:
#         p = expand(IN_PATH + "/NGS/Variant/Structure/struct.{popu}.meanP", popu=POPS),
#     output:
#         chooseK = IN_PATH + "/NGS/Variant/Structure/chooseK.txt",
#     params:
#         outPrefix = IN_PATH + "/NGS/Variant/Structure/struct",
#     run:
#         cmd = "source activate Assembly && python /home/wuzhikun/anaconda3/envs/Assembly/bin/chooseK.py --input %s > %s" % (params.outPrefix, output.chooseK)
#         os.system(cmd)


### DISPLAY
# rule PlotStructe:
#     input:
#         chooseK = IN_PATH + "/NGS/Variant/Structure/chooseK.txt",
#         profile = IN_PATH + "/NGS/Variant/Structure/profile.txt",
#     output:
#         chooseK = IN_PATH + "/NGS/Variant/Structure/chooseK3.svg",
#     params:
#         outPrefix = IN_PATH + "/NGS/Variant/Structure/struct",
#     run:
#         cmd = "source activate Assembly && python /home/wuzhikun/anaconda3/envs/Assembly/bin/distruct.py -K 3 --input %s --output %s --popfile %s " % (params.outPrefix, output.chooseK, input.profile)
#         print(cmd)
#         os.system(cmd)
####################################################



########################## admixture ###############
### https://zhuanlan.zhihu.com/p/370444841



# rule PruneChr:
#     input:
#         Ped = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.ped",
#         Map = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.map",
#     output:
#         prune1 = IN_PATH + "/NGS/Variant/Sample_variant_PruneLD_r01.prune.in",
#         prune2 = IN_PATH + "/NGS/Variant/Sample_variant_PruneLD_r01.prune.in.Chr",
#         bed = IN_PATH + "/NGS/Variant/Sample_variant_PruneLD_r01_Chr.bed",
#     params:
#         outPrefix = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005",
#         prunePrefix = IN_PATH + "/NGS/Variant/Sample_variant_PruneLD_r01",
#         outPrefix2 = IN_PATH + "/NGS/Variant/Sample_variant_PruneLD_r01_Chr",
#     log:
#         IN_PATH + "/log/PruneIN.log",
#     run:
#         shell("plink --file {params.outPrefix} --indep-pairwise 50 10 0.1 --out {params.prunePrefix} --allow-extra-chr  > {log} 2>&1")
#         shell("grep  '^Chr' {output.prune1} > {output.prune2}")
#         shell("plink --file {params.outPrefix} --extract {output.prune2} -make-bed  --out {params.outPrefix2} --allow-extra-chr >> {log} 2>&1")


# rule admixture:
#     input:
#         bed = IN_PATH + "/NGS/Variant/Sample_variant_PruneLD_r01_Chr.bed",
#     output:
#         log = IN_PATH + "/NGS/Variant/admixture_r01/log_{popu}.out",
#     threads:
#         THREADS
#     run:
#         shell("admixture --cv {input.bed}  {wildcards.popu}  -j{threads} |tee {output.log}")

# rule admixtureMerge:
#     input:
#         log = expand(IN_PATH + "/NGS/Variant/admixture_r01/log_{popu}.out", popu=POPS),
#     output:
#         valid = IN_PATH + "/NGS/Variant/admixture_r01/admixture_cross_validation_error.txt", 
#     run:
#         Logs = " ".join(input.log)
#         shell("grep -h CV {Logs} | sort -nk4 -t ' ' > {output.valid} ")

# #################################################################################



# #################### ibd #################
# rule ibd:
#     input:
#         bed = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.bed",
#         fam = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.fam",
#     output:
#         gz = IN_PATH + "/NGS/Variant/ibd/Sample_ibd.grm.gz",
#         grm = IN_PATH + "/NGS/Variant/ibd/Sample_ibd.grm.txt",
#     params:
#         inPrefix = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005",
#         outPrefix = IN_PATH + "/NGS/Variant/ibd/Sample_ibd",
#     log:
#         IN_PATH + "/log/ibd.log",
#     run:
#         shell("plink  -bfile {params.inPrefix} --make-grm-gz  --out {params.outPrefix} --allow-extra-chr > {log} 2>&1")
#         shell("pigz -dc {output.gz} > {output.grm}")

# ###########################################



# #################### GWAS ##################
# rule makebed:
#     input:
#         Ped = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.ped",
#         Map = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.map",
#     output:
#         bed = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.bed",
#         fam = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.fam",
#     params:
#         outPrefix = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005",
#     log:
#         IN_PATH + "/log/makebed.log",
#     run:
#         shell("plink --file {params.outPrefix} --make-bed --out {params.outPrefix} --allow-extra-chr  > {log} 2>&1")



# rule modifyfam:
#     input:
#         fam = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.fam",
#     output:
#         fam = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005-1.fam",
#     run:
#         cmd = """awk '{print $1" "$2" "$3" "$4" "$5" "1}' %s > %s """ % (input.fam, output.fam)
#         os.system(cmd)
#         shell("cp {output.fam} {input.fam}")
        


# rule kinship:
#     input:
#         bed = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.bed",
#         fam = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005-1.fam",
#     output:
#         ks = IN_PATH + "/NGS/Variant/GWASGemma/output/geno_kS.cXX.txt",
#     params:
#         inPrefix = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005",
#         outDir = IN_PATH + "/NGS/Variant/GWASGemma",
#         outPrefix = "geno_kS",
#     log:
#         IN_PATH + "/log/kinship.log", 
#     run:
#         if not os.path.exists(params.outDir):
#             os.makedirs(params.outDir)
#         # cmd = "source activate WGS && cd %s && gemma -bfile %s -gk 1 -o %s > %s 2>&1" % (params.outDir, params.inPrefix, params.outPrefix, log)
#         # print(cmd)
#         # os.system(cmd)
#         shell("cd {params.outDir} && gemma -bfile {params.inPrefix} -gk 1 -o {params.outPrefix} > {log} 2>&1")



rule associations:
    input:
        bed = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.bed",
        ks = IN_PATH + "/NGS/Variant/GWASGemma/output/geno_kS.cXX.txt",
        # pheno1 = IN_PATH + "/NGS/phenotype/{trait}.txt",
        pheno1 = IN_PATH + "/NGS/Variant/phenotype/{trait}.txt",
    output:
        asso1 = IN_PATH + "/NGS/Variant/GWASGemma/output/{trait}.assoc.txt",
    threads:
        THREADS
    params:
        inPrefix = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005",
        outDir = IN_PATH + "/NGS/Variant/GWASGemma",
    log:
        IN_PATH + "/log/association_{trait}.log", 
    run:
        shell("cd {params.outDir} && gemma -bfile {params.inPrefix} -k {input.ks} -lmm 1 -o {wildcards.trait} -p {input.pheno1} > {log} 2>&1")


# rule TraitManhattanPlot0:
#     input:
#         assoc = IN_PATH + "/NGS/Variant/GWASGemma/output/{trait}.assoc.txt",
#     output:
#         assoc = IN_PATH + "/NGS/Variant/GWASGemma/output/plot/{trait}.assoc.txt",
#         assoc1 = IN_PATH + "/NGS/Variant/GWASGemma/output/plot/{trait}.assoc-1.txt",
#     run:
#         shell("echo 'SNP CHR BP P' | sed 's/ /\t/g'   > {output.assoc}")
#         cmd = """sed '1d' %s | awk '{print $2"\t"$1"\t"$3"\t"$12}' >> %s"""  % (input.assoc, output.assoc)
#         print(cmd)
#         os.system(cmd)
#         changeChrString(output.assoc, output.assoc1)

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
        assoc = IN_PATH + "/NGS/Variant/GWASGemma/output/{trait}.assoc.txt",
    output:
        assoc = IN_PATH + "/NGS/Variant/GWASGemma/output/plot/{trait}.assoc.txt",
        assoc1 = IN_PATH + "/NGS/Variant/GWASGemma/output/plot/{trait}.assoc-1.txt",
    run:
        shell("echo 'SNP CHR BP P' | sed 's/ /\t/g'   > {output.assoc}")
        cmd = """sed '1d' %s | grep -v "_un" | awk '{print $2"\t"$1"\t"$3"\t"$12}' >> %s"""  % (input.assoc, output.assoc)
        print(cmd)
        os.system(cmd)
        changeChrString(output.assoc, output.assoc1)




rule TraitManhattanPlot1:
    input:
        assoc = IN_PATH + "/NGS/Variant/GWASGemma/output/plot/{trait}.assoc-1.txt",
    output:
        manhattan = IN_PATH + "/NGS/Variant/GWASGemma/output/plot/{trait}.assoc_mahattan.jpeg",
        qq = IN_PATH + "/NGS/Variant/GWASGemma/output/plot/{trait}.assoc_qq.jpeg",
    params:
        GWASPlot = SCRIPT_DIR + "/GWASPlot.R",
    log:
        IN_PATH + "/log/TraitManhattanPlot_{trait}.log", 
    run:
        shell("Rscript {params.GWASPlot} --input {input.assoc} --manhattan {output.manhattan} --qq {output.qq} > {log} 2>&1")




rule SigLoci2:
    input:
        asso = IN_PATH + "/NGS/Variant/GWASGemma/output/plot/{trait}.assoc-1.txt",
    output:
        sig = IN_PATH + "/NGS/Variant/GWASGemma/output/plot/{trait}.assoc_sig.txt",
    run:
        # cmd = "awk '{if ($4 < 2.3e-7){print $0}}'  %s > %s" % (input.asso, output.sig)
        cmd = "awk '{if ($4 < 2.57e-7){print $0}}'  %s > %s" % (input.asso, output.sig)
        print(cmd)
        os.system(cmd)



# rule association:
#     input:
#         bed = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.bed",
#         ks = IN_PATH + "/NGS/Variant/GWASGemma/output/geno_kS.cXX.txt",
#         pheno1 = IN_PATH + "/NGS/phenotype/Glucose_pheno.txt",
#         pheno2 = IN_PATH + "/NGS/phenotype/Levulose_pheno.txt",
#         pheno3 = IN_PATH + "/NGS/phenotype/Sucrose_pheno.txt",
#     output:
#         asso1 = IN_PATH + "/NGS/Variant/GWASGemma/output/Glucose.assoc.txt",
#         asso2 = IN_PATH + "/NGS/Variant/GWASGemma/output/Levulose.assoc.txt",
#         asso3 = IN_PATH + "/NGS/Variant/GWASGemma/output/Sucrose.assoc.txt",
#     threads:
#         THREADS
#     params:
#         inPrefix = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005",
#         outDir = IN_PATH + "/NGS/Variant/GWASGemma",
#     log:
#         IN_PATH + "/log/association.log", 
#     run:
#         ### gemma -bfile /home/wuzhikun/Project/Population/population/plink/Sample_SV_geno -k /home/wuzhikun/Project/Population/population/GWAS/output/SV_geno_kS.cXX.txt -lmm 1 -o pheno1 -p /home/wuzhikun/Project/Population/phenotype/CN329_phenotype_value-1.txt
#         # shell("cut -f {wildcards.trait} {input.pheno} > {output.pheno}")
#         shell("cd {params.outDir} && gemma -bfile {params.inPrefix} -k {input.ks} -lmm 1 -o Glucose -p {input.pheno1} > {log} 2>&1")
#         shell("cd {params.outDir} && gemma -bfile {params.inPrefix} -k {input.ks} -lmm 1 -o Levulose -p {input.pheno2} > {log} 2>&1")
#         shell("cd {params.outDir} && gemma -bfile {params.inPrefix} -k {input.ks} -lmm 1 -o Sucrose -p {input.pheno3} > {log} 2>&1")



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


# rule ManhattanPlot0:
#     input:
#         assoc = IN_PATH + "/NGS/Variant/GWASGemma/output/Glucose.assoc.txt",
#     output:
#         assoc = IN_PATH + "/NGS/Variant/GWASGemma/output/plot/Glucose.assoc.txt",
#         assoc1 = IN_PATH + "/NGS/Variant/GWASGemma/output/plot/Glucose.assoc-1.txt",
#     run:
#         shell("echo 'SNP CHR BP P' | sed 's/ /\t/g'   > {output.assoc}")
#         cmd = """sed '1d' %s | awk '{print $2"\t"$1"\t"$3"\t"$12}' >> %s"""  % (input.assoc, output.assoc)
#         print(cmd)
#         os.system(cmd)
#         changeChrString(output.assoc, output.assoc1)



# rule ManhattanPlot1:
#     input:
#         assoc = IN_PATH + "/NGS/Variant/GWASGemma/output/plot/Glucose.assoc-1.txt",
#     output:
#         manhattan = IN_PATH + "/NGS/Variant/GWASGemma/output/plot/Glucose.assoc_mahattan.jpeg",
#         qq = IN_PATH + "/NGS/Variant/GWASGemma/output/plot/Glucose.assoc_qq.jpeg",
#     params:
#         GWASPlot = SCRIPT_DIR + "/GWASPlot.R",
#     log:
#         IN_PATH + "/log/ManhattanPlot.log", 
#     run:
#         shell("Rscript {params.GWASPlot} --input {input.assoc} --manhattan {output.manhattan} --qq {output.qq} > {log} 2>&1")




# rule SigLoci2:
#     input:
#         asso = IN_PATH + "/NGS/Variant/GWASGemma/output/plot/Glucose.assoc-1.txt",
#     output:
#         sig = IN_PATH + "/NGS/Variant/GWASGemma/output/plot/Glucose.assoc_sig.txt",
#     run:
#         # cmd = "awk '{if ($12< 5e-8){print $0}}'  %s > %s" % (input.asso, output.sig)
#         cmd = "awk '{if ($4 < 3.27e-7){print $0}}'  %s > %s" % (input.asso, output.sig)
#         print(cmd)
#         os.system(cmd)




# >>> 0.05/3054377
# 1.6369950402324273e-08
# >>> 1/3054377
# 3.273990080464854e-07

rule ManhattanPlot01:
    input:
        assoc = IN_PATH + "/NGS/Variant/GWASGemma/output/Sucrose.assoc.txt",
    output:
        assoc = IN_PATH + "/NGS/Variant/GWASGemma/output/plot/Sucrose.assoc.txt",
        assoc1 = IN_PATH + "/NGS/Variant/GWASGemma/output/plot/Sucrose.assoc-1.txt",
    run:
        shell("echo 'SNP CHR BP P' | sed 's/ /\t/g'   > {output.assoc}")
        cmd = """sed '1d' %s | awk '{print $2"\t"$1"\t"$3"\t"$12}' >> %s"""  % (input.assoc, output.assoc)
        print(cmd)
        os.system(cmd)
        changeChrString(output.assoc, output.assoc1)

rule ManhattanPlot11:
    input:
        assoc = IN_PATH + "/NGS/Variant/GWASGemma/output/plot/Sucrose.assoc-1.txt",
    output:
        manhattan = IN_PATH + "/NGS/Variant/GWASGemma/output/plot/Sucrose.assoc_mahattan.jpeg",
        qq = IN_PATH + "/NGS/Variant/GWASGemma/output/plot/Sucrose.assoc_qq.jpeg",
    params:
        GWASPlot = SCRIPT_DIR + "/GWASPlot.R",
    log:
        IN_PATH + "/log/ManhattanPlot.log", 
    run:
        shell("Rscript {params.GWASPlot} --input {input.assoc} --manhattan {output.manhattan} --qq {output.qq} > {log} 2>&1")






rule ManhattanPlot02:
    input:
        assoc = IN_PATH + "/NGS/Variant/GWASGemma/output/Levulose.assoc.txt",
    output:
        assoc = IN_PATH + "/NGS/Variant/GWASGemma/output/plot/Levulose.assoc.txt",
        assoc1 = IN_PATH + "/NGS/Variant/GWASGemma/output/plot/Levulose.assoc-1.txt",
    run:
        shell("echo 'SNP CHR BP P' | sed 's/ /\t/g'   > {output.assoc}")
        cmd = """sed '1d' %s | awk '{print $2"\t"$1"\t"$3"\t"$12}' >> %s"""  % (input.assoc, output.assoc)
        print(cmd)
        os.system(cmd)
        changeChrString(output.assoc, output.assoc1)

rule ManhattanPlot12:
    input:
        assoc = IN_PATH + "/NGS/Variant/GWASGemma/output/plot/Levulose.assoc-1.txt",
    output:
        manhattan = IN_PATH + "/NGS/Variant/GWASGemma/output/plot/Levulose.assoc_mahattan.jpeg",
        qq = IN_PATH + "/NGS/Variant/GWASGemma/output/plot/Levulose.assoc_qq.jpeg",
    params:
        GWASPlot = SCRIPT_DIR + "/GWASPlot.R",
    log:
        IN_PATH + "/log/ManhattanPlot.log", 
    run:
        shell("Rscript {params.GWASPlot} --input {input.assoc} --manhattan {output.manhattan} --qq {output.qq} > {log} 2>&1")



############################################
