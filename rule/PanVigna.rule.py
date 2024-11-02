rule filtVariant:
    input:
        vcf = IN_PATH + "/NGS/Joint/VCF/Chrom_Fengchan6_joint_variant.vcf.gz",
    output:
        vcf = IN_PATH + "/NGS/Population/Chr_joint_variant.recode.vcf",
    params:
        outPrefix = IN_PATH + "/NGS/Population/Chr_joint_variant",
    log:
        IN_PATH + "/log/filtVariant.log"
    run:
        shell("vcftools --gzvcf {input.vcf} --max-missing 0.1 --min-alleles 2  --max-alleles 2 --min-meanDP 10 --recode  --out {params.outPrefix} > {log}")


rule vcfStat:
    input:
        vcf = IN_PATH + "/NGS/Joint/VCF/Chrom_Fengchan6_joint_variant.vcf.gz",
    output:
        stat = IN_PATH + "/NGS/Joint/VCF/Chrom_Fengchan6_joint_variant.stat.txt",
    run:
        shell("rtg vcfstats {input.vcf} > {output.stat}")


rule vcfStat2:
    input:
        vcf = IN_PATH + "/NGS/Population/Chr_joint_variant.recode.vcf",
    output:
        stat = IN_PATH + "/NGS/Population/Chr_joint_variant.stat.txt",
    run:
        shell("rtg vcfstats {input.vcf} > {output.stat}")




rule vcfGT:
    input:
        vcf = IN_PATH + "/NGS/Population/Chr_joint_variant.recode.vcf",
    output:
        vcf = temp(IN_PATH + "/NGS/Population/Chr_joint_variant_GT.vcf"),
    run:
        shell("bcftools annotate -x ^FORMAT/GT {input.vcf} > {output.vcf}")



rule vcfGTgz:
    input:
        vcf = IN_PATH + "/NGS/Population/Chr_joint_variant_GT.vcf",
    output:
        vcf = IN_PATH + "/NGS/Population/Chr_joint_variant_GT.vcf.gz",
    run:
        shell("bgzip {input.vcf}")
        shell("tabix -p vcf {output.vcf}")


rule vcf2bed:
    input:
        vcf = IN_PATH + "/NGS/Population/Chr_joint_variant_GT.vcf.gz",
    output:
        bed = IN_PATH + "/NGS/Population/Chr_joint_variant_GT.bed",
    params:
        outPrefix = IN_PATH + "/NGS/Population/Chr_joint_variant_GT",
    run:
        shell("plink2 --vcf {input.vcf} --make-bed --allow-extra-chr --out {params.outPrefix} ")


rule filtMaf:
    input:
        bed = IN_PATH + "/NGS/Population/Chr_joint_variant_GT.bed",
    output:
        bed = IN_PATH + "/NGS/Population/Chr_joint_variant_maf05.bed",
    params:
        inPrefix = IN_PATH + "/NGS/Population/Chr_joint_variant_GT",
        outPrefix = IN_PATH + "/NGS/Population/Chr_joint_variant_maf05",
    run:
        shell("plink2 --bfile {params.inPrefix} --maf 0.05  --make-bed --allow-extra-chr --out {params.outPrefix} ")




'''
rule SNVprune:
    input:
        vcf = IN_PATH + "/NGS/Population/Chr_joint_variant_GT.bed",
    output:
        prune = IN_PATH + "/NGS/Population/Chr_joint_variant_GT.prune.in",
    params:
        inPrefix = IN_PATH + "/NGS/Population/Chr_joint_variant_GT",
        outPrefix = IN_PATH + "/NGS/Population/Chr_joint_variant_GT.prune",
    run:
        shell("plink2 --bfile {params.inPrefix} --make-bed --allow-extra-chr --indep-pairwise 50 10 0.1  --out {params.outPrefix} ")



rule extractSNV:
    input:
        prune = IN_PATH + "/NGS/Population/Chr_joint_variant_GT.prune.in",
    output:
        bed = IN_PATH + "/NGS/Population/Chr_joint_variant_GT.LDfilt.bed", 
    params:
        inPrefix = IN_PATH + "/NGS/Population/Chr_joint_variant_GT",
        outPrefix = IN_PATH + "/NGS/Population/Chr_joint_variant_GT.LDfilt",
    run:
        shell("plink2 --bfile {params.inPrefix} --extract {input.prune} --make-bed --out {params.outPrefix}")

'''





rule PCA:
    input:
        bed = IN_PATH + "/NGS/Population/Chr_joint_variant_GT.bed",
    output:
        eigenval = IN_PATH + "/NGS/Population/PCA/Samples_PCA.eigenval",
        eigenvec = IN_PATH + "/NGS/Population/PCA/Samples_PCA.eigenvec",
    params:
        inPrefix = IN_PATH + "/NGS/Population/Chr_joint_variant_GT",
        outPrefix = IN_PATH + "/NGS/Population/PCA/Samples_PCA",
    log:
        IN_PATH + "/log/PCA.log",
    run:
        shell("plink --bfile {params.inPrefix} --pca --out {params.outPrefix} --allow-extra-chr > {log} 2>&1")


rule PCA12:
    input:
        eigenvec = IN_PATH + "/NGS/Population/PCA/Samples_PCA.eigenvec", 
    output:
        PCA = IN_PATH + "/NGS/Population/PCA/Samples_PCA.txt",
    run:
        out_h = open(output.PCA, "w")
        out_h.write("FID\tIID\tPC1\tPC2\tPC3\tPC4\tPC5\tPC6\tPC7\tPC8\tPC9\tPC10\n")
        out_h.close()
        shell("sed 's/ /\t/g' {input.eigenvec} >> {output.PCA}")


rule PCAPlot:
    input:
        pca = IN_PATH + "/NGS/Population/PCA/Samples_PCA.txt",
    output:
        pdf12 = IN_PATH + "/NGS/Population/PCA/Sample_group_PCA12.pdf",
        pdf34 = IN_PATH + "/NGS/Population/PCA/Sample_group_PCA34.pdf",
    params:
        plinkGenoPCA = SCRIPT_DIR + "/plinkGenoPCA.R",
        width = 5,
        height = 4,
    log:
        IN_PATH + "/log/PCAPlot.log", 
    run:
        cmd2 = "source activate Rmeta && Rscript %s --input %s --pdf %s --width %s --height %s > %s 2>&1" % (params.plinkGenoPCA, input.pca, output.pdf12, params.width, params.height, log)
        print(cmd2)
        os.system(cmd2)
        cmd2 = "source activate Rmeta && Rscript %s --input %s --pdf %s --width %s --height %s > %s 2>&1" % (params.plinkGenoPCA, input.pca, output.pdf34, params.width, params.height, log)
        print(cmd2)
        os.system(cmd2)

############################ extract sampel ####################


rule extractSample:
    input:
        sample = IN_PATH + "/NGS/Population/sample34.txt",
        bed = IN_PATH + "/NGS/Population/Chr_joint_variant_GT.bed",
    output:
        bed = IN_PATH + "/NGS/Population/Sample34_joint_variant.bed",
    params:
        inPrefix = IN_PATH + "/NGS/Population/Chr_joint_variant_GT",
        outPrefix = IN_PATH + "/NGS/Population/Sample34_joint_variant",
    log:
        IN_PATH + "/log/extractSample.log",
    run:
        shell("plink2 --allow-extra-chr --maf 0.01  --bfile {params.inPrefix} --keep {input.sample} --make-bed --out {params.outPrefix} > {log} 2>&1")





rule PCASample34:
    input:
        IN_PATH + "/NGS/Population/Sample34_joint_variant.bed",
    output:
        eigenval = IN_PATH + "/NGS/Population/PCA_sample34/Samples34_PCA.eigenval",
        eigenvec = IN_PATH + "/NGS/Population/PCA_sample34/Samples34_PCA.eigenvec",
    params:
        inPrefix = IN_PATH + "/NGS/Population/Sample34_joint_variant",
        outPrefix = IN_PATH + "/NGS/Population/PCA_sample34/Samples34_PCA",
    log:
        IN_PATH + "/log/PCASample34.log",
    run:
        shell("plink --bfile {params.inPrefix} --pca --out {params.outPrefix} --allow-extra-chr > {log} 2>&1")


rule Sample34PCA12:
    input:
        eigenvec = IN_PATH + "/NGS/Population/PCA_sample34/Samples34_PCA.eigenvec", 
    output:
        PCA = IN_PATH + "/NGS/Population/PCA_sample34/Samples34_PCA.txt",
    run:
        out_h = open(output.PCA, "w")
        out_h.write("FID\tIID\tPC1\tPC2\tPC3\tPC4\tPC5\tPC6\tPC7\tPC8\tPC9\tPC10\n")
        out_h.close()
        shell("sed 's/ /\t/g' {input.eigenvec} >> {output.PCA}")



rule Sample34PCAPlot:
    input:
        pca = IN_PATH + "/NGS/Population/PCA_sample34/Samples34_PCA.txt",
    output:
        pdf12 = IN_PATH + "/NGS/Population/PCA_sample34/Samples34_group_PCA12.pdf",
        pdf34 = IN_PATH + "/NGS/Population/PCA_sample34/Samples34_group_PCA34.pdf",
    params:
        plinkGenoPCA = SCRIPT_DIR + "/plinkGenoPCA.R",
        width = 5,
        height = 4,
    log:
        IN_PATH + "/log/Sample34PCAPlot.log", 
    run:
        cmd2 = "source activate Rmeta && Rscript %s --input %s --pdf %s --width %s --height %s > %s 2>&1" % (params.plinkGenoPCA, input.pca, output.pdf12, params.width, params.height, log)
        print(cmd2)
        os.system(cmd2)
        cmd2 = "source activate Rmeta && Rscript %s --input %s --pdf %s --width %s --height %s > %s 2>&1" % (params.plinkGenoPCA, input.pca, output.pdf34, params.width, params.height, log)
        print(cmd2)
        os.system(cmd2)


################################################################
