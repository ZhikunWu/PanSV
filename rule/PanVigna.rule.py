rule filtVariant:
    input:
        vcf = IN_PATH + "/NGS/Joint/VCF/Chrom_Fengchan6_joint_variant.vcf.gz",
    output:
        vcf = IN_PATH + "/NGS/Population/Chr_joint_variant.vcf",
    log:
        IN_PATH + "/log/filtVariant.log"
    run:
        shell("vcftools --gzvcf {input.vcf} --max-missing 0.1 --min-alleles 2  --max-alleles 2 --min-meanDP 10 --out {output.vcf} > {log}")

rule vcfStat:
    input:
        vcf = IN_PATH + "/NGS/Joint/VCF/Chrom_Fengchan6_joint_variant.vcf.gz",
    output:
        stat = IN_PATH + "/NGS/Joint/VCF/Chrom_Fengchan6_joint_variant.stat.txt",
    run:
        shell("rtg vcfstats {input.vcf} > {output.stat}")




rule vcfGT:
    input:
        vcf = IN_PATH + "/NGS/Population/Chr_joint_variant.vcf",
    output:
        vcf = IN_PATH + "/NGS/Population/Chr_joint_variant_GT.vcf",
    run:
        shell("bcftools annotate -x ^FORMAT/GT {input.vcf} > {output.vcf}")


