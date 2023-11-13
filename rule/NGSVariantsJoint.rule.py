#################### joint call ###########
rule tabix:
    input:
        vcf = IN_PATH + '/NGS/VCF/{sample}/{sample}_{Chr}.g.vcf.gz',
    output:
        tbi = IN_PATH + '/NGS/VCF/{sample}/{sample}_{Chr}.g.vcf.gz.tbi',
    run:
        shell("tabix -p vcf {input.vcf}")


rule GenomicsDBImport:
    input:
        gvcf = expand(IN_PATH + '/NGS/VCF/{sample}/{sample}_{Chr}.g.vcf.gz', sample=SAMPLES, Chr=CHRS),
    output:
        genomicsdb = directory(IN_PATH + '/NGS/Joint/genomicsdb/{Chr}_genomicsdb'),
    params:
        GATK4 = config['GATK4'],
    threads:
        THREADS
    log:
        IN_PATH + "/log/{Chr}.GenomicsDBImport.log",
    run:
        Files = input.gvcf
        Targets = []
        for f in Files:
            if wildcards.Chr in f:
                Targets.append(f)
        GVCFs = " -V ".join(sorted(Targets))
        print(sorted(Targets))
        shell("gatk GenomicsDBImport   -V {GVCFs} --genomicsdb-workspace-path {output.genomicsdb} -L {wildcards.Chr} > {log} 2>&1")
    

rule jointCall:
    input:
        genomicsdb = IN_PATH + '/NGS/Joint/genomicsdb/{Chr}_genomicsdb',
    output:
        vcf = IN_PATH + '/NGS/Joint/VCF/{Chr}_Fengchan6_joint_variant.vcf',
    params:
        GATK4 = config['GATK4'],
        REF = config['RefGenome'],
    threads:
        THREADS
    log:
        IN_PATH + "/log/{Chr}.jointCall.log",
    run:
        shell("gatk GenotypeGVCFs  -R {params.REF} -V gendb://{input.genomicsdb} -O {output.vcf} --allow-old-rms-mapping-quality-annotation-data > {log} 2>&1")

##########################################################
