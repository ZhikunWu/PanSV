
rule ConvertPed2:
    input:
        vcf = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf005.filtmiss.recode.filtHeter.vcf",
    output:
        Ped = IN_PATH + "/NGS/Structure2/Sample_jointcall.twoallel.maf005.ped",
        Map = IN_PATH + "/NGS/Structure2/Sample_jointcall.twoallel.maf005.map",
    params:
        outPrefix = IN_PATH + "/NGS/Structure2/Sample_jointcall.twoallel.maf005",
        Map = IN_PATH + "/NGS/Structure2/Sample_jointcall.twoallel.maf005-1.map",
    log:
        IN_PATH + "/log/ConvertPed2.log",
    run:
        shell("vcftools --vcf {input.vcf} --plink --out {params.outPrefix} > {log} 2>&1")
        cmd = """ cut -f 2 %s | sed 's/:/\t/' | awk '{print $1"\t"$1":"$2"\t"0"\t"$2}' > %s """ % (output.Map, params.Map)
        os.system(cmd)
        shell("mv {params.Map} {output.Map}")




rule PruneIN2:
    input:
        Ped = IN_PATH + "/NGS/Structure2/Sample_jointcall.twoallel.maf005.ped",
        Map = IN_PATH + "/NGS/Structure2/Sample_jointcall.twoallel.maf005.map",
    output:
        prune = IN_PATH + "/NGS/Structure2/Sample_variant_PruneLD.prune.in",
        bed = IN_PATH + "/NGS/Structure2/Sample_variant_PruneLD.bed",
    params:
        outPrefix = IN_PATH + "/NGS/Structure2/Sample_jointcall.twoallel.maf005",
        prunePrefix = IN_PATH + "/NGS/Structure2/Sample_variant_PruneLD",
    log:
        IN_PATH + "/log/PruneIN2.log",
    run:
        shell("plink --file {params.outPrefix} --indep-pairwise 50 10 0.5 --out {params.prunePrefix} --allow-extra-chr  > {log} 2>&1")
        shell("plink --file {params.outPrefix} --extract {output.prune} -make-bed  --out {params.prunePrefix} --allow-extra-chr >> {log} 2>&1")



rule Structure2:
    input:
        bed = IN_PATH + "/NGS/Structure2/Sample_variant_PruneLD.bed",
    output:
        p = IN_PATH + "/NGS/Structure2/Structure/struct.{popu}.meanP",
    params:
        prunePrefix = IN_PATH + "/NGS/Structure2/Sample_variant_PruneLD",
        outPrefix = IN_PATH + "/NGS/Structure2/Structure/struct",
    log:
        IN_PATH + "/log/Structure2_{popu}.log",
    run:
        cmd = "source activate Assembly && python /home/wuzhikun/anaconda3/envs/Assembly/bin/structure.py -K %s --input %s --out %s --full --seed 100 " % (wildcards.popu, params.prunePrefix, params.outPrefix)
        os.system(cmd)
    

rule chooseK2:
    input:
        p = expand(IN_PATH + "/NGS/Structure2/Structure/struct.{popu}.meanP", popu=POPS),
    output:
        chooseK = IN_PATH + "/NGS/Structure2/Structure/chooseK.txt",
    params:
        outPrefix = IN_PATH + "/NGS/Structure2/Structure/struct",
    run:
        cmd = "source activate Assembly && python /home/wuzhikun/anaconda3/envs/Assembly/bin/chooseK.py --input %s > %s" % (params.outPrefix, output.chooseK)
        os.system(cmd)



######################### LD  ###########################
rule groupVCF2:
    input:
        vcf = IN_PATH + "/NGS/Variant/Sample_jointcall.twoallel.maf001.GT.filtHeter.vcf",
    output:
        vcf1 = IN_PATH + "/NGS/Structure2/MAF001/Q1_samples.vcf",
        vcf2 = IN_PATH + "/NGS/Structure2/MAF001/Q2_samples.vcf",
        vcf3 = IN_PATH + "/NGS/Structure2/MAF001/Q3_samples.vcf",
        vcf4 = IN_PATH + "/NGS/Structure2/MAF001/Q4_samples.vcf",
    threads:
        THREADS
    params:
        q1 = IN_PATH + "/NGS/Structure2/Q1_sample.txt",
        q2 = IN_PATH + "/NGS/Structure2/Q2_sample.txt",
        q3 = IN_PATH + "/NGS/Structure2/Q3_sample.txt",
        q4 = IN_PATH + "/NGS/Structure2/Q4_sample.txt",
    run:
        shell("bcftools view --threads {threads} --samples-file {params.q1} -o {output.vcf1} {input.vcf} ")
        shell("bcftools view --threads {threads} --samples-file {params.q2} -o {output.vcf2} {input.vcf} ")
        shell("bcftools view --threads {threads} --samples-file {params.q3} -o {output.vcf3} {input.vcf} ")
        shell("bcftools view --threads {threads} --samples-file {params.q4} -o {output.vcf4} {input.vcf} ")



rule LD2:
    input:
        vcf1 = IN_PATH + "/NGS/Structure2/MAF001/Q1_samples.vcf",
        vcf2 = IN_PATH + "/NGS/Structure2/MAF001/Q2_samples.vcf",
        vcf3 = IN_PATH + "/NGS/Structure2/MAF001/Q3_samples.vcf",
        vcf4 = IN_PATH + "/NGS/Structure2/MAF001/Q4_samples.vcf",
    output:
        ld1 = IN_PATH + "/NGS/Structure2/MAF001/Q1/LD_{Chr}.ld.gz", 
        ld2 = IN_PATH + "/NGS/Structure2/MAF001/Q2/LD_{Chr}.ld.gz", 
        ld3 = IN_PATH + "/NGS/Structure2/MAF001/Q3/LD_{Chr}.ld.gz", 
        ld4 = IN_PATH + "/NGS/Structure2/MAF001/Q4/LD_{Chr}.ld.gz", 
    params:
        q1 = IN_PATH + "/NGS/Structure2/MAF001/Q1/LD_{Chr}",
        q2 = IN_PATH + "/NGS/Structure2/MAF001/Q2/LD_{Chr}",
        q3 = IN_PATH + "/NGS/Structure2/MAF001/Q3/LD_{Chr}",
        q4 = IN_PATH + "/NGS/Structure2/MAF001/Q4/LD_{Chr}",
    run:
        shell("plink --vcf {input.vcf1} --double-id --allow-extra-chr --set-missing-var-ids @:# --maf 0.01 --geno 0.1 --mind 0.5 --chr {wildcards.Chr} --thin 0.1 -r2 gz --ld-window 100 --ld-window-kb 1000 --ld-window-r2 0 --make-bed --out {params.q1}")
        shell("plink --vcf {input.vcf2} --double-id --allow-extra-chr --set-missing-var-ids @:# --maf 0.01 --geno 0.1 --mind 0.5 --chr {wildcards.Chr} --thin 0.1 -r2 gz --ld-window 100 --ld-window-kb 1000 --ld-window-r2 0 --make-bed --out {params.q2}")
        shell("plink --vcf {input.vcf3} --double-id --allow-extra-chr --set-missing-var-ids @:# --maf 0.01 --geno 0.1 --mind 0.5 --chr {wildcards.Chr} --thin 0.1 -r2 gz --ld-window 100 --ld-window-kb 1000 --ld-window-r2 0 --make-bed --out {params.q3}")
        shell("plink --vcf {input.vcf4} --double-id --allow-extra-chr --set-missing-var-ids @:# --maf 0.01 --geno 0.1 --mind 0.5 --chr {wildcards.Chr} --thin 0.1 -r2 gz --ld-window 100 --ld-window-kb 1000 --ld-window-r2 0 --make-bed --out {params.q4}")





rule MergeChr2:
    input:
        ld1 = expand(IN_PATH + "/NGS/Structure2/MAF001/Q1/LD_{Chr}.ld.gz", Chr = CHRS),
        ld2 = expand(IN_PATH + "/NGS/Structure2/MAF001/Q2/LD_{Chr}.ld.gz",  Chr = CHRS),
        ld3 = expand(IN_PATH + "/NGS/Structure2/MAF001/Q3/LD_{Chr}.ld.gz",  Chr = CHRS),
        ld4 = expand(IN_PATH + "/NGS/Structure2/MAF001/Q4/LD_{Chr}.ld.gz",  Chr = CHRS),
    output:
        ld1 = IN_PATH + "/NGS/Structure2/MAF001/Q1/All_Chr_merge.ld.gz",    
        ld2 = IN_PATH + "/NGS/Structure2/MAF001/Q2/All_Chr_merge.ld.gz",    
        ld3 = IN_PATH + "/NGS/Structure2/MAF001/Q3/All_Chr_merge.ld.gz",    
        ld4 = IN_PATH + "/NGS/Structure2/MAF001/Q4/All_Chr_merge.ld.gz",    
    params:
        ld1 = IN_PATH + "/NGS/Structure2/MAF001/Q1/All_Chr_merge.ld",
        ld2 = IN_PATH + "/NGS/Structure2/MAF001/Q2/All_Chr_merge.ld",
        ld3 = IN_PATH + "/NGS/Structure2/MAF001/Q3/All_Chr_merge.ld",
        ld4 = IN_PATH + "/NGS/Structure2/MAF001/Q4/All_Chr_merge.ld",
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
rule LDStats2:
    input:
        ld1 = IN_PATH + "/NGS/Structure2/MAF001/Q1/All_Chr_merge.ld.gz",   
        ld2 = IN_PATH + "/NGS/Structure2/MAF001/Q2/All_Chr_merge.ld.gz",    
        ld3 = IN_PATH + "/NGS/Structure2/MAF001/Q3/All_Chr_merge.ld.gz",    
        ld4 = IN_PATH + "/NGS/Structure2/MAF001/Q4/All_Chr_merge.ld.gz",  
    output:
        q1 = IN_PATH + "/NGS/Structure2/MAF001/Q1/All_Chr_merge.ld_decay_bins",  
        q2 = IN_PATH + "/NGS/Structure2/MAF001/Q2/All_Chr_merge.ld_decay_bins",   
        q3 = IN_PATH + "/NGS/Structure2/MAF001/Q3/All_Chr_merge.ld_decay_bins",     
        q4 = IN_PATH + "/NGS/Structure2/MAF001/Q4/All_Chr_merge.ld_decay_bins",   
    params:
        ld_decay_calc = "/home/wuzhikun/github/speciationgenomics/ld_decay_calc.py",
        q1 = IN_PATH + "/NGS/Structure2/MAF001/Q1/All_Chr_merge",
        q2 = IN_PATH + "/NGS/Structure2/MAF001/Q2/All_Chr_merge",
        q3 = IN_PATH + "/NGS/Structure2/MAF001/Q3/All_Chr_merge",
        q4 = IN_PATH + "/NGS/Structure2/MAF001/Q4/All_Chr_merge",
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



rule ChrAve2:
    input:
        q1 = IN_PATH + "/NGS/Structure2/MAF001/Q1/All_Chr_merge.ld_decay_bins",  
        q2 = IN_PATH + "/NGS/Structure2/MAF001/Q2/All_Chr_merge.ld_decay_bins",   
        q3 = IN_PATH + "/NGS/Structure2/MAF001/Q3/All_Chr_merge.ld_decay_bins",  
        q4 = IN_PATH + "/NGS/Structure2/MAF001/Q4/All_Chr_merge.ld_decay_bins",  
        #a = IN_PATH + "/NGS/Structure2/MAF001/all/All_Chr_merge.ld_decay_bins",
    output:
        q1 = IN_PATH + "/NGS/Structure2/MAF001/All_group_all_Chr_merge.ld_decay_bins_all.txt",  
    run:
        #files = ",".join([input.a, input.q1, input.q2, input.q3, input.q4])
        files = ",".join([input.q1, input.q2, input.q3, input.q4])
        Chr_LD_average(files, output.q1)


rule LDPlot2:
    input:
        q1 = IN_PATH + "/NGS/Structure2/MAF001/All_group_all_Chr_merge.ld_decay_bins_all.txt",  
    output:
        q1 = IN_PATH + "/NGS/Structure2/MAF001/All_group_all_Chr_merge.ld_decay_bins_all.pdf",
    params:
        LDdecay = SCRIPT_DIR + "/LDdecay.R",
    run:
        shell("Rscript {params.LDdecay} --input {input.q1}  --pdf {output.q1} --width 5 --height 4")


################################################################

