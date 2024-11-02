'''
################### merge for each individual ##################
rule ddID:
    input:
        cuteSV = IN_PATH + "/SVCall/SV/{sample}.cutesv.vcf.gz",
        sniffles = IN_PATH + "/SVCall/SV/{sample}.sniffles.vcf.gz",
        svision = IN_PATH + "/SVCall/SV/{sample}.svision.s3.graph.vcf.gz",
    output:
        cuteSV = IN_PATH + "/SVCall/Merge/{sample}_cutesv.vcf",
        sniffles = temp(IN_PATH + "/SVCall/Merge/{sample}_sniffles.vcf"),
        svision = temp(IN_PATH + "/SVCall/Merge/{sample}_svision.vcf"),
        sniffles2 = IN_PATH + "/SVCall/Merge/{sample}_sniffles.addID.vcf",
        svision2 = IN_PATH + "/SVCall/Merge/{sample}_svision.addID.vcf",
    params:
        addID = SRC_DIR + "/addID.py",
    threads:
        THREADS
    run:
        shell("pigz -dc -p {threads} {input.cuteSV} > {output.cuteSV}")
        shell("pigz -dc -p {threads} {input.sniffles} > {output.sniffles}")
        shell("pigz -dc -p {threads} {input.svision} > {output.svision}")
        shell("python {params.addID} --vcf {output.sniffles} --out {output.sniffles2} --name sniffles ")
        shell("python {params.addID} --vcf {output.svision} --out {output.svision2} --name svision ")




rule Mergemultiple2:
    input:
        cuteSV = IN_PATH + "/SVCall/Merge/{sample}_cutesv.vcf",    
        sniffles = IN_PATH + "/SVCall/Merge/{sample}_sniffles.addID.vcf",
        svision = IN_PATH + "/SVCall/Merge/{sample}_svision.addID.vcf",
    output:
        file = IN_PATH + "/SVCall/Merge/{sample}_file_list.txt",
        vcf = IN_PATH + "/SVCall/Merge/{sample}_merge.vcf",
    threads:
        THREADS
    log:
        IN_PATH + "/log/mergeSV2_{sample}.log"
    run:
        shell("echo {input.cuteSV} > {output.file}")
        shell("echo {input.sniffles} >> {output.file}")
        shell("echo {input.svision} >> {output.file}")
        shell("jasmine --output_genotypes threads={threads} file_list={output.file} out_file={output.vcf} > {log} 2>&1")



rule MergeFilter2:
    input:
        vcf = IN_PATH + "/SVCall/Merge/{sample}_merge.vcf",
    output:
        vcf = IN_PATH + "/SVCall/Merge/{sample}_merge_filter.vcf",
    params:
        MergeSVTarget = SCRIPT_DIR + "/MergeSVTarget.py"
    run:
        shell("python {params.MergeSVTarget} --vcf {input.vcf} --out {output.vcf}")


rule Mergemultiple3:
    input:
        vcf = expand(IN_PATH + "/SVCall/Merge/{sample}_merge_filter.vcf", sample=SAMPLES),
    output:
        file = IN_PATH + "/SVCall/Merge/All_samples_file_list.txt",
        vcf = IN_PATH + "/SVCall/Merge/All_samples_SV_merge.vcf",
    threads:
        THREADS
    log:
        IN_PATH + "/log/mergeSV3.log"
    run:
        Files = sorted(input.vcf)
        for i in range(len(Files)):
            f = Files[i]
            if i == 0:
                cmd = "echo %s > %s" % (f, output.file)
            else:
                cmd = "echo %s >> %s" % (f, output.file)
            os.system(cmd)
        shell("jasmine --output_genotypes threads={threads} file_list={output.file} out_file={output.vcf} > {log} 2>&1")


###################################################################
'''


'''
############ graph pan SV ########

rule sortedVCF:
    input:
        vcf = IN_PATH + "/SVCall/Merge/All_samples_SV_merge.vcf",
    output:
        vcf = IN_PATH + "/SVCall/Merge/All_samples_SV_merge_sorted.vcf",
    run:
        shell("grep '^#' {input.vcf} > {output.vcf} ")
        shell("grep -v '^#' {input.vcf} | sort -k 1,1 -k 2,2n >> {output.vcf} ")




rule VCF2vgFormatTest:
    input:
        vcf = IN_PATH + "/SVCall/Merge/All_samples_SV_merge_sorted.vcf",
    output:
        vcf = IN_PATH + "/SVCall/PanGraph3/All_samples_SV_merge_vgformat.vcf",
    params:
        RefGenome = "/home/wuzhikun/Project/Vigna/Final/Final/Vigna_unguiculata_assembly.fasta",
        #vcf = IN_PATH + "/SVCall/PanGraph3/All_samples_SV_merge_vgformat.vcf",
        VCF2vgFormat = SCRIPT_DIR + "/VCF2vgFormat.py",
    threads:
        THREADS
    log:
        IN_PATH + "/log/VCF2vgFormat.log"
    run:
        shell("python {params.VCF2vgFormat} --vcf {input.vcf} --out {output.vcf} --fasta {params.RefGenome} > {log} 2>&1")
        #shell("bgzip {params.vcf}")
        #shell("tabix -p vcf {output.vcf}")

rule vcfSimple:
    input:
        vcf = IN_PATH + "/SVCall/PanGraph3/All_samples_SV_merge_vgformat.vcf",
    output:
        vcf = IN_PATH + "/SVCall/PanGraph3/All_samples_SV_merge_vgformat_simple.vcf",
    params:
        VcfSimple = SCRIPT_DIR + "/VcfSimple.py",
    run:
        shell("python {params.VcfSimple} --vcf {input.vcf} --out {output.vcf}")



rule simpleVGtest:
    input:
        vcf = IN_PATH + "/SVCall/PanGraph3/All_samples_SV_merge_vgformat_simple.vcf",
    output:
        vg = IN_PATH + "/SVCall/PanGraph3/All_samples_SV_merge_construct.vg",
    threads:
        THREADS
    params:
        RefGenome = "/home/wuzhikun/Project/Vigna/Final/Final/Vigna_unguiculata_assembly.fasta",
    log:
        IN_PATH + "/log/simpleVG.log"
    run:
        shell("vg construct -p  --threads {threads}  -r  {params.RefGenome} -v {input.vcf} > {output.vg} 2>> {log}")


##################################


################ call NGS SV ##############


rule GiraffeIndex:
    input:
        vcf = IN_PATH + "/SVCall/PanGraph3/All_samples_SV_merge_vgformat_simple.vcf",
    output:
        dist = IN_PATH + "/SVCall/PanGraph3/All_samples_SV_merge_giraffe.dist",
        gbz = IN_PATH + "/SVCall/PanGraph3/All_samples_SV_merge_giraffe.giraffe.gbz",
        m = IN_PATH + "/SVCall/PanGraph3/All_samples_SV_merge_giraffe.min",
    params:
        RefGenome = "/home/wuzhikun/Project/Vigna/Final/Final/Vigna_unguiculata_assembly.fasta",
        outPrefix = IN_PATH + "/SVCall/PanGraph3/All_samples_SV_merge_giraffe",
    log:
        IN_PATH + "/log/GiraffeIndex.log"
    run:
        shell("/home/wuzhikun/software/vg autoindex --workflow giraffe -r {params.RefGenome} -v {input.vcf} -p {params.outPrefix} > {log} 2>&1")


rule vgConvert:
    input:
        gbz = IN_PATH + "/SVCall/PanGraph3/All_samples_SV_merge_giraffe.giraffe.gbz",
    output:
        xg = IN_PATH + "/SVCall/PanGraph3/All_samples_SV_merge_giraffe.giraffe.xg",
    run:
        shell("/home/wuzhikun/software/vg convert -x {input.gbz}  > {output.xg}")



rule snarls:
    input:
        gbz = IN_PATH + "/SVCall/PanGraph3/All_samples_SV_merge_giraffe.giraffe.gbz",
    output:
        snarls = IN_PATH + "/SVCall/PanGraph3/All_samples_SV_merge_giraffe.giraffe.snarls.pb",
    threads:
        THREADS
    run:
        shell("/home/wuzhikun/software/vg snarls -t {threads} {input.gbz} > {output.snarls}")
'''





rule Giraffe:
    input:
        dist = IN_PATH + "/SVCall/PanGraph3/All_samples_SV_merge_giraffe.dist",
        gbz = IN_PATH + "/SVCall/PanGraph3/All_samples_SV_merge_giraffe.giraffe.gbz",
        m = IN_PATH + "/SVCall/PanGraph3/All_samples_SV_merge_giraffe.min",
        xg = IN_PATH + "/SVCall/PanGraph3/All_samples_SV_merge_giraffe.giraffe.xg",
        R1 = IN_PATH + "/clean/{sample}_NGS.R1.fastq.gz",
        R2 = IN_PATH + "/clean/{sample}_NGS.R2.fastq.gz",
    output:
        gam = IN_PATH + "/SVCall/Giraffe3/{sample}.giraffe.gam",
    threads:
        THREADS
    log:
        IN_PATH + "/log/Giraffe_{sample}.log"
    run:
        shell("/home/wuzhikun/software/vg giraffe --gbz-name {input.gbz} --dist-name {input.dist} --minimizer-name {input.m} --xg-name {input.xg} --fastq-in {input.R1} --fastq-in {input.R2}   --threads {threads} --sample {wildcards.sample}  -o gam  >  {output.gam} 2>{log}")
        








rule pack:
    input:
        gbz = IN_PATH + "/SVCall/PanGraph3/All_samples_SV_merge_giraffe.giraffe.gbz",
        gam = IN_PATH + "/SVCall/Giraffe3/{sample}.giraffe.gam",
    output:
        pack = IN_PATH + "/SVCall/Giraffe3/{sample}.giraffe.pack",
    threads:
        THREADS
    log:
        IN_PATH + "/log/pack_{sample}.log"
    run:
        shell("/home/wuzhikun/software/vg pack -t {threads} -x {input.gbz} -g {input.gam} -Q5 -o {output.pack} > {log} 2>&1")




rule gamstat:
    input:
        gam = IN_PATH + "/SVCall/Giraffe3/{sample}.giraffe.gam",
    output:
        stats = IN_PATH + "/SVCall/Giraffe3/{sample}.giraffe.stats.txt",
    run:
        shell("/home/wuzhikun/software/vg stats -a {input.gam} > {output.stats}")



rule call:
    input:
        gbz = IN_PATH + "/SVCall/PanGraph3/All_samples_SV_merge_giraffe.giraffe.gbz",
        snarls = IN_PATH + "/SVCall/PanGraph3/All_samples_SV_merge_giraffe.giraffe.snarls.pb",
        pack = IN_PATH + "/SVCall/Giraffe3/{sample}.giraffe.pack",
    output:
        vcf = IN_PATH + "/SVCall/Giraffe3/{sample}.giraffe.vcf.gz",
    params:
        vcf = IN_PATH + "/SVCall/Giraffe3/{sample}.giraffe.vcf",
    threads:
        THREADS
    run:
        shell("/home/wuzhikun/software/vg call -t {threads} -k {input.pack} --min-support 2,2 -z -a -r {input.snarls} -z -s {wildcards.sample} {input.gbz} > {params.vcf}")
        shell("bgzip {params.vcf}")
        shell("tabix -p vcf {output.vcf}")




rule mergeSV:
    input:
        vcf = expand(IN_PATH + "/SVCall/Giraffe3/{sample}.giraffe.vcf.gz", sample=SAMPLES),
    output:
        vcf = IN_PATH + "/SVCall/Giraffe3/Sample884.giraffe.merge.vcf",
    run:
        files = " ".join(sorted(input.vcf))
        shell("bcftools merge --threads 24  {files} > {output.vcf} ")


###########################################
