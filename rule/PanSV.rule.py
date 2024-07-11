rule ddID:
    input:
        cuteSV = IN_PATH + "/SVCall/SV/{sample}.cutesv.vcf",
        svision = IN_PATH + "/SVCall/SV/{sample}.svision.s3.graph.vcf",
    output:
        cuteSV = IN_PATH + "/SVCall/addID/{sample}_cuteSV.vcf",
        svision = IN_PATH + "/SVCall/addID/{sample}_svision.vcf",
    params:
        addID = SCRIPT_DIR + "/addID.py",
    run:
        shell("python {params.addID} --vcf {input.cuteSV} --out {output.cuteSV} --name None ")
        shell("python {params.addID} --vcf {input.svision} --out {output.svision} --name svision ")



rule Mergemultiple2:
    input:
        cuteSV = IN_PATH + "/SVCall/addID/{sample}_cuteSV.vcf",   
        sniffles = IN_PATH + "/SVCall/SV/{sample}.sniffles.vcf",
        svision = IN_PATH + "/SVCall/addID/{sample}_svision.vcf",
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


