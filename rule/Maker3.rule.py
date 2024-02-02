### [使用MAKER进行注释: 学习MAKER的配置参数](https://blog.csdn.net/u012110870/article/details/107080165)

def modify_opts(in_file, out_file, genome_file, est_file, pro_file):
    in_h = open(in_file, "r")
    out_h = open(out_file, "w")
    for line in in_h:
        line = line.strip()
        if line.startswith("genome="):
            newLine = "genome=%s" % genome_file
        elif line.startswith("est="):
            newLine = "est=%s" % est_file
        elif line.startswith("protein="):
            newLine = "protein=%s" % pro_file
        else:
            newLine = line
        out_h.write("%s\n" % newLine)
    in_h.close()
    out_h.close()

rule createConfig:
    input:
        genome = 
        est = 
        protein = 
        est_gff =  #cufflinks/stringtie组装结果
    output:
        opts = IN_PATH + "/GenePrediction/Maker3/{sample}/maker_opts.ctl",
        out = IN_PATH + "/GenePrediction/Maker3/{sample}/{sample}.maker.output/{sample}_master_datastore_index.log",
    threads:
        THREADS
    params:
        bopts = IN_PATH + "/Maker3/config/config/maker_bopts.ctl",
        evm = IN_PATH + "/Maker3/config/config/maker_evm.ctl",
        exe = IN_PATH + "/Maker3/config/config/maker_exe.ctl",
        opts = IN_PATH + "/Maker3/config/config/maker_opts.ctl",
        outDir = IN_PATH + "/GenePrediction/Maker3/{sample}",
    log:
        IN_PATH + "/log/createConfig_{sample}.log",
    run:
        if not os.path.exists(params.outDir):
            os.makedirs(params.outDir)
        modify_opts(params.opts, output.opts, input.genome, input.est, input.protein)
        shell("cp {params.evm} {params.exe} {params.bopts} {params.outDir}")
        # shell("cd {params.outDir} && maker  -CTL")
        shell("cd {params.outDir} && mpiexec -n {threads} maker -base {wildcards.sample}   > {log} 2>&1")



rule mergefile:
    input:
        index = IN_PATH + "/GenePrediction/Maker3/{sample}/{sample}.maker.output/{sample}_master_datastore_index.log",
    output:
        gff = IN_PATH + "/GenePrediction/Maker3/{sample}/{sample}.all.gff",
        protein = IN_PATH + "/GenePrediction/Maker3/{sample}/{sample}.all.maker.protein.fasta",
        transcript = IN_PATH + "/GenePrediction/Maker3/{sample}/{sample}.all.maker.transcripts.fasta",
    params:
        outDir = IN_PATH + "/GenePrediction/Maker3/{sample}",
    run:
        shell("cd {params.outDir} && fasta_merge -d {input.index}")
        shell("cd {params.outDir} && gff3_merge -d {input.index}")




