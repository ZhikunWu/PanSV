
rule miniprot:
    input:
        assembly = IN_PATH + "/Scaffold/{sample}.genome.fasta",
        protein = "/home/wuzhikun/data/Vigna/Protein/Faba_10genome_protein_cdhit_0.9",
    output:
        gff = IN_PATH + "/GenePrediction/homo/miniprot/{sample}_align.gff",
    threads:
        THREADS
    log:
        IN_PATH + "/log/miniprot_{sample}.log",
    run:
        shell("miniprot --gff  -t {threads} {input.assembly} {input.protein} > {output.gff}")


rule gffFilt:
    input:
        gff = IN_PATH + "/GenePrediction/homo/miniprot/{sample}_align.gff",
    output:
        gff = IN_PATH + "/GenePrediction/homo/miniprot/{sample}_align_filt.gff",
    params:
        miniProtIdenFilter = SRC_DIR + "/miniProtIdenFilter.py"
    run:
        shell("python {params.miniProtIdenFilter} --gff {input.gff} --out {output.gff} --idenThreshold 0.9") 



