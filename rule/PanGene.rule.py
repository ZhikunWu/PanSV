rule diamondIndex:
    input:
        Ref = IN_PATH + "/GenePrediction/Gene/Fengchan6/Fengchan6.pep.fasta",
    output:
        dmnd = IN_PATH + "/GenePrediction/Gene/Fengchan6/Fengchan6.pep.dmnd",
    params:
        outPrefix = IN_PATH + "/GenePrediction/Gene/Fengchan6/Fengchan6.pep",
    run:
        shell("diamond makedb --in {input.Ref}  -d {params.outPrefix}")


rule diamond:
    input:
        dmnd = IN_PATH + "/GenePrediction/Gene/Fengchan6/Fengchan6.pep.dmnd",
        protein = IN_PATH + "/GenePrediction/Gene/{sample}/{sample}.pep.fasta",
    output:
        blast = IN_PATH + "/Protein/Diamond/{sample}/{sample}_blast.txt",
    threads:
        THREADS
    run:
        shell("diamond blastp --query {input.protein} --db {input.dmnd} --out {output.blast} --outfmt 6 --threads {threads} ")


rule bestHomo:
    input:
        blast = IN_PATH + "/Protein/Diamond/{sample}/{sample}_blast.txt",
        protein = IN_PATH + "/GenePrediction/Gene/{sample}/{sample}.pep.fasta",
    output:
        homo = IN_PATH + "/Protein/Diamond/{sample}/{sample}_homologous.txt",
        protein = IN_PATH + "/Protein/Diamond/{sample}/{sample}_special_protein.fasta",
    params:
        diamondHomoGene = SCRIPT_DIR + "/diamondHomoGene.py",
    run:
        shell("python {params.diamondHomoGene} --blast {input.blast} --protein {input.protein}  --out {output.protein} --homologous {output.homo}  --idenThreshold 80 ")



rule mergeSpecial:
    input:
        protein = expand(IN_PATH + "/Protein/Diamond/{sample}/{sample}_special_protein.fasta", sample=SAMPLES),
    output:
        protein = IN_PATH + "/Protein/Diamond/Merge_special_protein.fasta",
    run:
        Proteins = " ".join(sorted(input.protein))
        shell("cat {Proteins} > {output.protein}")



rule mergecluster:
    input:
        protein = IN_PATH + "/Protein/Diamond/Merge_special_protein.fasta",
    output:
        cluster = IN_PATH + "/Protein/Diamond/Merge_special_protein.culster",
    log:
        IN_PATH + "/log/cdhit.log",
    run:
        shell("cd-hit -i {input.protein} -o {output.cluster}  -c 0.8 -s 0.8   -n 5 -M 100000 -d 0 -T 20 > {log} 2>&1")

