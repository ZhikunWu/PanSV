rule RenameGenome:
    input:
        assembly = IN_PATH + "/Assembly/RagTag/{sample}/ragtag.scaffold.fasta",
    output:
        assembly = IN_PATH + "/Cactus/Assembly/{sample}.assembly.fasta",
    params:
        RenameScaffold = SCRIPT_DIR + "/RenameScaffold.py",
    threads:
        THREADS
    run:
        shell("python {params.RenameScaffold} --assembly {input.assembly} --out {output.assembly}")


rule fai:
    input:
        assembly = IN_PATH + "/Cactus/Assembly/{sample}.assembly.fasta",
    output:
        fai = IN_PATH + "/Cactus/Assembly/{sample}.assembly.fasta.fai",
    run:
        shell("samtools faidx {input.assembly}")

rule faiMerge:
    input:
        fai = expand(IN_PATH + "/Cactus/Assembly/{sample}.assembly.fasta.fai", sample=SAMPLES),
    output:
        fai = IN_PATH + "/Cactus/Assembly/All.assembly.length.txt",
    run:
        FAI = " ".join(sorted(input.fai))
        shell("cat {FAI} > {output.fai}")


rule minigraph:
    input:
        genome = IN_PATH + "/Cactus/Assembly/Fengchan6.assembly.fasta",
        assembly = expand(IN_PATH + "/Cactus/Assembly/{sample}.assembly.fasta", sample=SAMPLES),
    output:
        gfa = IN_PATH + "/Cactus/minigraph/Assembly.minigraph.gfa",
    threads:
        THREADS
    run:
        Assembly = " ".join(sorted(input.assembly))
        shell("minigraph -cxggs -t {threads} {input.genome} {Assembly} > {output.gfa}")

rule gaf2paf:
    input:
        length = IN_PATH + "/Cactus/Assembly/All.assembly.length.txt",
        gfa = IN_PATH + "/Cactus/minigraph/Assembly.minigraph.gfa",
        assembly = IN_PATH + "/Cactus/Assembly/{sample}.assembly.fasta",
    output:
        gaf = IN_PATH + "/Cactus/minigraph/{sample}.assembly.gaf",
        paf = IN_PATH + "/Cactus/minigraph/{sample}.assembly.paf",
    run:
        shell("minigraph -cxasm {input.gfa} {input.assembly}  > {output.gaf}")
        shell("gaf2paf {output.gaf} -l {input.length} > {output.paf}")


rule gfa2fa:
    input:
        gfa = IN_PATH + "/Cactus/minigraph/Assembly.minigraph.gfa",
    output:
        fa = IN_PATH + "/Cactus/minigraph/Assembly.minigraph.fa",
    run:
        shell("gfatools gfa2fa -s {input.gfa} > {output.fa}")
