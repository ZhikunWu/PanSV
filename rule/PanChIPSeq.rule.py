import os 

SAMPLES = ["input", "VuCENH3"]
THREADS = 10
IN_PATH = "/home/wuzhikun/Project/Faba"

rule all:
    input:
        expand(IN_PATH + '/ChIPSeq/mapping/{sample}_deduplicated.bam', sample=SAMPLES),
        IN_PATH + '/ChIPSeq/macs'


rule mapping:
    input:
        fq = IN_PATH + "/raw/{sample}_combine.fastq.gz",
    output:
        sam = temp(IN_PATH + "/ChIPSeq/mapping/{sample}.map.sam"),
    params:
        RefGenome = "/home/wuzhikun/Project/PanSV/BACKUP/Assembly/Fengchan6/Vigna_unguiculata_assembly.fasta",
    threads:
        THREADS
    log:
        IN_PATH + "/log/{sample}.mapping.log",
    run:
        shell("bwa mem -t {threads} {params.RefGenome} {input.fq} > {output.sam} 2>{log}")


rule SAM2BAM:
    input:
        sam = IN_PATH + "/ChIPSeq/mapping/{sample}.map.sam",
    output:
        tempbam = temp(IN_PATH + '/ChIPSeq/mapping/{sample}_temp.bam'),
        bam = temp(IN_PATH + '/ChIPSeq/mapping/{sample}_align.bam'),
    threads:
        THREADS
    log:
        IN_PATH + "/log/{sample}.SAM2BAM.log",
    run:
        shell('sambamba view --nthreads {threads} --sam-input --format bam  -o {output.tempbam} {input.sam} >{log} 2>&1' )
        shell('samtools sort --threads {threads} -o {output.bam} {output.tempbam} 2>>{log}')



rule MarkDuplicates:
    input:
        bam = IN_PATH + '/ChIPSeq/mapping/{sample}_align.bam',
    output:
        bam = IN_PATH + '/ChIPSeq/mapping/{sample}_deduplicated.bam',
        metrics = IN_PATH + '/ChIPSeq/mapping/{sample}_dup_metrics.txt',
    params:
        picard = "/home/wuzhikun/anaconda3/envs/WGS/bin/picard.jar",
    threads:
        THREADS
    log:
        IN_PATH + "/log/{sample}.markDuplicate.log",
    run:
        ### -Xmx30g
        shell('java -Xmx10g -jar {params.picard}  MarkDuplicates INPUT={input.bam} OUTPUT={output.bam} METRICS_FILE={output.metrics} REMOVE_DUPLICATES=true ASSUME_SORTED=true >{log} 2>&1')
        shell("samtools index {output.bam}")



rule calleak:
    input:
        control = IN_PATH + '/ChIPSeq/mapping/input_deduplicated.bam',
        treat = IN_PATH + '/ChIPSeq/mapping/VuCENH3_deduplicated.bam',
    output:
        bed = directory(IN_PATH + '/ChIPSeq/macs'),
        peak = IN_PATH + '/ChIPSeq/macs/VuCENH3_control_peaks.broadPeak',
    params:
        outDir = IN_PATH + '/ChIPSeq/macs',
    log:
        IN_PATH + "/log/callPeak.log",
    run:
        # shell("macs2 callpeak -t {input.treat} -c {input.control} -f BAM -g 520000000 --outdir {params.outDir} --name VuCENH3_control  --broad --broad-cutoff 0.1  &>{log}")
        cmd = "source activate macs && macs2 callpeak -t %s -c %s -f BAM -g 520000000 --outdir %s --name VuCENH3_control  --broad --broad-cutoff 0.1  >%s 2>&1" % (input.treat, input.control, params.outDir, log)
        print(cmd)
        os.system(cmd)

### https://github.com/KoesGroup/Snakemake_ChIPseq_PE/blob/master/rules/deeptools_post_processing.smk
# rule bamCov:
#     input:
#         bam = IN_PATH + '/ChIPSeq/mapping/{sample}_deduplicated.bam',
#     output:
#         bw = IN_PATH + '/ChIPSeq/mapping/{sample}_deduplicated.bw',
#     run:
#         shell("bamCoverage  --bam {input.bam} -o {output.bw} ")
