

################################### Methylation #####################
rule mergeFast5:
    input:
        fast5 = IN_PATH + "/fast5_pass/{sample}.fast5_pass",
    output:
        fast5 = IN_PATH + "/Methylation/fast5/{sample}/filename_mapping.txt",
    params:
        multiple2single = "/home/wuzhikun/anaconda3/envs/Methy/bin/multi_to_single_fast5",
        outDir = IN_PATH + "/Methylation/fast5/{sample}",
    log:
        IN_PATH + '/log/mergeFast5_{sample}.log' 
    run:
        ### /home/jiangzehang/biosoft/ont-guppy/bin
        ### Methy
        ### /home/wuzhikun/anaconda3/envs/Methy/bin/multi_to_single_fast5 -i /home/wuzhikun/data/Vigna/ONTLong/fast5_test -s /home/wuzhikun/data/Vigna/ONTLong/fast5_test2
        shell("{params.multiple2single} -i {input.fast5} -s {params.outDir} --threads {threads} > {log} 2>&1")


rule MergeFastq:
    input:
        fastq = IN_PATH + "/raw/{sample}_ONT.fastq.gz",
    output:
        fastq = IN_PATH + "/fastq_pass/{sample}/pass_merge.fastq",
    log:
        IN_PATH + '/log/MergeFastq_{sample}.log' 
    threads:
        THREADS
    run:
        shell("pigz -p {threads} -dc {input.fastq} > {output.fastq} 2>{log}")


rule nanopolishIndex:
    input:
        fast5 = IN_PATH + "/fast5_pass/{sample}.fast5_pass",
        fastq = IN_PATH + "/fastq_pass/{sample}/pass_merge.fastq",
    output:
        fai = IN_PATH + "/fastq_pass/{sample}/pass_merge.fastq.index.fai",
    log:
        IN_PATH + '/log/nanopolishIndex_{sample}.log' 
    run:
        ### source activate /home/jiangzehang/miniconda2/envs/methylation
        ### output: .index, .index.fai, .index.gzi
        # shell("nanopolish index -d  {input.fast5}  {input.fastq}  > {log} 2>&1")
        cmd = "source activate /home/jiangzehang/miniconda2/envs/methylation && nanopolish index -d  %s  %s  > %s 2>&1" % (input.fast5, input.fastq, log)
        print(cmd)
        os.system(cmd)


rule alignFastq:
    input:
        fastq = IN_PATH + "/fastq_pass/{sample}/pass_merge.fastq",
    output:
        sam = temp(IN_PATH + "/Methylation/mapping/{sample}_fastq.sam"),
        bam = IN_PATH + "/Methylation/mapping/{sample}_fastq.bam",
    threads:
        THREADS
    params:
        Ref = IN_PATH + "/Final/Final/Vigna_unguiculata_assembly.fasta",
    log:
        IN_PATH + '/log/alignFastq_{sample}.log' 
    run:
        shell("minimap2 --MD -a -x map-ont -t {threads} {params.Ref} {input.fastq} > {output.sam} 2>{log}")
        shell("samtools view -Sb --threads {threads} {output.sam} | samtools sort --threads {threads} -o {output.bam}  -  2>>{log}")
        shell("samtools index {output.bam}")


rule callMethy:
    input:
        fastq = IN_PATH + "/fastq_pass/{sample}/pass_merge.fastq",
        fai = IN_PATH + "/fastq_pass/{sample}/pass_merge.fastq.index.fai",
        bam = IN_PATH + "/Methylation/mapping/{sample}_fastq.bam",
    output:
        calls = IN_PATH + "/Methylation/methy/{sample}/methylation_calls_{Chr}.tsv",
    params:
        # region = "Vu1:1,000,000-5,000,000",
        # region = "Vu01",
        outDir = IN_PATH + "/Methylation/methy",
        Ref = IN_PATH + "/Final/Final/Vigna_unguiculata_assembly.fasta",
    threads:
        THREADS
    log:
        IN_PATH + '/log/callMethy_{sample}_{Chr}.log'
    run:
        ### source activate /home/jiangzehang/miniconda2/envs/methylation
        # shell("nanopolish call-methylation -t {threads} -r {input.fastq} -b {input.bam} -g {input.Ref} -w {params.region} > {output.calls} 2>{log}")
        if not os.path.exists(params.outDir):
            os.makedirs(params.outDir)
        cmd = "source activate /home/jiangzehang/miniconda2/envs/methylation && nanopolish call-methylation -t %s -r %s -b %s -g %s -w %s > %s 2>%s" % (threads, input.fastq, input.bam, params.Ref, wildcards.Chr, output.calls, log)
        print(cmd)
        os.system(cmd)



# rule methFreq:
#     input:
#         calls = IN_PATH + "/Methylation/methy/methylation_calls_all_chrs.tsv",
#     output:
#         calls = IN_PATH + "/Methylation/methy/methylation_calls_all_chrs_meth_freq.tsv",
#         sig = IN_PATH + "/Methylation/methy/methylation_calls_all_chrs_meth_freq_sig.tsv",
#     params:
#         calculate_methylation_frequency = "/home/wuzhikun/software/nanopolish/scripts/calculate_methylation_frequency.py",
#     run:
#         shell("python {params.calculate_methylation_frequency} {input.calls} > {output.calls} 2>{log}")
#         cmd = "awk '{if($7 >=0.8){print $0}}' %s > %s" % (input.calls, input.sig)
#         print(cmd)
#         os.system(cmd)




# rule MethWinNumber:
#     input:
#         sig = IN_PATH + "/Methylation/methy/methylation_calls_all_chrs_meth_freq_sig.tsv",
#         genomeLen = IN_PATH + "/scaffold/final3/Vigna_unguiculata_assembly.fasta.split/Vigna_unguiculata_assembly.id.chr_length.txt",
#     output:
#         win = IN_PATH + "/Methylation/methy/methylation_calls_freq_sig_win_number.txt",
#     params:
#         WindowSVNumber = SRC_DIR + "/WindowSVNumber.py",
#         window = 50000,
#     log:
#         IN_PATH + '/log/WinNumber.log'
#     run:
#         shell("python {params.WindowSVNumber} --genome {input.genomeLen} --input {input.sig} --out {output.win} --window {params.window} --sliding 0 >{log} 2>&1")

# #####################################################################



# ################# region intersection #############

# rule overlapGene:
#     input:
#         meth = IN_PATH + "/Methylation/methy/methylation_calls_all_chrs_meth_freq_0.5.bed",
#         gene = IN_PATH + "/GenePrediction/gene/Vigna_unguiculata_assembly.all.gene.bed",
#     output:
#         gene = IN_PATH + "/Methylation/methy/gene/meth_gene_overlap.txt",
#     run:
#         shell("bedtools intersect -a {input.meth} -b {input.gene} -wa -wb > {output.gene}")




# # rule OverlapFilt:
# #     input:
# #         gene = IN_PATH + "/Methylation/methy/gene/meth_gene_overlap.txt",
# #     output:
# #         gene = IN_PATH + "/Methylation/methy/gene/meth_gene_overlap_filt.txt",
# #     params:
# #         bedOverlapRecord = SRC_DIR + "/bedOverlapRecord.py",
# #         ratioThreshold = 0.5,
# #     log:
# #         IN_PATH + "/log/OverlapFilt.log",
# #     run:
# #         shell("python {params.bedOverlapRecord} --bed {input.gene} --ratioThreshold {params.ratioThreshold} --out {output.gene} --method first > {log} 2>&1")

# rule overlapCentromere:
#     input:
#         meth = IN_PATH + "/Methylation/methy/methylation_calls_all_chrs_meth_freq_0.5.bed",
#         centro = IN_PATH + "/Centromere/CEN1600_repeatmask/Vigna_unguiculata_assembly.fasta.out.CENTROMERE.bed",
#     output:
#         centro = IN_PATH + "/Methylation/methy/centromere/meth_centromere_overlap.txt",
#     run:
#         shell("bedtools intersect -a {input.meth} -b {input.centro} -wa -wb > {output.centro}")


# rule overlapGyspy:
#     input:
#         meth = IN_PATH + "/Methylation/methy/methylation_calls_all_chrs_meth_freq_0.5.bed",
#         gyspy = IN_PATH + "/Centromere/CEN1600_repeatmask/Centromere_outside/Vigna_unguiculata_assembly.fasta.masked.out.Gypsy.bed",
#     output:
#         gyspy = IN_PATH + "/Methylation/methy/repeat/meth_gyspy_overlap.txt",
#     run:
#         shell("bedtools intersect -a {input.meth} -b {input.gyspy} -wa -wb > {output.gyspy}")


# rule overlapCopia:
#     input:
#         meth = IN_PATH + "/Methylation/methy/methylation_calls_all_chrs_meth_freq_0.5.bed",
#         copia = IN_PATH + "/Centromere/CEN1600_repeatmask/Centromere_outside/Vigna_unguiculata_assembly.fasta.masked.out.Copia.bed",
#     output:
#         copia = IN_PATH + "/Methylation/methy/repeat/meth_copia_overlap.txt",
#     run:
#         shell("bedtools intersect -a {input.meth} -b {input.copia} -wa -wb > {output.copia}")
        

# ###################################################