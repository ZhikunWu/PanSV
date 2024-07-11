rule GeneRename:
    input:
        pep = IN_PATH + "/GenePrediction/evidencemodeler/{sample}/{sample}.EVM.pep",
        gff = IN_PATH + "/GenePrediction/evidencemodeler/{sample}/{sample}.EVM.gff3",
        cds = IN_PATH + "/GenePrediction/evidencemodeler/{sample}/{sample}.EVM.cds",
    output:
        pep = IN_PATH + "/GenePrediction/Gene/{sample}/{sample}.pep.fasta",
        gff = IN_PATH + "/GenePrediction/Gene/{sample}/{sample}.gene.gff3",
        cds = IN_PATH + "/GenePrediction/Gene/{sample}/{sample}.cds.fasta",
    run:
        shell("sed 's/evm.model/{wildcards.sample}/g' {input.pep}  | sed 's/evm.TU/{wildcards.sample}/g'  >  {output.pep}")
        shell("sed 's/evm.model/{wildcards.sample}/g' {input.cds}  | sed 's/evm.TU/{wildcards.sample}/g'  >  {output.cds}")
        shell("sed 's/evm.model/{wildcards.sample}/g' {input.gff}  | sed 's/evm.TU/{wildcards.sample}/g'  >  {output.gff}")






rule ArabidopsisHomo:
    input:
        fa = IN_PATH + "/GenePrediction/Gene/{sample}/{sample}.pep.fasta",
    output:
        ara = IN_PATH + "/Protein/arabidopsis/{sample}_arabidopsis_diamond_blast.txt",
    threads:
        THREADS
    params:
        ara_dmnd = "/home/wuzhikun/database/genome/arabidopsis/TAIR10_pep_20101214.dmnd",
    log:
        IN_PATH + "/log/AThomo_{sample}.log"
    run:
        shell("diamond blastp --query {input.fa} --db {params.ara_dmnd} --out {output.ara} --outfmt 6 --threads {threads} > {log} 2>&1")






rule NR:
    input:
        fa = IN_PATH + "/GenePrediction/Gene/{sample}/{sample}.pep.fasta",
    output:
        nr = IN_PATH + "/Protein/NR/{sample}_nr_diamond_blast.txt",
    threads:
        THREADS
    params:
        # nr = "/home/wuzhikun/database/NCBI/nr2020/nr",
        nr_dmnd = "/home/wuzhikun/database/NCBI/NR/nr.dmnd",
    log:
        IN_PATH + "/log/nr_blast_{sample}.log"
    run:
        ### -max_target_seqs 1
        # shell("blastp -query {input.fa} -db {params.nr} -out {output.blast} -num_threads {threads}  -evalue 1e-5 -outfmt 6 > {log} 2>&1")
        shell("diamond blastp --query {input.fa} --db {params.nr_dmnd} --out {output.nr} --outfmt 6 --threads {threads} > {log} 2>&1")



rule swissProt:
    ### UniProtKB/Swiss-Prot  which is manually annotated and is reviewed and
    ### UniProtKB/TrEMBL  which is automatically annotated and is not reviewed.
    input:
        fa = IN_PATH + "/GenePrediction/Gene/{sample}/{sample}.pep.fasta",
    output:
        trembl = IN_PATH + "/Protein/SwissProt/{sample}_trembl_diamond_blast.txt",
        swiss = IN_PATH + "/Protein/SwissProt/{sample}_SwissProt_diamond_blast.txt",
    threads:
        THREADS
    params:
        # swiss = "/home/wuzhikun/database/uniprot/Swiss_Prot_review/swiss_prot.fasta",
        # swissPre = "/home/wuzhikun/database/uniprot/Swiss_Prot_review/swiss_prot",
        trembl_dmnd = "/home/wuzhikun/database/uniprot/uniprot_trembl.dmnd",
        swiss_dmnd = "/home/wuzhikun/database/uniprot/Swiss_Prot_review/swiss_prot.dmnd",
    log:
        IN_PATH + "/log/swissProt_{sample}.log"
    run:
        ### makeblastdb -in /home/wuzhikun/database/uniprot/Swiss_Prot_review/swiss_prot.fasta -out /home/wuzhikun/database/uniprot/Swiss_Prot_review/swiss_prot -dbtype prot
        # shell("makeblastdb -in {params.swiss} -out {params.swissPre} -dbtype prot > {log} 2>&1")
        ### diamond blastp --query /home/wuzhikun/Project/Vigna/RNA/NGS/Transcript/EVidenceModeler_1231/result/EVidenceModeler_out.fasta  --db /home/wuzhikun/database/uniprot/uniprot_trembl.dmnd --out /home/wuzhikun/Project/Vigna/RNA/NGS/Transcript/EVidenceModeler_1231/result/EVidenceModeler_out_nr_diamond.txt  --outfmt 6 --threads 20
        #shell("blastp -query {input.fa} -db {params.swissPre} -out {output.blast} -num_threads {threads}  -evalue 1e-5 -outfmt 6 > {log} 2>&1")
        shell("diamond blastp --query {input.fa} --db {params.trembl_dmnd} --out {output.trembl} --outfmt 6 --threads {threads} > {log} 2>&1")
        shell("diamond blastp --query {input.fa} --db {params.swiss_dmnd} --out {output.swiss} --outfmt 6 --threads {threads} 2>>{log}")


rule pfam:
    input:
        fa = IN_PATH + "/GenePrediction/Gene/{sample}/{sample}.pep.fasta",
    output:
        hmm = IN_PATH + "/Protein/Pfam/{sample}_Pfam_hmm_out.txt",
        seq = IN_PATH + "/Protein/Pfam/{sample}_Pfam_hmm_seq.txt",
        domain = IN_PATH + "/Protein/Pfam/{sample}_Pfam_hmm_domain.txt",
        seqDomain = IN_PATH + "/Protein/Pfam/{sample}_Pfam_hmm_seq_domain.txt",
    params:
        hmm = "/home/wuzhikun/database/Pfam/Pfam-A.hmm",
    threads:
        THREADS
    log:
        IN_PATH + "/log/pfam_{sample}.log"
    run:
        ### hmmsearch: multiple threads
        ### index
        # shell("hmmpress {params.hmm} > {log}")
        # shell("hmmscan --cpu {threads} -o {output.hmm} --tblout {output.seq} --domtblout {output.domain} --pfamtblout {output.seqDomain} {params.hmm} {input.fa} > {log} 2>&1")
        shell("hmmsearch  --cpu {threads} --tblout {output.seq} --domtblout {output.domain} --pfamtblout {output.seqDomain}  -o {output.hmm} {params.hmm} {input.fa} > {log} 2>&1")



rule bestMatch:
    input:
        ara = IN_PATH + "/Protein/arabidopsis/{sample}_arabidopsis_diamond_blast.txt",
        nr = IN_PATH + "/Protein/NR/{sample}_nr_diamond_blast.txt",
        trembl = IN_PATH + "/Protein/SwissProt/{sample}_trembl_diamond_blast.txt",
        swiss = IN_PATH + "/Protein/SwissProt/{sample}_SwissProt_diamond_blast.txt",
    output:
        ara = IN_PATH + "/Protein/arabidopsis/{sample}_arabidopsis_diamond_blast_best.txt",
        nr = IN_PATH + "/Protein/NR/{sample}_nr_diamond_blast_best.txt",
        trembl = IN_PATH + "/Protein/SwissProt/{sample}_trembl_diamond_blast_best.txt",
        swiss = IN_PATH + "/Protein/SwissProt/{sample}_SwissProt_diamond_blast_best.txt",
    params:
        blastBestMatch = SCRIPT_DIR + "/blastBestMatch.py",
    log:
        IN_PATH + "/log/bestMatch_{sample}.log"
    run:
        shell("python {params.blastBestMatch} --blast {input.ara} --out {output.ara} > {log} 2>&1")
        shell("python {params.blastBestMatch} --blast {input.nr} --out {output.nr} > {log} 2>&1")
        shell("python {params.blastBestMatch} --blast {input.trembl} --out {output.trembl} > {log} 2>&1")
        shell("python {params.blastBestMatch} --blast {input.swiss} --out {output.swiss} > {log} 2>&1")



rule interproscan:
    input:
        fa = IN_PATH + "/GenePrediction/Gene/{sample}/{sample}.pep.fasta",
    output:
        fa = temp(IN_PATH + "/GenePrediction/Gene/{sample}/{sample}.pep_temp.fasta"),
        interpro = IN_PATH + "/Protein/interproscan/{sample}_protein.iprscan",
    params:
        interproscan = "/home/wuzhikun/anaconda3/envs/PanSV/bin/interproscan.sh",
    threads:
        THREADS
    log:
        IN_PATH + "/log/interproscan_{sample}.log"
    run:
        shell("sed  's/*//g' {input.fa} > {output.fa}")
        shell("{params.interproscan} -cpu {threads}  -dp -f TSV -goterms -iprlookup -pa -t p -i {output.fa} -o {output.interpro} > {log} 2>&1")




rule pfamLRR:
    input:
        fa = IN_PATH + "/GenePrediction/Gene/{sample}/{sample}.pep.fasta",
    output:
        hmm = IN_PATH + "/Protein/NBS_LRR/{sample}_Pfam_hmm_out.txt",
        seq = IN_PATH + "/Protein/NBS_LRR/{sample}_Pfam_hmm_seq.txt",
        domain = IN_PATH + "/Protein/NBS_LRR/{sample}_Pfam_hmm_domain.txt",
        seqDomain = IN_PATH + "/Protein/NBS_LRR/{sample}_Pfam_hmm_seq_domain.txt",
    params:
        hmm = "/home/wuzhikun/database/Pfam/NB-ARC.hmm",
    threads:
        THREADS
    log:
        IN_PATH + "/log/NBS_LRR_{sample}.log"
    run:
        ### hmmsearch: multiple threads
        ### index
        # shell("hmmpress {params.hmm} > {log}")
        # shell("hmmscan --cpu {threads} -o {output.hmm} --tblout {output.seq} --domtblout {output.domain} --pfamtblout {output.seqDomain} {params.hmm} {input.fa} > {log} 2>&1")
        shell("hmmsearch  --incE 0.01 -E 0.01 --cpu {threads} --tblout {output.seq} --domtblout {output.domain} --pfamtblout {output.seqDomain}  -o {output.hmm} {params.hmm} {input.fa} > {log} 2>&1")
        ### seqkit grep -f NBS_LRR_Vu_gene.txt  /home/wuzhikun/Project/Vigna/GenePrediction/gene/Vigna_unguiculata_assembly.all.protein.fasta > NBS_LRR_Vu_gene_seq.fa



rule LRRGene:
    input:
        fa = IN_PATH + "/GenePrediction/Gene/{sample}/{sample}.pep.fasta",
        domain = IN_PATH + "/Protein/NBS_LRR/{sample}_Pfam_hmm_seq.txt",
    output:
        gene = IN_PATH + "/Protein/NBS_LRR/{sample}_NBS_LRR_hmm_gene.txt",
        fa = IN_PATH + "/Protein/NBS_LRR/{sample}_NBS_LRR_hmm_gene.fa",
    run:
        shell("grep -v '^#' {input.domain} | cut -f 1 -d ' ' | sort | uniq > {output.gene}")
        shell("seqkit grep -f {output.gene} {input.fa} > {output.fa} ")


