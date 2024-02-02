
### environment TE2

############## tandem repeats ###################
rule trf:
    input:
        assembly = IN_PATH + "/Assembly/Scaffold/Fengchan6.scaffold.fasta",
    output:
        mask = IN_PATH + "/Repeats/TRF/Fengchan6.scaffold.fasta.2.7.7.80.10.50.2000.mask",
        dat = IN_PATH + "/Repeats/TRF/Fengchan6.scaffold.fasta.fasta.2.7.7.80.10.50.2000.dat",
    params:
        outDir = IN_PATH + "/Repeats/TRF/",
    threads:
        THREADS
    log:
        IN_PATH + "/log/trf.log",
    run:
        ### 2 7 7 80 10 50 500 -f -d -m
        if not os.path.exists(params.outDir):
            os.makedirs(params.outDir)
        shell("cd {params.outDir} && trf {input.assembly} 2 7 7 80 10 50 2000 -f -d -m -h > {log} 2>&1")

#################################################################

#################### long terminal repeat (LTR) retrotransposons #####
rule ltrfinder:
    input:
        assembly = IN_PATH + "/Assembly/Scaffold/Fengchan6.scaffold.fasta",
    output:
        scn = IN_PATH + "/Repeats/ltrfinder/Fengchan6_ltrfinder.scn",
    threads:
        THREADS
    log:
        IN_PATH + "/log/ltrfinder.log",
    run:
        ## -w 2 -C -D 15000 -d 1000 -L 7000 -l 100 -p 20 -M 0.85
        shell("ltr_finder -w 2 -C -D 15000 -d 1000 -L 7000 -l 100 -p 20 -M 0.85 {input.assembly} > {output.scn} 2>{log}")


rule ltrharvest:
    input:
        assembly = IN_PATH + "/Assembly/Scaffold/Fengchan6.scaffold.fasta",
    output:
        scn = IN_PATH + "/Repeats/ltrharvest/Fengchan6_ltrharvest.scn",
    params:
        index = IN_PATH + "/Repeats/ltrharvest/Fengchan6_ltrharvest",
    threads:
        THREADS
    log:
        IN_PATH + "/log/ltrharvest.log",
    run:
        shell("gt suffixerator -db {input.assembly} -indexname {params.index} -tis -suf -lcp -des -ssp -sds -dna")
        shell("gt ltrharvest -index {params.index} -similar 90 -vic 10 -seed 20 -seqids yes   -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6  -motif TGCA -motifmis 1  > {output.scn} 2>{log}")



rule ltrretriever:
    input:
        assembly = IN_PATH + "/Assembly/Scaffold/Fengchan6.scaffold.fasta",
        ltrfinder = IN_PATH + "/Repeats/ltrfinder/Fengchan6_ltrfinder.scn",
        ltrharvest = IN_PATH + "/Repeats/ltrharvest/Fengchan6_ltrharvest.scn",
    output:
        gff = IN_PATH + "/Repeats/ltrretriever/Fengchan6.scaffold.fasta.LTR.gff3",
        LTRlib = IN_PATH + "/Repeats/ltrretriever/Fengchan6.scaffold.fasta.LTRlib.fa",
        scn = IN_PATH + "/Repeats/ltrretriever/Fengchan6.scaffold.fasta.retriever.all.scn",
    params:
        outDir = IN_PATH + "/Repeats/ltrretriever",
    threads:
        THREADS
    log:
        IN_PATH + "/log/ltrretriever.log",
    run:
        ### output files:
        # Table output for intact LTR-RTs (detailed info)
        #     Vigna_unguiculata_contig.fasta.mod.pass.list (All LTR-RTs)
        #     Vigna_unguiculata_contig.fasta.mod.nmtf.pass.list (Non-TGCA LTR-RTs)
        #     Vigna_unguiculata_contig.fasta.mod.pass.list.gff3 (GFF3 format for intact LTR-RTs)

        # LTR-RT library
        #     Vigna_unguiculata_contig.fasta.mod.LTRlib.redundant.fa (All LTR-RTs with redundancy)
        #     Vigna_unguiculata_contig.fasta.mod.LTRlib.fa (All non-redundant LTR-RTs)
        #     Vigna_unguiculata_contig.fasta.mod.nmtf.LTRlib.fa (Non-TGCA LTR-RTs)

        # Whole-genome LTR-RT annotation by the non-redundant library
        #     Vigna_unguiculata_contig.fasta.mod.LTR.gff3 (GFF3 format)
        #     Vigna_unguiculata_contig.fasta.mod.out.fam.size.list (LTR family summary)
        #     Vigna_unguiculata_contig.fasta.mod.out.superfam.size.list (LTR superfamily summary)

        # LTR Assembly Index (LAI)
        #     Vigna_unguiculata_contig.fasta.mod.out.LAI
        if not os.path.exists(params.outDir):
            os.mkdir(params.outDir)
        shell("cd {params.outDir}  &&  LTR_retriever -genome {input.assembly} -infinder {input.ltrfinder}  -inharvest {input.ltrharvest} -threads {threads} >{log} 2>&1")

##################################################################################





# ################### MITETracker #############
# rule MITETracker:
#     input:
#         assembly = IN_PATH + "/scaffold/final3/Vigna_unguiculata_assembly.fasta",
#     output:
#         faimly = IN_PATH + "/Repeats/MITETracker/results/Vigna/families_nr.fasta",
#         candidate = IN_PATH + "/Repeats/MITETracker/results/Vigna/candidates.fasta",
#     params:
#         MITETracker = "/home/wuzhikun/software/MITE-Tracker/MITETracker.py",
#         outDir = IN_PATH + "/Repeats/MITETracker",
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/MITETracker.log",
#     run:
#         if not os.path.exists(params.outDir):
#             os.mkdir(params.outDir)
#         shell("cd {params.outDir} && python {params.MITETracker} --workers {threads} --jobname Vigna --genome {input.assembly} > {log} 2>&1")


# #############################################



# ############################################################################
# rule BuildDatabase:
#     input:
#         assembly = IN_PATH + "/scaffold/final3/Vigna_unguiculata_assembly.fasta",
#     output:
#         nhr = IN_PATH + "/Repeats/repeatModeler/Vigna_unguiculata_contig.nhr",
#     params:
#         outPrefix = IN_PATH + "/Repeats/repeatModeler/Vigna_unguiculata_contig",
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/BuildDatabase.log",
#     run:
#         shell("BuildDatabase -name {params.outPrefix} -engine ncbi  {input.assembly} > {log} 2>&1")



# rule repeatModeler:
#     ### http://xuzhougeng.top/archives/Repeat-annotation-with-RepeatModeler
#     input:
#         nhr = IN_PATH + "/Repeats/repeatModeler/Vigna_unguiculata_contig.nhr",
#     output:
#         lib = IN_PATH + "/Repeats/repeatModeler/Vigna_unguiculata_contig-families.fa",
#     params:
#         outPrefix = IN_PATH + "/Repeats/repeatModeler/Vigna_unguiculata_contig",
#         outDir = IN_PATH + "/Repeats/repeatModeler",
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/repeatModeler.log",
#     run:
#         if not os.path.exists(params.outDir):
#             os.mkdir(params.outDir)
#         shell("cd {params.outDir} &&  RepeatModeler -pa {threads} -database {params.outPrefix}  -engine ncbi > {log} 2>&1")
#         # output file:
#         # unaligned.fa
#         # RM_113399.ThuJul81254042021
#         #     ├── consensi.fa
#         #     ├── consensi.fa.classified
#         #     ├── consensi.fa.masked
#         #     ├── families-classified.stk
#         #     ├── families.stk




# rule repeatMasker:
#     input:
#         assembly = IN_PATH + "/scaffold/final3/Vigna_unguiculata_assembly.fasta",
#         lib = IN_PATH + "/Repeats/repeatModeler/Vigna_unguiculata_contig-families.fa",
#     output:
#         # log = IN_PATH + "/log/repeatMasker.log",
#         mask = IN_PATH + "/Repeats/repeatModeler/Vigna_unguiculata_contig/Vigna_unguiculata_assembly.fasta.masked",
#         out = IN_PATH + "/Repeats/repeatModeler/Vigna_unguiculata_contig/Vigna_unguiculata_assembly.fasta.out",
#     threads:
#         THREADS
#     params:
#         outPrefix = IN_PATH + "/Repeats/repeatModeler/Vigna_unguiculata_contig",
#     log:
#         IN_PATH + "/log/repeatMasker.log",
#     run:
#         # -libdir <string>
#         #     Path to the RepeatMasker libraries directory.
#         shell("RepeatMasker -pa {threads} -e ncbi  -gff -lib {input.lib}  -dir {params.outPrefix}  {input.assembly} > {log} 2>&1")
    



# rule repeatMasker2:
#     input:
#         assembly = IN_PATH + "/scaffold/final3/Vigna_unguiculata_assembly.fasta",
#         lib = IN_PATH + "/Repeats/repeatModeler/Vigna_unguiculata_contig-families.fa",
#         lib2 = IN_PATH + "/Repeats/ltrretriever/Vigna_unguiculata_assembly.fasta.LTRlib.fa",
#     output:
#         # log = IN_PATH + "/log/repeatMasker.log",
#         outDir = IN_PATH + "/Repeats/repeatMasker/Vigna_unguiculata_contig.fasta.out",
#     threads:
#         THREADS
#     params:
#         outPrefix = IN_PATH + "/Repeats/repeatMasker/Vigna_unguiculata_contig",
#     log:
#         IN_PATH + "/log/repeatMasker2.log",
#     run:
#         shell("RepeatMasker -pa {threads} -e ncbi  -gff -lib {input.lib}  -dir {params.outPrefix}  {input.assembly} > {log} 2>&1")


# rule repeatWin:
#     input:
#         repeat = IN_PATH + "/Repeats/repeatModeler/Vigna_unguiculata_contig/Vigna_unguiculata_assembly.fasta.out",
#     output:
#         win50 = IN_PATH + "/Repeats/repeatModeler/Vigna_unguiculata_contig/Vigna_unguiculata_assembly.repeat_window50k.txt",
#         win10 = IN_PATH + "/Repeats/repeatModeler/Vigna_unguiculata_contig/Vigna_unguiculata_assembly.repeat_window10k.txt",
#     params:
#         WindowTRNumber = SRC_DIR + "/WindowTRNumber.py",
#         genomeLen = IN_PATH +  "/scaffold/final3/Vigna_unguiculata_assembly.fasta.split/Vigna_unguiculata_assembly.id.chr_length.txt",
#         # window = 100000,
#         # sliding = 0,
#     log:
#         IN_PATH + "/log/repeatWin.log",
#     run:
#         shell("python {params.WindowTRNumber} --genome {params.genomeLen} --input {input.repeat} --out {output.win50} --window 50000 --sliding 0 > {log} 2>&1")
#         shell("python {params.WindowTRNumber} --genome {params.genomeLen} --input {input.repeat} --out {output.win10} --window 10000 --sliding 0 > {log} 2>&1")




# ##################################################







# ################### repeatscount #############
# rule buildTable:
#     input:
#         assembly = IN_PATH + "/scaffold/final3/Vigna_unguiculata_assembly.fasta",
#     output:
#         freq = IN_PATH + "/Repeats/repeatscount/Vigna_unguiculata_contig_repeatscount_ferquency.txt",
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/buildTable.log",
#     run:
#         shell("build_lmer_table -sequence {input.assembly} -freq {output.freq} > {log} 2>&1")

# rule RepeatScout:
#     input:
#         assembly = IN_PATH + "/scaffold/final3/Vigna_unguiculata_assembly.fasta",
#         freq = IN_PATH + "/Repeats/repeatscount/Vigna_unguiculata_contig_repeatscount_ferquency.txt",
#     output:
#         fa = IN_PATH + "/Repeats/repeatscount/Vigna_unguiculata_contig_repeatscount.fasta",
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/RepeatScout.log",
#     run:
#         shell("RepeatScout -sequence {input.assembly} -freq {input.freq} -output {output.fa} > {log} 2>&1")


# rule RepeatFilt:
#     input:
#         fa = IN_PATH + "/Repeats/repeatscount/Vigna_unguiculata_contig_repeatscount.fasta",
#     output:
#         fa = IN_PATH + "/Repeats/repeatscount/Vigna_unguiculata_contig_repeatscount_filt.fasta",    
#     threads:
#         THREADS
#     log:
#         IN_PATH + "/log/RepeatFilt.log",
#     run:
#         ### /home/wuzhikun/anaconda3/envs/Assembly3/bin/filter-stage-1.prl
#         shell("filter-stage-1.prl  {input.fa} > {output.fa} 2>{log}")


# rule repeatMaskerScout:
#     input:
#         assembly = IN_PATH + "/scaffold/final3/Vigna_unguiculata_assembly.fasta",
#         lib = IN_PATH + "/Repeats/repeatscount/Vigna_unguiculata_contig_repeatscount_filt.fasta",
#     output:
#         # log = IN_PATH + "/log/repeatMasker.log",
#         outDir = IN_PATH + "/Repeats/repeatscount/Vigna_unguiculata_contig.fasta.out",
#     threads:
#         THREADS
#     params:
#         outPrefix = IN_PATH + "/Repeats/repeatscount/Vigna_unguiculata_contig",
#     log:
#         IN_PATH + "/log/repeatMaskerScout.log",
#     run:
#         shell("RepeatMasker -pa {threads} -e ncbi  -gff -lib {input.lib}  -dir {params.outPrefix}  {input.assembly} > {log} 2>&1")

# ##############################################





















# ######################################################################
# rule rfamInfernal:
#     input:
#         assembly = IN_PATH + "/scaffold/final3/Vigna_unguiculata_assembly.fasta",
#     output:
#         table = IN_PATH + "/Feature/Rfam/Vigna_unguiculata_Rfam_table.txt",
#         hit = IN_PATH + "/Feature/Rfam/Vigna_unguiculata_Rfam_predict.txt",
#     params:
#         Rfam = "/home/wuzhikun/database/Rfam/Rfam.cm"
#     threads:
#         THREADS 
#     log:
#         IN_PATH + '/log/rfam.log' 
#     run:
#         ### infernal
#         # shell("cmpress {params.Rfam}")
#         shell("cmscan --cpu {threads}  --tblout {output.table}  -o {output.hit}   {params.Rfam} {input.assembly} > {log} 2>&1")


# rule tRNAscan:
#     input:
#         assembly = IN_PATH + "/scaffold/final3/Vigna_unguiculata_assembly.fasta",
#     output:
#         table = IN_PATH + "/Feature/tRNAscan/Vigna_unguiculata_tRNAscan_final.txt",
#         stats = IN_PATH + "/Feature/tRNAscan/Vigna_unguiculata_tRNAscan_stats.txt",
#         structure = IN_PATH + "/Feature/tRNAscan/Vigna_unguiculata_tRNAscan_structure.txt",
#         bed = IN_PATH + "/Feature/tRNAscan/Vigna_unguiculata_tRNAscan_bed.txt",
#     threads:
#         THREADS 
#     log:
#         IN_PATH + '/log/tRNAscan.log' 
#     run:
#         ### -g ~/soft/miniconda3/pkgs/trnascan-se-2.0.6-pl526h516909a_0/lib/tRNAscan-SE/gcode/gcode.vertmito
#         shell("tRNAscan-SE -qQ --detail --thread {threads} -o {output.table} -m {output.stats} -f {output.structure} -b {output.bed}   {input.assembly} > {log} 2>&1")
# ############################################################################



# ############ RNAmmer ################
# rule RNAmmer:
#     input:
#         assembly = IN_PATH + "/scaffold/final3/Vigna_unguiculata_assembly.fasta",
#     output:
#         gff = IN_PATH + "/Feature/RNAmmer/RNAmmer_rna.gff",
#     params:
#         outDir = IN_PATH + "/Feature/RNAmmer",
#     threads:
#         THREADS 
#     log:
#         IN_PATH + '/log/RNAmmer.log' 
#     run:
#         ### /home/wuzhikun/software/rnammer-1.2/rnammer  -S euk -m tsu,ssu,lsu  -gff RNAmmer_rna.gff  -f RNAmmer_rna.fasta -h RNAmmer_report.html    /home/wuzhikun/Project/Vigna/Assembly/Polish/Racon/Vigna_unguiculata_ONT_polish_racon3.fasta
#         if not os.path.exists(params.outDir):
#             os.mkdir(params.outDir)
#         shell("cd {params.outDir} &&  rnammer  -S euk -m tsu,ssu,lsu  -gff RNAmmer_rna.gff  -f RNAmmer_rna.fasta -h RNAmmer_report.html   {input.assembly} > {log} 2>&1")

# ######################################
