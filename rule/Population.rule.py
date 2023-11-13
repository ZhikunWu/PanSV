rule PopStat:
    input:
        vcf = IN_PATH + "/SVCall/Sniffles2/Samples_SV_merge_filt.vcf",
    output:
        stat1 = IN_PATH + "/Population/Stats/Samples_SV_stat1.txt",
        stat2 = IN_PATH + "/Population/Stats/Samples_SV_stat2.txt",
    params:
        PopGenoStats = SCRIPT_DIR + "/PopGenoStats.py",
    log:
        IN_PATH + "/log/PopStat.log",
    run:
        shell("python {params.PopGenoStats} --vcf {input.vcf} --stat1 {output.stat1} --stat2 {output.stat2} > {log} 2>&1")
