


rule orthoStat:
    input:
        stat = "/bioinfor/mulinlab/wuzhikun/Project/Evolution/orthology/orthoResult/Results_Aug28/Orthogroups/ortho_cat_stats.txt",
    output:
        cate = "/bioinfor/mulinlab/wuzhikun/Project/Evolution/orthology/orthoResult/Results_Aug28/Orthogroups/ortho_category.txt",
        stat = "/bioinfor/mulinlab/wuzhikun/Project/Evolution/orthology/orthoResult/Results_Aug28/Orthogroups/ortho_cat_stats.txt",
    params:
        PangeneClass = "/bioinfor/mulinlab/wuzhikun/github/PanSV/script/PangeneClass.py",
    run:
        shell("python {params.PangeneClass} --input {input.stat} --category {output.cate} --stats {output.stat}")



rule categoryNum:
    input:
        cate = "/bioinfor/mulinlab/wuzhikun/Project/Evolution/orthology/orthoResult/Results_Aug28/Orthogroups/ortho_category.txt",
        stat = "/bioinfor/mulinlab/wuzhikun/Project/Evolution/orthology/orthoResult/Results_Aug28/Orthogroups/ortho_cat_stats.txt",
    output:
        panstat = "/bioinfor/mulinlab/wuzhikun/Project/Evolution/orthology/orthoResult/Results_Aug28/Orthogroups/sample_ortho_number.xls",
    params:
        SampleCateNumber = "/bioinfor/mulinlab/wuzhikun/github/PanSV/script/SampleCateNumber.py",
    run:
        shell("python {params.SampleCateNumber} --category {input.cate} --orthogroup {input.stat} --out {output.panstat}")




rule orthogroupBar:
    input:
        panstat = "/bioinfor/mulinlab/wuzhikun/Project/Evolution/orthology/orthoResult/Results_Aug28/Orthogroups/sample_ortho_number.xls",
    output:
        pd =  "/bioinfor/mulinlab/wuzhikun/Project/Evolution/orthology/orthoResult/Results_Aug28/Orthogroups/sample_ortho_number.pdf",
    params:
        PangeneBar = "/bioinfor/mulinlab/wuzhikun/github/PanSV/script/PangeneBar.R",
    run:
        shell("source activate Rmeta && Rscript {params.PangeneBar} --input {input.panstat} --pdf {output.pdf} --width 5 --height 4")




rule orthogroupPie:
    input:
        stat = "/bioinfor/mulinlab/wuzhikun/Project/Evolution/orthology/orthoResult/Results_Aug28/Orthogroups/ortho_cat_stats.txt",
    output:
        pdf =  "/bioinfor/mulinlab/wuzhikun/Project/Evolution/orthology/orthoResult/Results_Aug28/Orthogroups/ortho_cat_stats.pdf",
    params:
        orthoGroupPie = "/bioinfor/mulinlab/wuzhikun/github/PanSV/script/OrthoGroupPie.R",
    run:
        shell("source activate Rmeta && Rscript {params.orthoGroupPie}  --input {input.stat} --pdf {output.pdf} --width 5 --height 4")


