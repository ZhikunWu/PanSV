### [Consolidate GVCFs for joint calling with GenotypeGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/360035889971--How-to-Consolidate-GVCFs-for-joint-calling-with-GenotypeGVCFs)

```
$ java -jar /home/wuzhikun/software/gatk-4.0.6.0/gatk-package-4.0.6.0-local.jar  GenomicsDBImport   -V /home/wuzhikun/Project/PanSV/NGS/VCF/E918/E918_NGS.g.vcf -V /home/wuzhikun/Project/PanSV/NGS/VCF/E919/E919_NGS.g.vcf  --genomicsdb-workspace-path /home/wuzhikun/Project/PanSV/NGS/genomicsdb  -L Chr11
```

out files:
```
$ lk /home/wuzhikun/Project/PanSV/NGS/genomicsdb
total 116K
-rwx------ 1 wuzhikun wuzhikun    0 Nov  1 15:33 __tiledb_workspace.tdb
-rw-rw-r-- 1 wuzhikun wuzhikun  182 Nov  1 15:33 callset.json
-rw-rw-r-- 1 wuzhikun wuzhikun  68K Nov  1 15:33 vidmap.json
-rw-rw-r-- 1 wuzhikun wuzhikun  47K Nov  1 15:33 vcfheader.vcf
drwx------ 1 wuzhikun wuzhikun 4.0K Nov  1 15:50 Chr11$1$48714482

```



```
java -jar /home/wuzhikun/software/gatk-4.0.6.0/gatk-package-4.0.6.0-local.jar  GenotypeGVCFs  -R /home/wuzhikun/database/genome/Vigna_unguiculata/G98/Lachesis_assembly_changed.fa   -V gendb:///home/wuzhikun/Project/PanSV/NGS/genomicsdb  -O /home/wuzhikun/Project/PanSV/NGS/Chr11_joint_SNV.vcf
```



