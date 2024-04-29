
```

cd-hit -i Faba_10genome_protein.fasta  -o Faba_10genome_protein_cdhit  -c 0.9 -s 0.9  -n 5 -M 100000 -d 0 -T 20


cd-hit -i Faba_10genome_protein.fasta  -o Faba_10genome_protein_cdhit  -c 0.8 -s 0.8  -n 5 -M 200000 -d 0 -T 40



miniprot -t 20  --gff /home/wuzhikun/Project/PanSV/BACKUP/Assembly/Fengchan6/Vigna_unguiculata_assembly.fasta  /home/wuzhikun/data/Vigna/Protein/Faba_10genome_protein_cdhit_0.8 > aln.gff


```
