
### A telomere-to-telomere genome assembly of Zhonghuang 13, a widely-grown soybean variety from the original center of Glycine max


### The genome of a wild Medicago species provides insights into the tolerant mechanisms of legume forage to environmental stress

De novo predictions, homolog-based, and RNA-seq-based predictions were employed to annotate the protein-coding genes. Five ab initio gene prediction programs were used to predict genes, including Augustus v3.0.2, Genescan v1.0, Geneid, GlimmerHMM v3.0.2, and SNAP. Protein sequences of ten homologous species A. thaliana, M. truncatula, G. max, T. pratense, C. arietinum, Vigna angularis, V. radiate, C. cajan, P. vulgaris, and Arachis ipaensis were downloaded from the Ensembl or NCBI. Homologous sequences were aligned against to the repeat-masked M. ruthenica genome using TblastN (e value ≤1e−5). Genewise v2.2.0 was employed to predict gene models based on the alignment sequences [49]. There were two ways to assemble the RNA-seq data into the unique sequences of transcripts. One was mapping the RNA-seq data to the M. ruthenica genome using Tophat v2.0.8 [50] and Cufflinks v2.1.1 [51] for transcript assembly. The other was applying Trinity [52] to assemble the RNA-seq data, and then PASA software [53] improved the gene structures. A weighted and non-redundant gene set was generated by EVidenceModeler [54] which merged all genes models predicted by the above three approaches. PASA adjusted the gene models generated by EVM. Lastly, the gene sets were filtered according to the following standards: coding region lengths of amino acids ≤50, supported only by de novo methods and with FPKM<5.


Functional annotation of protein-coding genes was obtained according to best Blast hit by BlastP (e value ≤1e−5) against SwissProt [55] and NCBI non-redundant (NR) protein databases. Motifs and domains were annotated by using InterProScan v4.7 [56] to search against InterPro v29.0 databases [56], including Pfam, Prints, Prosite, ProDom, and Smart. The tRNA genes were predicted by tRNAscan-SE software [57]. The miRNA and snRNA fragments were identified by INFERNAL software [58]against the Rfam database [59]. The rRNAs were found by using BlastN (e value ≤1e−10) against invertebrate rRNA database. The structure figure was drawn along eight chromosomes of M. ruthenica genome using Circos program [60]. To estimate the assembly of genome, transcriptome data. from roots, stems, leaves, flowers, and pods were mapped to the genome assembly using Hisat2 [61], and transcripts were assembled using Trinity [52]. The unigene, which was the longest transcript selected from Trinity, was aligned to the genome assembly by Blat [62].





### Chromosome-scale assembly of the African yam bean genome
Gene prediction and functional annotation of genome We combined transcript and homology evidence to annotate the gene content of the AYB genome. The transcript evidence was generated from 17,117,377 ONT longread RNA reads totalling 7.1 Gbp of sequencing data used for de novo assembly of 60,249 transcripts. Briefly, Minimap239 (v2.22) was used to index the AYB genome assembly and the RNA reads were mapped to the indexed assembly. Samtools16 (v1.9) was used to sort mapped reads by coordinates that were used to assemble transcripts with Stringtie240 (v2.0.1). Transdecoder41 (v2.0.1) was then used to identify candidate CDS regions and select transcripts with a minimum protein length of 100 amino acids.
We combined the de novo transcripts with protein homology evidence from four wellannotated plant genomes (Arabidopsis thaliana TAIR10, Phaseolus vulgaris v1.0, Glycine max v2.1, Vigna_angularis v1.1) together with a soft-masked (for repeats) AYB genome as inputs into Funannotate42 (v1.8.11) to identify protein coding genes. Funannotate ‘predict’ uses ab initio gene predictors Augustus43, PASA44, SNAP45 and GlimmerHMM46 together with protein sequences as evidence to predict genes. Gene predictions from all four ab initio predictors are passed to EVidenceModeler47 with various weights for integration. This resulted in 30,840 coding gene models totalling 31,614 transcripts with a median exon length of 231 bp and a median of three exons per transcript. Additionally, we detected 774 non-overlapping tRNAencoding
genes using tRNAscan-SE48 for tRNA prediction. The gene and
transposable element distribution across the genome are inversely correlated (Fig. 5).
Protein domains were annotated using InterProScan-5.25-64.049 based on InterPro protein databases, including TIGRFAM, SUPERFAMILY, PANTHER, Pfam, PRINTS and ProDom. We also used eggNOG-mapper50 (v2.1.7) to annotate predicted gene models. Funannotate ‘annotate’ uses results of InterProScan and eggNOG-mapper to annotate putative functions of protein sequences using PFAM51, UniProtKB52 and Gene Ontologies53 databases. In total, 25,241 (81.85%) of the genes.

Gene family analysis To delineate gene families, AYB was compared to other legumes L. purpureus, P. vulgaris, Vigna unguiculata, and Mycrotyloma uniflorum (with Solanum tuberosum as an outgroup) using OrthoFinder54 (v2.5.4). This analysis placed 26,038 (84.4%) of the 30,840 AYB proteins into orthogroups. Clustering using Venn Diagrams55 revealed 1,296 AYB proteins (4.2%) segregated in 384 species-specific orthogroups (Fig. 7).


### Chromosome-level genome assembly of bean flower thrips Megalurothrips usitatus
Gene and functional predictions. 
Genes in the assembled genome were predicted using a combination of homology-based, transcriptome-based, and ab initio methods. Homology-based predictions involved downloaded sequences of peptides and transcripts from Aptinothrips rufus (http://v2.insect-genome.com/ Organism/87), Frankliniella occidentalis (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/697/945/ GCF_000697945.3_Focc_3.1), and Thrips palmi (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/012/932/325/ GCF_012932325.1_TpBJ-2018v1). The IsoSeq version 3.4.0 workflow was utilized to generate 28,608 high-quality transcripts from the full-length transcriptome data, with quality parameters of 0.99 (https://github.com/ PacificBiosciences/IsoSeq). Next, RNA-seq short data were mapped to the reference genome using HISAT2 version 2.2.134 with the parameter ‘-k 2’. The mapped reads were then assembled into transcripts using StringTie version 2.4.035 with default parameters. Homologous proteins and transcripts were aligned using Exonerate version 2.4.0 with default parameters to train the gene sets. Meanwhile, a sorted and mapped bam file of RNA-seq data was transferred to a hints file using the bam2hints program in AUGUSTUS version 3.2.336 with the parameter ‘–intronsonly’. The trained gene sets and hint files were combined as inputs for AUGUSTUS version 3.2.336, which predicted coding genes from the assembled genome with default parameters. Finally, homology-based, de novo-derived, and transcript genes were merged in MAKER version 2.31.10 to generate a high-confidence gene set37. It resulted in the annotation of 14,000 M. usitatus genes. The average transcript length was 2,243.30 bp with an average length of coding sequence (CDS) of 1,588.94 bp. The average exon number per gene was 7.38, and the average exon length was 303.85 bp (Table 4).




## Genome assemblies of Vigna reflexo-pilosa (créole bean) and its progenitors, Vigna hirtella and Vigna trinervia , revealed homoeolog expression bias and expression-level dominance in the allotetraploid

To identify re petiti ve element families in the genome assembly, RepeatModeler version 2.0.3 ( http://www.repeatmasker.org/Rep eatModeler/) was used to construct a de novo r epeat libr ary [ 33 ]. This pipeline emplo y ed 2 distinct repeat discovery algorithms, RECON (version 1.08) and RepeatScout (version 1.0.6), to identify the boundaries of re petiti ve elements and build consensus models of interspersed repeats [ 34 , 35 ]. We aligned repeat sequences in the library to GenBank’s nr protein database using BLASTX (e-value cutoff = 10 −6 ) to ensure that they did not contain large families of protein-coding sequences.
To identify protein-coding sequences in the unmasked assembly, we used EvidenceModeler (EVM) version 1.1.1 to combine evidence from RNA-based prediction, homology-based prediction, and ab initio prediction [ 36 ]. For RNA-based prediction, we used e vidence fr om P acBio Iso-seq data obtained fr om leaf, stem, and flo w er tissues. Full-length transcripts were mapped to the final assembly using the genomic mapping and alignment pr ogr am (GMAP) version 2020–09-12 [ 37 ]. Protein sequences from Glycine max , Phaseolus vulgaris , Vigna unguiculata , Vigna angularis , Vigna mungo , and Arabidopsis thaliana available on the public databases were aligned to the unmasked genome using analysis and annotation tool (AAT) [ 38 ]. Protein-coding gene predictions were obtained with Augustus version 3.2.1 [ 39 ] trained with G. max , P. vulgaris , V. unguiculata , V. angularis , V. mungo , V. reflexo-pilosa , and A. thaliana PASA transcriptome alignment assembly using V. reflexopilosa alignment files as inputs. All gene pr edictions wer e integrated by EVM to generate consensus gene models using the follo wing w eights for eac h e vidence type: PASA2–1, GMAP–0.5, AAT–0.3, and Augustus–0.3. Any predicted genes that had more than 20% ov erla pping sequence with r e petiti ve sequences or had no RNA-seq support were excluded from the list of annotated genes.



## Genome assemblies of Vigna reflexo-pilosa (créole bean) and its progenitors, Vigna hirtella and Vigna trinervia , revealed homoeolog expression bias and expression-level dominance in the allotetraploid
### Phylogenetic analyses and comparati v e genomics 
OrthoFinder version 2.4.0 [ 40 ] was used to identify orthologous groups in A. thaliana , Citrullus lanatus , Cucumis melo , Cucumis sativus , G. max , Oryza sativa , P. vulgaris , V. hirtella , V. mungo , Vigna radiata , V. reflexo-pilosa , V. trinervia , and V. unguiculata . We constructed a phylogenetic tree based on protein sequences from single-copy orthologous groups using RAxML-NG pr ogr am v ersion 1.0.2 [ 41 ]. We first aligned protein sequences in each single-copy orthologous group with MUSCLE version 3.8.1551 [ 42 ] and r emov ed alignment gaps with trimAI version 1.4 rev15 [ 43 ] using the automated1 heuristic method. We subsequently concatenated alignment blocks using the catsequences pr ogr am ( https: //github.com/Chr isCr eevey/catsequences ), and the substitution model for each block was estimated using the ModelTest-NG progr am v ersion 0.1.7 [ 44 ]. The outputs wer e used to compute a maximum likelihood phylogenetic tree. Divergence times were estimated using the MCMCtree software version 4.0 (PAML 4 package) [ 45 ] using the r elaxed-cloc k model with the known div er gence time between C. melo and C. sativus , which was estimated at 8.4 to 11.8 million years ago (MYA) [ 46 ].

### Genome synteny analysis 
MCscanX [ 47 ] was used to analyze the colinearity within the V. reflexo-pilosa genome and between V. reflexo-pilosa–G. max , V. reflexo-pilosa–P. vulgaris , V. reflexo-pilosa–V. hirtella , V. reflexo-pilosa– V . mungo , V . reflexo-pilosa–V . radiata , V . reflexo-pilosa–V . trinervia , and V . reflexo-pilosa–V . unguiculata genomes. V . reflexo-pilosa amino acid sequences were aligned against themselves , G . max , P. vulgaris , V . hirtella , V . mungo , V . radiata , V . trinervia , or V . unguiculata using BLASTP (with an e-value cutoff of 10 −10 ) in order to identify putativ e par alogs. Intr a genic homeologous bloc ks wer e defined as regions of ten or more genes with colinear or nearly colinear runs of paralogs elsewhere in the genome with fewer than six intervening genes . T hese intragenic homeologous blocks were visualized using CIRCOS version 0.69.8 [ 48 ]. Similarly, we also performed pairwise comparisons of input protein sequences from V . reflexo-pilosa , V . mungo , and V . radiata . Clustering was carried out using OrthoMCL software version 2.0.9 [ 49 ] based on a Markov clustering algorithm. Syntenic blocks between V. reflexo-pilosa , V. mungo , and V. radiata were identified by MCscanX and plotted with CIRCOS using the criteria mentioned above (at least ten colinear genes and fewer than six intervening genes allo w ed).
