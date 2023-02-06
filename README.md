# tx_star-bPAC_brain; Time-series zebrafish whole brain mRNA-seq analyses.

This project for the data generation and analyses of the study by Choi et al., 2023 titled: "Developmental Glucocorticoid exposure alters adult behavior and the brain-wide transcriptional landscape ". 

BAckgrounds of sequencing.
paired-end TruSeq Stranded mRNA libraries (Illumina, CA, USA) were constructed and over 20M of 50 bp reads/sample were sequenced. For this study, a total of 60 samples were sequenced, consisting of 5 biological replicates at four different time points (6, 13, 120 dpf, and acute-stressed at 120dpf) in every three groups (wildtype, bPAC+ and bPAC-). Raw sequenced reads(.fastq) from Illumina NovaSeq 6000 are avaiable to acceces through ENA (PRJEB53713). 

Read Mapping and quantification.
shell scripts(.sh) in Supplementary folder were used for mapping and quantification of reads in HPC environment. 
Briefly, after the quality control procedures of fastq files using FASTQC (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/), sequenced reads of each sample were processed by fastp{Chen, 2018} for removal of adapter sequences, low quality reads, and polyG in 3’ ends. Cleaned reads were mapped to the zebrafish reference genome assembly (GRCz11) with Ensemble annotation version 107 using HISAT2{Kim, 2015}. Samtools{Li, 2009} processed sam to bam file conversion. And then, quantification of expressed transcripts was achieved using Stringtie{Pertea, 2016}. Outputfiles from Stringtie were transformed for DESeq2 and edgeR by using prepDE.py (https://ccb.jhu.edu/software/stringtie/dl/prepDE.py3). "gtf.file.list" file contains sample IDs and paths (eg. https://ccb.jhu.edu/software/stringtie/dl/sample_lst.txt) 

> python prepDE.py -i gtf.file.list -l 50

DEG identification and downstream analysis.
Downstream and statistical analyzes were performed with R scripts, and files; "gene_count_matrix.csv", and "sample_info.csv".

Simply, sourced the R script file, "mRNA-seq_analysis.R" in the code folder.

R package edgeR{Robinson, 2010} was used to quantify the number of reads. DEGs were determined by which showed more than absolute 2-fold changes (FC) in their expression with under 0.01 of the False discovery rates (FDR). We used a lower fold change threshold (absolute 1.5 FC) for 120 dpf samples because there is no elevation of cortisol in bPAC+ compared to wildtype at 120 dpf. 

For the functional analyzes, R packages including clusterprofiler{Wu, 2021}, gProfiler2{Raudvere, 2019}, enrichGO{Wu, 2021} and annotation platforms (Chea3{Keenan, 2019}, DisGeNET{Piñero, 2021}, Enrichr {Chen, 2013}, Tabula_Muris{Schaum, 2018}) were used. String database{Szklarczyk, 2019} and Cytoscape{Shannon, 2003} were utilized for the protein-protein interaction network. The list of gene including social genes (101 genes) and epigenetic modifiers which are associated with DNA modification (52 genes) and histone modification (365 genes) as well as RNA modification and processing (1886 genes) was obtained from the Gene Ontology Annotation (GOA) database and they are available in the Supplementary folder. 


References

1	Chen, S., Zhou, Y., Chen, Y. & Gu, J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34, i884-i890 (2018). https://doi.org:10.1093/bioinformatics/bty560
2	Kim, D., Langmead, B. & Salzberg, S. L. HISAT: a fast spliced aligner with low memory requirements. Nat Methods 12, 357-360 (2015). https://doi.org:10.1038/nmeth.3317
3	Li, H. et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics 25, 2078-2079 (2009). https://doi.org:10.1093/bioinformatics/btp352
4	Pertea, M., Kim, D., Pertea, G. M., Leek, J. T. & Salzberg, S. L. Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown. Nature Protocols 11, 1650-1667 (2016). https://doi.org:10.1038/nprot.2016.095
5	Team, R. RStudio: Integrated Development for R.  (2020). 
6	Robinson, M. D., McCarthy, D. J. & Smyth, G. K. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139-140 (2010). https://doi.org:10.1093/bioinformatics/btp616
7	Wu, T. et al. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation 2, 100141 (2021). https://doi.org:https://doi.org/10.1016/j.xinn.2021.100141
8	Raudvere, U. et al. g:Profiler: a web server for functional enrichment analysis and conversions of gene lists (2019 update). Nucleic Acids Research 47, W191-W198 (2019). https://doi.org:10.1093/nar/gkz369
9	Keenan, A. B. et al. ChEA3: transcription factor enrichment analysis by orthogonal omics integration. Nucleic Acids Research 47, W212-W224 (2019). https://doi.org:10.1093/nar/gkz446
10	Piñero, J., Saüch, J., Sanz, F. & Furlong, L. I. The DisGeNET cytoscape app: Exploring and visualizing disease genomics data. Computational and Structural Biotechnology Journal 19, 2960-2967 (2021). https://doi.org:https://doi.org/10.1016/j.csbj.2021.05.015
11	Szklarczyk, D. et al. STRING v11: protein–protein association networks with increased coverage, supporting functional discovery in genome-wide experimental datasets. Nucleic Acids Research 47, D607-D613 (2019). https://doi.org:10.1093/nar/gky1131
12	Shannon, P. et al. Cytoscape: a software environment for integrated models of biomolecular interaction networks. Genome Res 13, 2498-2504 (2003). https://doi.org:10.1101/gr.1239303
13  Chen, E. Y. et al. Enrichr: interactive and collaborative HTML5 gene list enrichment analysis tool. BMC Bioinformatics 14, 128 (2013). https://doi.org:10.1186/1471-2105-14-128
14	Schaum, N. et al. Single-cell transcriptomics of 20 mouse organs creates a Tabula Muris. Nature 562, 367-372 (2018). https://doi.org:10.1038/s41586-018-0590-4

