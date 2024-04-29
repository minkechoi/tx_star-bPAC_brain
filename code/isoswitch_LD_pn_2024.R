suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(IsoformSwitchAnalyzeR))
suppressPackageStartupMessages(library(ggplot2))

options(readr.show_progress=FALSE)

#### A newly created switchAnalyzeRlist + switch analysis

# part 1 ------------------------------------------------------------------


### Import data into R via one of:
#myQantifications <- importIsoformExpression() # RSEM/Kallisto/Salmon/StringTie

### Please note
# The way of importing files in the following example with
# "system.file('pathToFile', package="IsoformSwitchAnalyzeR") is
# specialized way of accessing the example data in the IsoformSwitchAnalyzeR package
# and not something  you need to do - just supply the string e.g.
# parentDir = "/myLD.stringtieQuantifications/" pointing to the parent directory (where 
# each sample is a seperate sub-directory) to the function.

### Import stringtie example data in R package
library(BSgenome.Drerio.UCSC.danRer11)
genome <- BSgenome.Drerio.UCSC.danRer11
# change genome to ensembl style
seqnames(genome) = gsub("chr", "", seqnames(genome))
seqnames(genome)[26] = "MT"

#load stringtie files
LD.stringtieQuant <- importIsoformExpression(
  parentDir = file.path("string_2nd","LD.string"),
  pattern = "bg.0256_",addIsofomIdAsColumn = FALSE, readLength = 50
)


#> Step 1 of 3: Identifying which algorithm was used...
#>     The quantification algorithm used was: Salmon
#>     Found 6 quantification file(s) of interest
#> Step 2 of 3: Reading data...
#> reading in files with read_tsv
#> 1 2 3 4 5 6 
#> Step 3 of 3: Normalizing abundance values (not counts) via edgeR...
#> Done

### Make design matrix
sample_info = read.table("./data/string.gtf.list",sep = "\t")
condi = c()
for (i in 1:15) {
  samp = str_split(sample_info$V1[grep(colnames(LD.stringtieQuant$abundance)[i],sample_info$V2)],"_")[[1]][2]
  condi = c(condi,samp)
}


myDesign <- data.frame(
  sampleID = colnames(LD.stringtieQuant$abundance),
  condition = condi
)

#for pos vs neg


pos.neg.myDesign <-myDesign[c(6:15),]

### Please note
# The way of importing files in the following example with
# "system.file('pathToFile', package="IsoformSwitchAnalyzeR") is
# specialized way of accessing the example data in the IsoformSwitchAnalyzeR package
# and not something  you need to do - just supply the string e.g.:
# isoformExonAnnoation = "/myAnnotations/annotation.gtf".

#merged gtf file from HISAT2 running

### Create switchAnalyzeRlist
switchAnalyzeRlist.pn <- importRdata(
  isoformCountMatrix   = LD.stringtieQuant$counts,
  isoformRepExpression = LD.stringtieQuant$abundance,
  designMatrix         = pos.neg.myDesign,
  isoformExonAnnoation = "./tx_star_merged.gtf",
  #addAnnotatedORFs = 
  #isoformNtFasta       = 
  removeNonConvensionalChr = TRUE,
  fixStringTieAnnotationProblem = TRUE,
  showProgress = FALSE
)

### Extract gene count matrix
geneCountMatrix.pn <- extractGeneExpression(
  switchAnalyzeRlist.pn,
  extractCounts = TRUE # set to FALSE for abundances
)


#> Step 1 of 7: Checking data...
#> Step 2 of 7: Obtaining annotation...
#>     importing GTF (this may take a while)...
#> Step 3 of 7: Fixing StringTie gene annoation problems...
#>     There were no need to rescue any annotation
#>     302 genes_id were assigned their original gene_id instead of the StringTie gene_id.
#>         This was only done when it could be done unambiguous.
#> Step 4 of 7: Calculating gene expression and isoform fractions...
#> Step 5 of 7: Merging gene and isoform expression...
#> Step 6 of 7: Making comparisons...
#> Step 7 of 7: Making switchAnalyzeRlist object...
#> The GUESSTIMATED number of genes with differential isoform usage are:
#>            comparison estimated_genes_with_dtu
#> 1 Fibroblasts vs hESC                  18 - 30
#> 2  Fibroblasts vs iPS                  49 - 82
#> 3         hESC vs iPS                  27 - 46
#> Done

aSwitchList.pn=switchAnalyzeRlist.pn

aSwitchList.pn=analyzeORF(switchAnalyzeRlist = aSwitchList.pn,genomeObject = genome)

### Filter
ft.aSwitchList.pn <- preFilter(aSwitchList.pn
                               # geneExpressionCutoff = 10
)

### Test for isoform switches
ft.SwitchList.pn <- isoformSwitchTestDEXSeq( ft.aSwitchList.pn,
                                             reduceToSwitchingGenes = TRUE )  

### If analysing (some) novel isoforms (else use CDS from ORF as explained in importRdata() )
#ft.SwitchList.pn <- addORFfromGTF( ft.SwitchList.pn, pathToGTF ="./tx_star_merged.gtf")
#ft.SwitchList.pn <- analyzeNovelIsoformORF( ft.SwitchList.pn,genomeObject = genome, analysisAllIsoformsWithoutORF = F)

### Extract Sequences

#####

ft.SwitchList.pn <- extractSequence(
  ft.SwitchList.pn, 
  pathToOutput = './',
  outputPrefix = "pn.LD",
  writeToFile=T # to avoid output when running this example data
)

### Summary
extractSwitchSummary(ft.SwitchList.pn)
#save
#save.image(file = "isoswitch_part1LD.RData")

#### out_data

#CPC2  => run;  python3 ./CPC2_standalone-1.0.1/bin/CPC2.py -i ../pn.LD_nt.fasta -o ./pn.LD.cpc2.txt
#pFAM  => run; cpan Moose => perl pfam/PfamScan/pfam_scan.pl -fasta ./pn.LD_AA.fasta -dir pfam/pfamdb/ > LD.pn.pfam.txt
#SignalP => run in online,https://services.healthtech.dtu.dk/services/SignalP-5.0/;  short output; download summary
#IUPred2A  =>  run in online, https://iupred2a.elte.hu/; download report from email.
#NetSurfP3 => local run; python nsp3.py -m models/nsp3.pth -i ../../../pn.LD_AA.fasta -o ../../nsp3.pn.LD.out
#DeepLoc2  => local run; deeploc2 -f pn.LD_AA.fasta -o ../../nsp3.pn.LD.out
#DeepTMHMM
'
analyzeCPAT() # OR
analyzeCPC2()
analyzePFAM()
analyzeSignalP()
analyzeIUPred2A() # OR
analyzeNetSurfP2()
analyzeDeepLoc2()
analyzeDeepTMHMM()
'

#analyzeCPC2()

ft.SwitchList.pn <- analyzeCPC2(
  switchAnalyzeRlist   = ft.SwitchList.pn,
  pathToCPC2resultFile = './data/pn.LD.cpc2.txt',
  removeNoncodinORFs   = TRUE   # because ORF was predicted de novo
)

### Add PFAM analysis
ft.SwitchList.pn <- analyzePFAM(
  switchAnalyzeRlist   = ft.SwitchList.pn,
  pathToPFAMresultFile = './data/LD.pn.pfam.txt',
  showProgress=FALSE
)
#> Converting AA coordinats to transcript and genomic coordinats...
#> Added domain information to 127 (78.4%) transcripts

### Add SignalP analysis
ft.SwitchList.pn <- analyzeSignalP(
  switchAnalyzeRlist       = ft.SwitchList.pn,
  pathToSignalPresultFile  = './data/signalIP_pn.LD.output_protein_type.txt'
)
#> Added signal peptide information to 17 (10.49%) transcripts

### Add IUPred2A analysis
ft.SwitchList.pn <- analyzeIUPred2A(
  switchAnalyzeRlist        = ft.SwitchList.pn,
  pathToIUPred2AresultFile = './data/IUPred2A_result.pn.LD.txt',
  showProgress = FALSE
)

### Add NetSurfP analysis

ft.SwitchList.pn <- analyzeNetSurfP2(
  switchAnalyzeRlist        = ft.SwitchList.pn,
  pathToNetSurfP2resultFile = "./data/netsurfp_pn.csv",
  showProgress = FALSE
)

'
### Add DeepLoc2 analysis
exampleSwitchListAnalyzed <- analyzeDeepLoc2(
  switchAnalyzeRlist = exampleSwitchListAnalyzed,
  pathToDeepLoc2resultFile = system.file("extdata/deeploc2.csv", package = "IsoformSwitchAnalyzeR"),
  quiet = FALSE
)
#> Added subcellular information to 135 (83.33%) transcripts

### Add DeepTMHMM analysis
exampleSwitchListAnalyzed <- analyzeDeepTMHMM(
  switchAnalyzeRlist   = exampleSwitchListAnalyzed,
  pathToDeepTMHMMresultFile = system.file("extdata/DeepTMHMM.gff3", package = "IsoformSwitchAnalyzeR"),
  showProgress=FALSE
)
#> Step 1 of 2: Reading results into R...
#> Step 2 of 2: Converting AA coordinats to transcript and genomic coordinats...
#> Added topology information to 141 transcripts (87.04%).
'

#alternative splicing


ft.SwitchList.pn <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = ft.SwitchList.pn,
  quiet=TRUE
)

alt_sp_tb.pn=ft.SwitchList.pn$AlternativeSplicingAnalysis



###
as_switch_tb.pn=ft.SwitchList.pn$isoformFeatures
tx.ID.pn = unique(as_switch_tb.pn$isoform_id)

###
library(biomaRt)
listMarts()
ensembl = useEnsembl(biomart="ensembl", mirror = "www") #if not working, try with "useast"
datasets <- listDatasets(ensembl)
head(datasets)
ensembl = useDataset("drerio_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

ft_txid.anno.pn=getBM(attributes =c("ensembl_gene_id","ensembl_transcript_id","external_transcript_name","transcript_biotype","description",
                                    "transcript_is_canonical","strand","chromosome_name","transcript_start","transcript_end","transcription_start_site","transcript_length"),
                      filters = data.frame("ensembl_transcript_id" = tx.ID.pn), mart = ensembl)

colnames(ft_txid.anno.pn)=c("ensembl_gene_id","isoform_id","external_transcript_name","transcript_biotype","description",
                            "transcript_is_canonical","strand","chromosome_name","transcript_start","transcript_end","transcription_start_site","transcript_length")

as_switch_tb.pn=dplyr::full_join(x=as_switch_tb.pn,y=ft_txid.anno.pn, copy=T, by="isoform_id")
write.csv(as_switch_tb.pn,"./output/as_switch_tb.pn.LD.csv")

ft_as_switch_tb.pn = as_switch_tb.pn %>% filter(gene_switch_q_value <0.05) %>% filter(abs(dIF) > 0.3)

ft_as_switch_tb.pn = ft_as_switch_tb.pn[c(which(ft_as_switch_tb.pn$gene_value_1 >=10|ft_as_switch_tb.pn$gene_value_2 >= 10)),]
#apr
ft_as_switch_tb.pn= as_switch_tb.pn %>% filter(gene_switch_q_value < 0.05) %>% filter(abs(dIF) > 0.4) 
ft_as_switch_tb.pn = ft_as_switch_tb.pn[c(which(ft_as_switch_tb.pn$iso_value_1 >=10|ft_as_switch_tb.pn$iso_value_2 >= 10)),]

unknown_ft = ft_as_switch_tb.pn$gene_name[c(grep("MSTRG",ft_as_switch_tb.pn$isoform_id))]

ft_as_switch_tb.pn =  ft_as_switch_tb.pn[-c(which(ft_as_switch_tb.pn$gene_name %in%unknown_ft )),]
d=duplicated(ft_as_switch_tb.pn$gene_name)
ft_as_switch_tb.pn= ft_as_switch_tb.pn[c(which(ft_as_switch_tb.pn$gene_name %in% ft_as_switch_tb.pn[d,"gene_name"])),]
write.csv(ft_as_switch_tb.pn,"./ft_as_switch_tb.pn.LD_apr.csv")

#save
save.image(file = "isoswitch_part1LD.pn.apr.RData")

#### analyzeSwitchesConsequences

consequencesOfInterest <- c('intron_retention',
                            'coding_potential',
                            'NMD_status',
                            'domains_identified',
                            'ORF_seq_similarity')



ft.SwitchList.pn <- analyzeSwitchConsequences(
  ft.SwitchList.pn,
  consequencesToAnalyze = consequencesOfInterest, 
  dIFcutoff = 0.4, # very high cutoff for fast runtimes - you should use the default (0.1)
  showProgress=FALSE
)

extractSwitchSummary(ft.SwitchList.pn, dIFcutoff = 0.4, filterForConsequences = FALSE)

sigLD_gene_as_switch <- subsetSwitchAnalyzeRlist(
  ft.SwitchList.pn, 
  ft.SwitchList.pn$isoformFeatures$gene_name %in% ft_as_switch_tb.pn$gene_name
)

####plot

#switch plot
sigLD_gene_as_switch
top.pn=extractTopSwitches(sigLD_gene_as_switch,dIFcutoff = 0.4,  n=100)
gene.list = na.omit(unique(top.pn$gene_name))
gene.list = top.pn$gene_name


pdf(file = "./figures/top_sig_switchplot_pn.LD_apr_rev.pdf",   # The directory you want to save the file in
    width = 11.69, # The width of the plot in inches
    height = 8.27 ) # The height of the plot in inches

par(mfrow=c(1,2))

for (i in 1:length(gene.list)) {
  
  switchPlot(ft.SwitchList.pn,gene =gene.list[i])
  print(i)
}

dev.off()
