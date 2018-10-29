## ProCOC analysis cohort 1 37pat 57 samples plus 2 controls search from 25.04.2018 ##

## Step 1: Prepare R session

#load required packages
library(SWATH2stats)
library(data.table)

# set working directory
setwd('Documents/Biomarker_Cancer/ProCOC/ProCOC_final/ProCOC_preprocessing_SWATH2stats/ProCOC1_37pat_2ctrl/')

#document information about versions for future reproducibility
sessionInfo()
R version 3.2.2 (2015-08-14)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 15.10

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=de_CH.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=de_CH.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=de_CH.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=de_CH.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] data.table_1.9.6  SWATH2stats_1.0.3

loaded via a namespace (and not attached):
[1] magrittr_1.5   plyr_1.8.3     tools_3.2.2    reshape2_1.4.1 Rcpp_0.12.7   
[6] stringi_1.0-1  grid_3.2.2     stringr_1.0.0  chron_2.3-47  
## Step 2 - Read in data ##

#read in data
data <- fread('/home/betty/Documents/Biomarker_Cancer/ProCOC/ProCOC_final/ProCOC_Data/ProCOC1_37pat_2ctrl_data/E1804251107_feature_alignment_requant.tsv', sep='\t', header=T)
Read 506515 rows and 58 (of 58) columns from 0.567 GB file in 00:00:07

#convert to dataframe
data <- as.data.frame(data)
dim(data)
[1] 506515     58


#read in annotation file
annotation.file <- read.delim2("", dec=".", sep=",", header=T)
dim(annotation.file)
[1] 59  4


head(annotation.file)
                      Filename Condition BioReplicate Run
1  sajict_L140907_SW_100-1      high        100-1   1
2 fbetty_L151027_SW_100-1C      high        100-2   2
3  sajict_L140907_SW_126-1       low        126-1   3
4  fbetty_L151027_SW_126-1       low        126-2   4
5  sajict_L140907_SW_137-1       low        137-1   5
6  fbetty_L151027_SW_137-1       low        137-2   6
                     

## Step 3 - Assess m-score/FDR ##

#assess decoy rate
assess_decoy_rate(data)
Number of non-decoy peptides: 3442
Number of decoy peptides: 251
Decoy rate: 0.0729


#assess a global decoy rate and estimate an FDR
assess_fdr_overall(data, FFT=0.4, filename="ProCOC_37pat_FDR_overall")
ProCOC_37pat_FDR_overall.pdf written to working folder

ProCOC_37pat_FDR_overalltable.csv written to working folder

#assess FDR for individual runs
assess_fdr_byrun(data, FFT=0.4, filename="ProCOC_37pat_FDR_run")
The average FDR by run on assay level is 0.007

The average FDR by run on peptide level is 0.014

The average FDR by run on protein level is 0.101

Individual run FDR qualities can be retrieved from ProCOC_37pat_FDR_run .csv

ProCOC_37pat_FDR_run.csv reports written to working folder

ProCOC_37pat_FDR_run.pdf report plots written to working folder


#find suitable m-score to obtain specific global assay FDR
mscore4assayfdr(data, FFT=0.4, fdr_target=0.01)
Target assay FDR: 0.01

Required overall m-score cutoff:0.01
achieving assay FDR =0.00625

[1] 0.01




#find suitable m-score to obtain specific global peptide FDR
mscore4pepfdr(data, FFT=0.4, fdr_target=0.01)
Target peptide FDR:0.01

Required overall m-score cutoff:3.5481e-05
achieving peptide FDR =0.00987

[1] 3.548134e-05

mscore4pepfdr(data, FFT=0.4, fdr_target=0.02)
Target peptide FDR:0.02

Required overall m-score cutoff:0.00015849
achieving peptide FDR =0.0195

[1] 0.0001584893

#find suitable m-score to obtain specific global protein FDR
mscore4protfdr(data, FFT=0.4, fdr_target=0.01)
Target protein FDR:0.01

Required overall m-score cutoff:2.8184e-08
achieving protein FDR =0.0093

[1] 2.818383e-08

mscore4protfdr(data, FFT=0.4, fdr_target=0.02)
Target protein FDR:0.02

Required overall m-score cutoff:8.9125e-08
achieving protein FDR =0.0197

[1] 8.912509e-08


## Step 4 - Reduce data ##

#remove unnecessary columns
data.reduced <- reduce_OpenSWATH_output(data)
dim(data.reduced)
[1] 506515     12


#remove unneeded proteins
data.reduced <- data.reduced[grep("iRT", data.reduced$ProteinName, invert=TRUE),]
dim(data.reduced)
[1] 506515     12


data.reduced <- data.reduced[grep("Glyco_LightPeptideMixture_171",data.reduced$ProteinName, invert=TRUE),]
dim(data.reduced)
[1] 501087     12

data.reduced <- data.reduced[grep("Glyco_HeavyPeptideMixture_171", data.reduced$ProteinName, invert=TRUE),]
dim(data.reduced)
[1] 496013     12


## Step 5 - Annotate data ##
nlevels(as.factor(data.reduced$align_origfilename))
[1] 57

#annotate the data
data.annotated <- sample_annotation(data.reduced, annotation.file)
dim(data.annotated)
[1] 496013     15


head(unique(data.annotated$ProteinName))
[1] "1/Q14624" "1/P01876" "1/P02790" "1/P05155" "1/P00751" "1/P08697"


#OPTIONAL, reduce non-unique information in protein ProteinName
#data.annotated$ProteinName <- gsub("sp\\|([[:alnum:]]+)\\|[[:alnum:]]*_HUMAN", "\\1", data.annotated$ProteinName)
#head(unique(data.annotated$ProteinName))
#[1] "1/Q14624" "1/P01876" "1/P02790" "1/P05155" "1/P00751" "1/P08697"

## Step 6 - Filter data based on obtained m-score ##

#number of re-quant values
length(which(data.annotated$m_score==2))
[1] 188019

#percentage of requant values is 37.9%
 188019/496013
[1] 0.3790606


#data is filtered by mscore and percentage of replicates which need to reach this mscore filter_mscore_requant(data, mscore, percentage, rm.decoy=T)
#
data.filtered.mscore01.requant60 <- filter_mscore_requant(data.annotated, 0.01, 0.6, rm.decoy=T)
Treshold, peptides need to have been quantified in more conditions than: 35.4
Fraction of peptides selected: 0.7
Dimension difference: 135700, 0

dim(data.filtered.mscore01.requant60)
[1] 345563     16

length(which(data.filtered.mscore01.requant60$m_score==2))
[1] 64945

# percentage of requant values reduced to 18.8%
64945/345563
[1] 0.1879397

#filter for proteotypic peptides
data.filtered.mscore01.requant60.proteotypic <- filter_proteotypic_peptides(data.filtered.mscore01.requant60)
Number of proteins detected: 232
Protein identifiers: P02765, P01008, Q08380, P13598, O75882, P08185
Number of proteins detected that are supported by a proteotypic peptide: 197
Number of proteotypic peptides detected: 2134

dim(data.filtered.mscore01.requant60.proteotypic)
[1] 315532     16

#requant values and percentage 19%
length(which(data.filtered.mscore01.requant60.proteotypic$m_score==2))
[1] 59837
59837/315532
[1] 0.1896385

## Step 7 - Output on various levels and Conversion various tools ##

#output on peptide levels
data.peptide <- data.filtered.mscore01.requant60.proteotypic
dim(data.peptide)
[1] 303297     16


data.peptide$aggr_Fragment_Annotation <- NULL
data.peptide$aggr_Peak_Area <- NULL
dim(data.peptide)
[1] 303297     14

write.csv(data.peptide, file="ProCOC_37pat_mscore01r60_pep_20180502.csv", row.names=F, quote=F)

#alternative peptide level converted in wide format, meaning samples to columns
#write_matrix_peptides(data.filtered.proteotypic, filename="ProCOC_151020_matrix_pep_level", rm.decoy=F)

#splitting transitions needed for MSstats, mapDIA, aLFQ
data.transition <- disaggregate(data.filtered.mscore01.requant60.proteotypic)
dim(data.transition)
[1] 1893192      10


#list all objects
ls()
#save intermediate steps
save(data.annotated, data.filtered.mscore01.requant60, data.filtered.mscore01.requant60.proteotypic, data.peptide, data.reduced, data.transition, file='ProCOC_37pat_preprocessing_workflow_20180502.Rdata')


#conversion for MSstats
MSstats.input <- convert4MSstats(data.transition)
One or several columns required by MSstats were not in the data and filled with NAs (except IsotopeLabelType is filled with light).
              Missing columns: ProductChargeIsotopeLabelType
Warning message:
In convert4MSstats(data.transition) :
  Intensity values that were 0, were replaced by NA


head(MSstats.input) 
  ProteinName               PeptideSequence PrecursorCharge
1      P02765 AAFNAQNN(UniMod_7)GSNFQLEEISR               2
2      P02765 AAFNAQNN(UniMod_7)GSNFQLEEISR               2
3      P02765 AAFNAQNN(UniMod_7)GSNFQLEEISR               2
4      P02765 AAFNAQNN(UniMod_7)GSNFQLEEISR               2
5      P02765 AAFNAQNN(UniMod_7)GSNFQLEEISR               2
6      P02765 AAFNAQNN(UniMod_7)GSNFQLEEISR               2
                          FragmentIon ProductCharge IsotopeLabelType Intensity
1 263_AAFNAQNN(UniMod_7)GSNFQLEEISR_2            NA            light       140
2 250_AAFNAQNN(UniMod_7)GSNFQLEEISR_2            NA            light       190
3 263_AAFNAQNN(UniMod_7)GSNFQLEEISR_2            NA            light       190
4 250_AAFNAQNN(UniMod_7)GSNFQLEEISR_2            NA            light      5417
5 263_AAFNAQNN(UniMod_7)GSNFQLEEISR_2            NA            light        20
6 263_AAFNAQNN(UniMod_7)GSNFQLEEISR_2            NA            light       262
  BioReplicate    Condition Run
1        100-2         high   2
2        137-1          low   5
3        137-1          low   5
4         57-2 intermediate  46
5         83-2 intermediate  49
6         83-1 intermediate  48

dim(MSstats.input)
[1] 1819782      10

save(MSstats.input, file="ProCOC_37pat_20180502_MSstats_input.Rdata")

#conversion for mapDIA
mapDIA.input <- convert4mapDIA(data.transition, RT=T)

dim(mapDIA.input)
[1] 32088    63


write.table(mapDIA.input, file="ProCOC_37pat_20180502_mapDIA_input.txt", quote=F, row.names=F, sep="\t")
