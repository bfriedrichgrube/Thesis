## ProCOC analysis cohort 1 56pat plus 2 controls search from 25.04.2018 ##

## Step 1: Prepare R session ##

library(SWATH2stats)

library(data.table)
data.table 1.9.6  For help type ?data.table or https://github.com/Rdatatable/data.table/wiki
The fastest way to learn (by data.table authors): https://www.datacamp.com/courses/data-analysis-the-data-table-way
library(preprocessCore)
library(pheatmap)
library(RColorBrewer)

setwd('Documents/Biomarker_Cancer/ProCOC/ProCOC_final/ProCOC_preprocessing_SWATH2stats/ProCOC2_56pat_2ctrl/')

#sessionInfo
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
[1] RColorBrewer_1.1-2    pheatmap_1.0.8        preprocessCore_1.32.0
[4] data.table_1.9.6      SWATH2stats_1.0.3    

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.7      chron_2.3-47     grid_3.2.2       plyr_1.8.3      
 [5] gtable_0.2.0     magrittr_1.5     scales_0.4.1     stringi_1.0-1   
 [9] reshape2_1.4.1   tools_3.2.2      stringr_1.0.0    munsell_0.4.3   
[13] colorspace_1.2-6

# Step 2: Load in data ##
data <- read.delim('ProCOC_37pat_20180502_mapDIA_input.txt', sep='\t', header=T, as.is=T)
dim(data)
[1] 30276    62

colnames(data)
 [1] "ProteinName"     "PeptideSequence" "FragmentIon"     "high_10"        
 [5] "high_159"        "high_165"        "high_166"        "high_176"       
 [9] "high_187"        "high_193"        "high_198"        "high_200"       
[13] "high_214"        "high_224"        "high_235"        "high_237"       
[17] "high_242"        "high_256"        "high_279"        "high_294"       
[21] "high_305"        "high_315"        "high_32"         "high_340"       
[25] "high_352"        "high_37"         "high_42"         "high_43"        
[29] "high_47"         "high_58"         "high_69"         "high_72"        
[33] "high_73"         "high_81"         "high_93"         "low_11"         
[37] "low_120"         "low_121"         "low_131"         "low_138"        
[41] "low_15"          "low_156"         "low_164"         "low_179"        
[45] "low_182"         "low_196"         "low_20"          "low_212"        
[49] "low_23"          "low_28"          "low_286"         "low_300"        
[53] "low_311"         "low_321"         "low_322"         "low_327"        
[57] "low_34"          "low_361"         "low_39"          "NA_R1"          
[61] "NA_R2"           "RT"  


#define replicates for visualization
replicate <- c(1, 0,0,0,0, 0,0,0,0, 0,2,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,1, 2)
length(replicate)
[1] 58

#log transform
ProCOC_log <- log2(as.matrix(data[,4:61]))
dim(ProCOC_log)
[1] 30276    58

ProCOC_log[which(is.infinite(ProCOC_log))]<-NA


anno_col <- data.frame(replicate)
rownames(anno_col)<- colnames(ProCOC_log)

#visualize samples before normalization and batch correction
pheatmap(ProCOC_log, show_rownames=F, show_colnames=F, filename='heatmap_data_before_norm_20180503.pdf',
border_color=NA, cluster_rows=F, cluster_cols=T, annotation_col=anno_col)

#normali
ProCOC_norm <- normalize.quantiles(ProCOC_log)
colnames(ProCOC_norm)<- colnames(ProCOC_log)

pheatmap(ProCOC_norm, show_rownames=F, show_colnames=F, filename='heatmap_after_norm_batchcorr_20180502.pdf',
 border_color=NA, cluster_rows=F, cluster_cols=T, annotation_col=anno_col)


ProCOC_mapDIA_input <- cbind(data[,1:3], ProCOC_norm, data[,62])
dim(ProCOC_mapDIA_input)
[1] 30276    62

colnames(ProCOC_mapDIA_input)[62] <- 'RT'
colnames(ProCOC_mapDIA_input)
                
write.table(ProCOC_mapDIA_input, file='ProCOC_56pat_mapDIA_input_norm_20180503.tsv', quote=F, row.names=F, sep='\t')

FETUB <- ProCOC_mapDIA_input[which(ProCOC_mapDIA_input$ProteinName=='sp|Q58D62|FETUB_BOVIN'),]
dim(FETUB)
[1] 54 62
FETUB_pep <- aggregate(FETUB[,4:61], by=list(FETUB$ProteinName, FETUB$PeptideSequence), FUN=mean)

A1AG <- ProCOC_mapDIA_input[which(ProCOC_mapDIA_input$ProteinName=='sp|Q3SZR3|A1AG_BOVIN'),]
dim(A1AG)
[1] 174  62
A1AG_pep <- aggregate(A1AG[,4:61], by=list(A1AG$ProteinName, A1AG$PeptideSequence), FUN=mean)

pdf('boxplot2_FETUB_180719.pdf', width=10, height=10)
par(cex.axis=2.5, las=1, mgp = c(0, 1.5, 0))
boxplot(as.matrix(FETUB[,4:61]), ylim=c(0,20), col='lightblue', xaxt='n')
axis(1, at=c(5,10,15,20,25,30,35,40,45,50,55,60))
dev.off()


pdf('boxplot2_A1AG_180719.pdf', width=10, height=10)
par(cex.axis=2.5, las=1, mgp = c(0, 1.5, 0))
boxplot(as.matrix(A1AG[,4:61]), ylim=c(0,20), col='lightblue', xaxt='n')
axis(1, at=c(5,10,15,20,25,30,35,40,45,50,55,60))
dev.off()

pdf('profile2_FETUB_180719.pdf', width=15, height=10)
par(cex.axis=2.5, las=1, mgp = c(0, 1.5, 0))
matplot(t(as.matrix(FETUB_pep[, 3:60])), type='l', ylim=c(0,20), xaxt='n', ylab='', lwd=3)
axis(1, at=c(5,10,15,20,25,30,35,40,45,50,55,60))
dev.off()
 
pdf('profile2_A1AG_180719.pdf', width=15, height=10)
par(cex.axis=2.5, las=1, mgp = c(0, 1.5, 0))
matplot(t(as.matrix(A1AG_pep[, 3:60])), type='l', ylim=c(0,20), xaxt='n', ylab='', lwd=3)
axis(1, at=c(5,10,15,20,25,30,35,40,45,50,55,60))
dev.off()
pdf('boxplot_FETUB_20180503.pdf', width=10, height=10)
boxplot(as.matrix(FETUB[,4:61]), ylim=c(0,20), col='lightblue')
dev.off()

pdf('boxplot_A1AG_20180503.pdf', width=10, height=10)
boxplot(as.matrix(A1AG[,4:61]), ylim=c(0,20), col='lightblue')
dev.off()

	

