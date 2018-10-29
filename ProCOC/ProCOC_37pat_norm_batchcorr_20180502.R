## ProCOC analysis cohort 1 37pat 56 samples plus 2 controls search from 25.04.2018 ##

## Step 1: Prepare R session ##

library(SWATH2stats)

library(data.table)
data.table 1.9.6  For help type ?data.table or https://github.com/Rdatatable/data.table/wiki
The fastest way to learn (by data.table authors): https://www.datacamp.com/courses/data-analysis-the-data-table-way
library(preprocessCore)
library(pheatmap)
library(RColorBrewer)

setwd('Documents/Biomarker_Cancer/ProCOC/ProCOC_final/ProCOC_preprocessing_SWATH2stats/ProCOC1_37pat_2ctrl/')

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
[1] 32088    63

colnames(data)
 [1] "ProteinName"        "PeptideSequence"    "FragmentIon"       
 [4] "high_100.1"         "high_100.2"         "high_201.1"        
 [7] "high_201.2"         "high_21.1"          "high_218.1"        
[10] "high_283.2"         "high_285.1"         "high_285.2"        
[13] "high_289.1"         "high_291.1"         "high_291.2"        
[16] "high_299.1"         "high_299.2"         "high_339.1"        
[19] "high_339.2"         "high_55.2"          "intermediate_22.1" 
[22] "intermediate_25.1"  "intermediate_26.1"  "intermediate_296.1"
[25] "intermediate_31.1"  "intermediate_31.2"  "intermediate_319.1"
[28] "intermediate_320.1" "intermediate_49.1"  "intermediate_57.1" 
[31] "intermediate_57.2"  "intermediate_79.1"  "intermediate_83.1" 
[34] "intermediate_83.2"  "intermediate_87.1"  "low_12.2"          
[37] "low_125.2"          "low_126.1"          "low_126.2"         
[40] "low_137.1"          "low_137.2"          "low_148.1"         
[43] "low_148.2"          "low_149.2"          "low_163.1"         
[46] "low_163.2"          "low_222.1"          "low_222.2"         
[49] "low_276.1"          "low_276.2"          "low_317.1"         
[52] "low_317.2"          "low_355.1"          "low_355.2"         
[55] "low_46.1"           "low_46.2"           "low_53.1"          
[58] "low_53.2"           "low_92.1"           "low_92.2"          
[61] "NA_.R1"             "NA_.R2"             "RT"       

#define batches and replicates for visualization
batch <- c(1,2,1,2,1,1,2,1,2,1,1,2,1,2,1,2,2,1,1,1,1,1,2,1,1,1,1,2,1,1,2,1,2,2,1,2,1,2,1,2,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,1)
length(batch)
[1] 59
replicate <- c(0,0,0,0,1,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0 ,0,0,0, 0,0,0, 0,0,0, 0,0,0,  0,0,0, 0,0,0,  0,2,0, 0,0,0, 0,0,0, 1,2)
length(replicate)
[1] 59

#log transform
ProCOC_log <- log2(as.matrix(data[,4:62]))
dim(ProCOC_log)
[1] 32088    59
ProCOC_log[which(is.infinite(ProCOC_log))]<-NA


anno_col <- data.frame(batch, replicate)
rownames(anno_col)<- colnames(ProCOC_log)

#visualize samples before normalization and batch correction
pheatmap(ProCOC_log, show_rownames=F, show_colnames=F, filename='heatmap_data_before_norm_batchcorr_20180502.pdf',
border_color=NA, cluster_rows=F, cluster_cols=T, annotation_col=anno_col)

#normali
ProCOC_norm <- normalize.quantiles(ProCOC_log)
colnames(ProCOC_norm)<- colnames(ProCOC_log)


batch1 <- ProCOC_norm[,which(batch==1)]
dim(batch1)
[1] 32088    34
batch2 <- ProCOC_norm[,which(batch==2)]
dim(batch2)
[1] 32088    25
batch1_corr <- batch1 - rowMeans(batch1, na.rm=T)
batch2_corr <- batch2 - rowMeans(batch2, na.rm=T)
ProCOC_batchcorr <- cbind(batch1_corr, batch2_corr)
dim(ProCOC_batchcorr)
[1] 32088    59

rm <- rowMeans(ProCOC_norm, na.rm=T)
ProCOC_batchcorr_final <-  ProCOC_batchcorr + rm

pheatmap(ProCOC_batchcorr_final, show_rownames=F, show_colnames=F, filename='heatmap_after_norm_batchcorr_20180502.pdf',
 border_color=NA, cluster_rows=F, cluster_cols=T, annotation_col=anno_col)

colnames(ProCOC_batchcorr_final)
 [1] "high_100.1"         "high_201.1"         "high_21.1"         
 [4] "high_218.1"         "high_285.1"         "high_289.1"        
 [7] "high_291.1"         "high_299.1"         "high_339.1"        
[10] "intermediate_22.1"  "intermediate_25.1"  "intermediate_26.1" 
[13] "intermediate_296.1" "intermediate_31.1"  "intermediate_319.1"
[16] "intermediate_320.1" "intermediate_49.1"  "intermediate_57.1" 
[19] "intermediate_79.1"  "intermediate_83.1"  "intermediate_87.1" 
[22] "low_126.1"          "low_137.1"          "low_148.1"         
[25] "low_163.1"          "low_222.1"          "low_276.1"         
[28] "low_317.1"          "low_355.1"          "low_46.1"          
[31] "low_53.1"           "low_92.1"           "NA_.R1"            
[34] "NA_.R2"             "high_100.2"         "high_201.2"        
[37] "high_283.2"         "high_285.2"         "high_291.2"        
[40] "high_299.2"         "high_339.2"         "high_55.2"         
[43] "intermediate_31.2"  "intermediate_57.2"  "intermediate_83.2" 
[46] "low_12.2"           "low_125.2"          "low_126.2"         
[49] "low_137.2"          "low_148.2"          "low_149.2"         
[52] "low_163.2"          "low_222.2"          "low_276.2"         
[55] "low_317.2"          "low_355.2"          "low_46.2"          
[58] "low_53.2"           "low_92.2" 

 
ProCOC_batchcorr_final <- ProCOC_batchcorr_final[, sort(colnames(ProCOC_batchcorr_final))]
colnames(ProCOC_batchcorr_final)
 [1] "high_100.1"         "high_100.2"         "high_201.1"        
 [4] "high_201.2"         "high_21.1"          "high_218.1"        
 [7] "high_283.2"         "high_285.1"         "high_285.2"        
[10] "high_289.1"         "high_291.1"         "high_291.2"        
[13] "high_299.1"         "high_299.2"         "high_339.1"        
[16] "high_339.2"         "high_55.2"          "intermediate_22.1" 
[19] "intermediate_25.1"  "intermediate_26.1"  "intermediate_296.1"
[22] "intermediate_31.1"  "intermediate_31.2"  "intermediate_319.1"
[25] "intermediate_320.1" "intermediate_49.1"  "intermediate_57.1" 
[28] "intermediate_57.2"  "intermediate_79.1"  "intermediate_83.1" 
[31] "intermediate_83.2"  "intermediate_87.1"  "low_12.2"          
[34] "low_125.2"          "low_126.1"          "low_126.2"         
[37] "low_137.1"          "low_137.2"          "low_148.1"         
[40] "low_148.2"          "low_149.2"          "low_163.1"         
[43] "low_163.2"          "low_222.1"          "low_222.2"         
[46] "low_276.1"          "low_276.2"          "low_317.1"         
[49] "low_317.2"          "low_355.1"          "low_355.2"         
[52] "low_46.1"           "low_46.2"           "low_53.1"          
[55] "low_53.2"           "low_92.1"           "low_92.2"          
[58] "NA_.R1"             "NA_.R2"  

ProCOC_mapDIA_input <- cbind(data[,1:3], ProCOC_batchcorr_final, data[,63])
dim(ProCOC_mapDIA_input)
[1] 32088    63

colnames(ProCOC_mapDIA_input)[63] <- 'RT'
colnames(ProCOC_mapDIA_input)
 [1] "ProteinName"        "PeptideSequence"    "FragmentIon"       
 [4] "high_100.1"         "high_100.2"         "high_201.1"        
 [7] "high_201.2"         "high_21.1"          "high_218.1"        
[10] "high_283.2"         "high_285.1"         "high_285.2"        
[13] "high_289.1"         "high_291.1"         "high_291.2"        
[16] "high_299.1"         "high_299.2"         "high_339.1"        
[19] "high_339.2"         "high_55.2"          "intermediate_22.1" 
[22] "intermediate_25.1"  "intermediate_26.1"  "intermediate_296.1"
[25] "intermediate_31.1"  "intermediate_31.2"  "intermediate_319.1"
[28] "intermediate_320.1" "intermediate_49.1"  "intermediate_57.1" 
[31] "intermediate_57.2"  "intermediate_79.1"  "intermediate_83.1" 
[34] "intermediate_83.2"  "intermediate_87.1"  "low_12.2"          
[37] "low_125.2"          "low_126.1"          "low_126.2"         
[40] "low_137.1"          "low_137.2"          "low_148.1"         
[43] "low_148.2"          "low_149.2"          "low_163.1"         
[46] "low_163.2"          "low_222.1"          "low_222.2"         
[49] "low_276.1"          "low_276.2"          "low_317.1"         
[52] "low_317.2"          "low_355.1"          "low_355.2"         
[55] "low_46.1"           "low_46.2"           "low_53.1"          
[58] "low_53.2"           "low_92.1"           "low_92.2"          
[61] "RT"                
write.table(ProCOC_mapDIA_input, file='ProCOC_37pat_mapDIA_input_norm_batchcorr_20180502.tsv', quote=F, row.names=F, sep='\t')

FETUB <- ProCOC_mapDIA_input[which(ProCOC_mapDIA_input$ProteinName=='sp|Q58D62|FETUB_BOVIN'),]
FETUB_pep <- aggregate(FETUB[,4:62], by=list(FETUB$ProteinName, FETUB$PeptideSequence), FUN=mean)

dim(FETUB)
[1] 54 63

A1AG <- ProCOC_mapDIA_input[which(ProCOC_mapDIA_input$ProteinName=='sp|Q3SZR3|A1AG_BOVIN'),]
dim(A1AG)
[1] 186  63
A1AG_pep <- aggregate(A1AG[,4:62], by=list(A1AG$ProteinName, A1AG$PeptideSequence), FUN=mean)

pdf('boxplot_FETUB_180719.pdf', width=10, height=10)
par(cex.axis=2.5, las=1, mgp = c(0, 1.5, 0))
boxplot(as.matrix(FETUB[,4:62]), ylim=c(0,20), col='lightblue', xaxt='n')
axis(1, at=c(5,10,15,20,25,30,35,40,45,50,55,60))
dev.off()


pdf('boxplot_A1AG_180719.pdf', width=10, height=10)
par(cex.axis=2.5, las=1, mgp = c(0, 1.5, 0))
boxplot(as.matrix(A1AG[,4:62]), ylim=c(0,20), col='lightblue', xaxt='n')
axis(1, at=c(5,10,15,20,25,30,35,40,45,50,55,60))
dev.off()

pdf('profile_FETUB_180719.pdf', width=15, height=10)
par(cex.axis=2.5, las=1, mgp = c(0, 1.5, 0))
matplot(t(as.matrix(FETUB_pep[, 3:61])), type='l', ylim=c(0,20), xaxt='n', ylab='', lwd=3)
axis(1, at=c(5,10,15,20,25,30,35,40,45,50,55,60))
dev.off()
 
pdf('profile_A1AG_180719.pdf', width=15, height=10)
par(cex.axis=2.5, las=1, mgp = c(0, 1.5, 0))
matplot(t(as.matrix(A1AG_pep[, 3:61])), type='l', ylim=c(0,20), xaxt='n', ylab='', lwd=3)
axis(1, at=c(5,10,15,20,25,30,35,40,45,50,55,60))
dev.off()


