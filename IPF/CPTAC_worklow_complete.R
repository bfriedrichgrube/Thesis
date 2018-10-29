## Step 1: Prepare R session ##


library(preprocessCore)
setwd('Documents/Biomarker_Cancer/CPTAC_IPF_final/')


## Step 2: Load data ##

data <- read.delim('feature_alignment_matrix.csv', sep=',', header=T)
annotation.file <- read.delim2('CPTAC_IPF_annotation.txt', dec='.', sep='\t', header=T)
batch_annotation <- read.delim2('CPTAC_IPF_batch_annotation.csv', sep=',', header=T)

dim(data)
[1] 36925   103

## Step 3: Remove non-proteotypic peptides ##

data_proteotypic <- data[grep(";", data$Protein, invert=TRUE),]
dim(data_proteotypic)
[1] 34332   103

## Step 4: Annotate samples

col_n <- paste(annotation.file[,2], annotation.file[,3], sep='_')
head(col_n)
colnames(data_proteotypic)[3:100]<- col_n
head(colnames(data_proteotypic))

## Step 5: log2 transform and Normalization ##

data_log <- log2(as.matrix(data_proteotypic[,3:100]))
data_log[which(is.infinite(data_log))] <- NA # replace the infinite values with NA

dim(data_log)
[1] 34332    98

data_norm <- normalize.quantiles(data_log)
colnames(data_norm) <- colnames(data_log)
data_norm[1:3, 1:4]

## Step 6: Batch correction ##



data_batchA <- data_norm[,c(3:6, 10, 14, 17, 21, 25, 26, 31, 35, 38, 40, 41, 47, 66, 68, 69, 82, 83, 86, 87, 94)]
data_batchB <- data_norm[,c(7, 15, 33, 36, 39, 51, 53, 55, 56, 58, 60, 61:63, 65, 74, 77, 81, 84, 85, 88, 90, 92, 98)]
data_batchC <- data_norm[,c(1, 8, 12, 20, 23, 24, 27, 32, 34, 42, 44, 49, 54, 59, 64, 67, 72, 73, 75, 76, 78, 80, 89, 96)]
data_batchD <- data_norm[,c(13, 18, 19, 22, 28:30, 37, 43, 45, 48, 50, 52, 95, 97)]
data_batchE <- data_norm[,c(2, 9, 11, 16, 46, 57, 70, 71, 79, 91, 93)]

data_batchA_corr <- as.matrix(data_batchA) - rowMeans(as.matrix(data_batchA), na.rm=T)
data_batchB_corr <- as.matrix(data_batchB) - rowMeans(as.matrix(data_batchB), na.rm=T)
data_batchC_corr <- as.matrix(data_batchC) - rowMeans(as.matrix(data_batchC), na.rm=T)
data_batchD_corr <- as.matrix(data_batchD) - rowMeans(as.matrix(data_batchD), na.rm=T)
data_batchE_corr <- as.matrix(data_batchE) - rowMeans(as.matrix(data_batchE), na.rm=T)
data_batchcorr <- cbind(data_batchA_corr, data_batchB_corr, data_batchC_corr, data_batchD_corr, data_batchE_corr)

colnames(data_batchcorr)
data_batchcorr_order <- data_batchcorr[,order(colnames(data_batchcorr))]
data_norm_order <- data_norm[,order(colnames(data_norm))]

r_m <- rowMeans(data_norm_order,na.rm=T)
data_batchcorr_order_final <- data_batchcorr_order + r_m

save(data_batchcorr, data_batchcorr_order, data_norm, data_norm_order, data_log, data_proteotypic, file='CPTAC_IPF_workflow_process_20181011.Rdata')

data_final <- cbind(data_proteotypic[,1:2], data_batchcorr_order_final)
save(data_final, file='CPTAC_data_processed.Rdata')


## Step 7: Report numbers of modified peptides in matrix ##

length(unique(data_final$Protein))
[1] 5459

dim(data_final)
[1] 34332   100
# 34332 proteotypic peptides

length(grep('\\(', data_final$Peptide))
[1] 15773
# 15773 modified peptides

length(grep('\\(Carbamidomethyl\\)', data_final$Peptide))
[1] 6546

length(grep('\\(Deamidated\\)', data_final$Peptide))
[1] 9179

length(grep('\\(Oxidation\\)', data_final$Peptide))
[1] 2594

length(grep('K\\(Acetyl\\)', data_final$Peptide))
[1] 45

length(grep('\\.\\(Acetyl\\)', data_final$Peptide))
[1] 130

length(grep('\\(Phospho\\)', data_final$Peptide))
[1] 66

##############################################
# make venn for figures (library numbers, see separate script) paste into other script later when uploading

library(VennDiagram)

#make venn for protein overlap of both libraries


pdf('venn_library_proteins.pdf', width=10, height=10)
draw.pairwise.venn(9021, 8699, 7624, category=c('standard search engines - Chapter 5', 'MSFraggerGUI - Modifications'), lty=rep('blank', 2), fill=c('lightblue', 'green'), cat.pos=c(-25,100), cex=3, cat.cex=3,margin=0.05)
dev.off()

pdf('venn_library_proteins_nolabels.pdf', width=10, height=10)
draw.pairwise.venn(9021, 8699, 7624, lty=rep('blank', 2), fill=c('lightblue', 'green'), cex=0, cat.cex=0,margin=0.05)
dev.off()

pdf('venn_library_peptides.pdf', width=10, height=10)
draw.pairwise.venn(57671, 54449, 36846, category=c('standard search engines - Chapter 5', 'MSFraggerGUI - Modifications'), lty=rep('blank', 2), fill=c('lightblue', 'green'), cat.pos=c(-25,100), cex=3, cat.cex=3,margin=0.05)
dev.off()

pdf('venn_library_peptides_nolabels.pdf', width=10, height=10)
draw.pairwise.venn(57671, 54449, 36846, lty=rep('blank', 2), fill=c('lightblue', 'green'), cex=0, cat.cex=0,margin=0.05)
dev.off()

###############################################

## Step 8: Filter quantitative matrix only (at least 20% measurements)

missing <- NULL
for (i in 1:nrow(data_final)) {
tmp <- length(which(is.na(as.matrix(data_final[i,3:100]))))
missing <- c(missing, tmp)
}

data_filtered <- data_final[which(missing<78.4),]
dim(data_filtered)
[1] 4435  100
length(unique(data_filtered$Protein))
[1] 1433

# save for mapDIA

mapDIA_input <- cbind(data_filtered$Protein, data_filtered$Peptide, data_filtered[,3:100])
colnames(mapDIA_input)[1:2] <- c('ProteinName', 'PeptideSequence')
write.table(mapDIA_input, file="IPF_mapDIA_input.tsv", quote=F, row.names=F, sep="\t")
# feed into mapDIA


## Step 9: Make systematic t-test for differential expression

ttest_input <- as.matrix(data_filtered[,3:100])

ttest_fc <-NULL
ttest_p <- NULL
for (i in 1:nrow(ttest_input)){
tmp1 <- ttest_input[i,1:54]
tmp2 <- ttest_input[i,55:98]
tmp3 <- t.test(tmp1, tmp2)$p.value
tmp4 <- mean(tmp1, na.rm=T)-mean(tmp2, na.rm=T)
ttest_p <- c(ttest_p, tmp3)
ttest_fc <- c(ttest_fc, tmp4)
}

length(which(ttest_p<0.05))
[1] 175
length(which(ttest_p<0.01))
[1] 31

ttest_p_adj <- p.adjust(ttest_p, method = 'fdr')
length(which(ttest_p_adj<0.05))
[1] 0

ttest_result <- data.frame(Peptide=data_filtered$Peptide,Protein=data_filtered$Protein, p_value=ttest_p, adj_p_value=ttest_p_adj, log2_fc=ttest_fc)

pdf('Volcano.pdf', width=12, height=12)
with(ttest_result, plot(log2_fc, -log10(p_value), pch=20, cex=2, cex.lab=0.1, yaxt='n', cex.axis=2.2, xlab='', ylab='', xlim=c(-1.3, 1.3)))
with(subset(ttest_result, p_value<0.05 & log2_fc>log2(1.3)), points(log2_fc, -log10(p_value), pch=20, cex=2.2, col='red'))
with(subset(ttest_result, p_value<0.05 & log2_fc<(-log2(1.3))), points(log2_fc, -log10(p_value), pch=20, cex=2.2, col='dodgerblue'))
axis(2, las=2, cex.axis=2.2)
abline(v=c(log2(1.3),-log2(1.3)), col=c("darkgray","darkgray"), lty=c(2,2))
abline(h=-log10(0.05), col="darkgray", lty=2)
dev.off()



write.table(ttest_result[which(ttest_result$p_value<0.05),], file='significant_peptides_ttest_20181012.tsv', quote=F, row.names=F, sep="\t")
save(data_filtered, ttest_result, file='CPTAC_IPF_quantitative_analysis_20181012.Rdata')


## Step 10: Make binary matrix and fisher's exact test for differential appearance of specific modifications

# make binary matrix
data_bi <- as.matrix(data_final[,3:100])
data_bi[which(!is.na(data_bi))]<- 1
data_bi[which(is.na(data_bi))]<- 0

rownames(data_bi) <- data_final[,1]

# make fishers test on subsets of data

condition <- c(rep('HRD', 54),rep('non-HRD', 44))

a <- data_bi[grep('\\.\\(Acetyl\\)', rownames(data_bi)),]
fishers_p_nacetyl <- NULL
for (i in 1:nrow(a)){
tmp <- a[i,]
tmp2 <- fisher.test(table(tmp,condition))$p.value
fishers_p_nacetyl <- c(fishers_p_nacetyl, tmp2)
}
names(fishers_p_nacetyl) <- rownames(a)

# .(Acetyl)MDIAIHHPWIR_119 -> P02511
                     

a <- data_bi[grep('K\\(Acetyl\\)', rownames(data_bi)),]
fishers_p_kacetyl <- NULL
for (i in 1:nrow(a)){
tmp <- a[i,]
tmp2 <- fisher.test(table(tmp,condition))$p.value
fishers_p_kacetyl <- c(fishers_p_kacetyl, tmp2)
}


a <- data_bi[grep('\\(Phospho\\)', rownames(data_bi)),]
fishers_p_phospho <- NULL
for (i in 1:nrow(a)){
tmp <- a[i,]
tmp2 <- fisher.test(table(tmp,condition))$p.value
fishers_p_phospho <- c(fishers_p_phospho, tmp2)
}
names(fishers_p_phospho) <- rownames(a)

# GFPGPPGPDGLPGS(Phospho)MGPPGTPSVDHGFLVTR_92120 -> P02462


save(data_bi, fishers_p_nacetyl, fishers_p_kacetyl, fishers_p_phospho, file='CPTAC_IPF_fisher_analysis.Rdata')



