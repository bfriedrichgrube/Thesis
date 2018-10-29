## Correlation analysis with 37pat and 2 ctrl ##

## Step 1: Prepare R session 

setwd('Documents/Biomarker_Cancer/ProCOC/ProCOC_final/ProCOC2_correlation_analysis_new/')

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


## Step 2: Load data 

data <- read.delim('../ProCOC_mapDIA/ProCOC2_56pat_mapDIA_1pep/fragments_for_protein_quantification.txt', sep='\t', dec='.', as.is=T, header=T)
dim(data)
[1] 30276    61

## Step 3: Infer proteins from mapDIA output

#remove empty lines of filtered fragments in mapDIA
data_filtered <- data[1:8040, ]
dim(data_filtered)
[1] 8040   61

# filter for top3 fragment
data_filtered$fragment_sum <- apply(as.matrix(data_filtered[,4:61]), 1, sum, na.rm=T)

l <- data_filtered[order(data_filtered$fragment_sum, decreasing=T),]
d <- by(l, l$Peptide, head, n=3)
data_top3 <- Reduce(rbind, d)
dim(data_top3)
[1] 4203   62

#aggregate fragments to peptides -> sum
length(unique(data_top3$Peptide))
[1] 1401
data_peptide <- aggregate(data_top3[,4:61], by=list(data_top3$Protein, data_top3$Peptide), FUN=sum)
dim(data_peptide)
[1] 1401   60

#aggregate peptides to proteins -> average
length(unique(data_peptide$Group.1))
[1] 177

data_protein <- aggregate(as.matrix(data_peptide[,3:60]), by=list(data_peptide$Group.1), FUN=mean, na.rm=T)
dim(data_protein)
[1] 177  59

rownames(data_protein) <- data_protein[,1]
data_protein[,1]<- NULL
dim(data_protein)
[1] 177  58
class(data_protein)
[1] "data.frame"

#log transformation
ProCOC_log <- log2(as.matrix(data_protein))
ProCOC_log[which(is.infinite(ProCOC_log))]<- NA



# check replicates 
pdf('replicate_correlation_ProCOC2_20180827.pdf', height=4.5, width=10)
	par(mfrow=c(1,2),mai=c(1,1,0.1,0.5),cex.axis=1.5, las=1, mgp = c(3, 1.5, 0), cex.lab=1.5)
	plot(ProCOC_log[,1], ProCOC_log[,57], xlim=c(8,21), ylim=c(8,21), xlab=colnames(ProCOC_log)[1], ylab=colnames(ProCOC_log)[57])
		text(18,9.5, labels=paste("R=",round(cor(ProCOC_log[,1], ProCOC_log[,57], use='pairwise.complete.obs'),digits=4), sep=""), cex=1.5)
	plot(ProCOC_log[,11], ProCOC_log[,58], xlim=c(8,21), ylim=c(8,21), xlab=colnames(ProCOC_log)[11], ylab=colnames(ProCOC_log)[58])
		text(18,9.5, labels=paste("R=",round(cor(ProCOC_log[,11], ProCOC_log[,58], use='pairwise.complete.obs'),digits=4), sep=""), cex=1.5)
dev.off()

#check replicates in clustering 
replicate <- c(1,0,0,0,0,0, 0,0,0, 0,2,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0 ,0,0,0, 0,0,0, 0,0,0, 0,0,0,  0,0,0, 0,0,0,  0,0,0, 0,0,0, 0,0, 1,2)
anno_col <- data.frame(replicate)
rownames(anno_col)<- colnames(ProCOC_log)
library(pheatmap)
pheatmap(ProCOC_log, show_rownames=F, show_colnames=F, filename='heatmap_protein2_data_20180828.pdf',
border_color=NA, cluster_rows=F, cluster_cols=T, annotation_col=anno_col)


## Step 4: Internal fold change calculation

ProCOC_mean  <- apply(ProCOC_log[,1:56], 2, mean, na.rm=T)
length(ProCOC_mean)
[1] 56

ProCOC_log2FC <- matrix(0, nrow(ProCOC_log),56)

for (i in 1:56){
	for (j in 1:nrow(ProCOC_log)){
		tmp <- ProCOC_log[j,i] - ProCOC_mean[i]

		ProCOC_log2FC[j, i] <- tmp
	}
}

rownames(ProCOC_log2FC) <- rownames(ProCOC_log)
colnames(ProCOC_log2FC) <- colnames(ProCOC_log)[1:56]

ProCOC_log2FC <- ProCOC_log2FC[1:175,]

## Step 5: Correlation analysis

gleason <- c(9,9,9,9,8,9,9,9,9,9,9,9,9,9,9,9,9,9,9,8,9,9,9,9,9,9,9,9,9,9,9,9,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6)

prot_gleason <- sapply(seq_len(nrow(ProCOC_log2FC)), function(i)
	cor.test(ProCOC_log2FC[i,],gleason, method="pearson"))
padjust.prot <- p.adjust(prot_gleason[3,], method='fdr')
prot_gleason <- rbind(prot_gleason, padjust.prot)

colnames(prot_gleason)<- rownames(ProCOC_log2FC)

#extract pearson phi
a <- as.numeric(prot_gleason[4,])
names(a)<- colnames(prot_gleason)

#get top 10 positive and negative correlating proteins

> head(sort(a, decreasing=T), 10)
   O95998    P04070    Q9UBG0    P01861    O95497    Q9UHG3    P11597    P35542 
0.4149158 0.3141852 0.2635677 0.2279125 0.2196640 0.2189273 0.2128925 0.1999726 
   P08697    Q92954 
0.1978848 0.1868601 
> tail(sort(a, decreasing=T), 10)
    P80108     Q6EMK4     Q9NUL3     P05164     Q8N6C8     P01042     P20160 
-0.2085973 -0.2206379 -0.2212457 -0.2304781 -0.2360330 -0.2553962 -0.2661831 
    P05160     P20851     P02788 
-0.2727079 -0.2769309 -0.3056685 
> names(head(sort(a, decreasing=T), 10))
 [1] "O95998" "P04070" "Q9UBG0" "P01861" "O95497" "Q9UHG3" "P11597" "P35542"
 [9] "P08697" "Q92954"
> names(tail(sort(a, decreasing=T), 10))
 [1] "P80108" "Q6EMK4" "Q9NUL3" "P05164" "Q8N6C8" "P01042" "P20160" "P05160"
 [9] "P20851" "P02788"


## Step 7: Check if correlation is better than random 

mean(abs(a))
[1] 0.1017227

# mix gleason scores randomly 100-times
mean_distribution <- NULL
for(i in 1:100){
	tmp_gleason <- sample(gleason)
	cor_tmp <- sapply(seq_len(nrow(ProCOC_log2FC)), function(i)
	cor.test(ProCOC_log2FC[i,],tmp_gleason, method="pearson"))
	tmp_mean <- mean(abs(as.numeric(cor_tmp[4,])))
	mean_distribution <- c(mean_distribution, tmp_mean)
}

length(mean_distribution)
[1] 100
boxplot(mean_distribution)


# test if distribution of corrleations is significantly different from real gleason assignments (one-sample t.test)

t.test(mean_distribution, mu=mean(abs(a)))

	One Sample t-test

data:  mean_distribution
t = 6.5863, df = 99, p-value = 2.188e-09
alternative hypothesis: true mean is not equal to 0.1017227
95 percent confidence interval:
 0.1081468 0.1136864
sample estimates:
mean of x 
0.1109166 

## yes random correlations are significantly lower than real gleason assignments 

> prot_gleason[c(3,4,10),names(head(sort(a, decreasing=T), 10))]
             O95998      P04070     Q9UBG0     P01861     O95497    Q9UHG3   
p.value      0.006291711 0.02964801 0.04967731 0.09114515 0.1104921 0.1307043
estimate     0.4149158   0.3141852  0.2635677  0.2279125  0.219664  0.2189273
padjust.prot 0.9907155   0.9907155  0.9907155  0.9907155  0.9907155 0.9907155
             P11597    P35542    P08697    Q92954   
p.value      0.1186328 0.147115  0.2149001 0.1679025
estimate     0.2128925 0.1999726 0.1978848 0.1868601
padjust.prot 0.9907155 0.9907155 0.9907155 0.9907155
> prot_gleason[c(3,4,10),names(tail(sort(a, decreasing=T), 10))]
             P80108     Q6EMK4     Q9NUL3     P05164     Q8N6C8    P01042    
p.value      0.1228747  0.1550913  0.1187018  0.10371    0.0954195 0.05984877
estimate     -0.2085973 -0.2206379 -0.2212457 -0.2304781 -0.236033 -0.2553962
padjust.prot 0.9907155  0.9907155  0.9907155  0.9907155  0.9907155 0.9907155 
             P20160     P05160     P20851     P02788    
p.value      0.07716072 0.04397174 0.03880756 0.0291616 
estimate     -0.2661831 -0.2727079 -0.2769309 -0.3056685
padjust.prot 0.9907155  0.9907155  0.9907155  0.9907155 


## Step 8: Save data

#save protein inference

save(data_protein, ProCOC_log, file='ProCOC2_protein_level_20180827.Rdata')
save(data_top3, data_filtered, data_peptide, file='ProCOC2_protein_inference_workflow_20180827.Rdata')

save(a, names_prot2, ProCOC_log2FC, mean_distribution, prot_gleason, file='ProCOC2_correlation_analysis_20180827.Rdata')

pdf('overlap_measured_prots_20180903.pdf', width=10, height=10)
draw.pairwise.venn(194, 177, 170, category=c('Sub-cohort 1', 'Sub-cohort 2'), lty=rep('blank', 2), fill=c('lightblue', 'green'), cat.cex=3, cat.pos=c(-140,140), cex=3)
dev.off()

pdf('overlap_measured_prots_20180903.pdf', width=10, height=10)
> draw.pairwise.venn(192, 175, 168, lty=rep('blank', 2), fill=c('lightblue', 'green'), cex=0)
(polygon[GRID.polygon.131], polygon[GRID.polygon.132], polygon[GRID.polygon.133], polygon[GRID.polygon.134], text[GRID.text.135], text[GRID.text.136], lines[GRID.lines.137], text[GRID.text.138], text[GRID.text.139], text[GRID.text.140]) 
> dev.off()



O95998
P04070
Q9UBG0
P01861
O95497
Q9UHG3
P11597
P35542
P08697
Q92954
P80108
Q6EMK4
Q9NUL3 
P05164 
Q8N6C8
P01042
P20160
P05160
P20851
P02788

