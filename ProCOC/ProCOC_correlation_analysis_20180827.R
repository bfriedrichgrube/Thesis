## Correlation analysis with 37pat and 2 ctrl ##

## Step 1: Prepare R session 

setwd('Documents/Biomarker_Cancer/ProCOC/ProCOC_final/ProCOC1_correlation_analysis_new/')

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

data <- read.delim('../ProCOC_mapDIA/ProCOC1_/fragments_for_protein_quantification.txt', sep='\t', dec='.', as.is=T, header=T)
dim(data)
[1] 32088    62
## Step 3: Infer proteins from mapDIA output

#remove empty lines of filtered fragments in mapDIA
data_filtered <- data[1:8673, ]
dim(data_filtered)
[1] 8673   62

# filter for top3 fragment
data_filtered$fragment_sum <- apply(as.matrix(data_filtered[,4:62]), 1, sum, na.rm=T)

l <- data_filtered[order(data_filtered$fragment_sum, decreasing=T),]
d <- by(l, l$Peptide, head, n=3)
data_top3 <- Reduce(rbind, d)
dim(data_top3)
[1] 4503   63

#aggregate fragments to peptides -> sum
length(unique(data_top3$Peptide))
[1] 1501
data_peptide <- aggregate(data_top3[,4:62], by=list(data_top3$Protein, data_top3$Peptide), FUN=sum)
dim(data_peptide)
[1] 1501   61

#aggregate peptides to proteins -> average
length(unique(data_peptide$Group.1))
[1] 194

data_protein <- aggregate(as.matrix(data_peptide[,3:61]), by=list(data_peptide$Group.1), FUN=mean, na.rm=T)
dim(data_protein)
[1] 194  60

rownames(data_protein) <- data_protein[,1]
data_protein[,1]<- NULL
dim(data_protein)

class(data_protein)
[1] "data.frame"

#log transformation
ProCOC_log <- log2(as.matrix(data_protein))
ProCOC_log[which(is.infinite(ProCOC_log))]<- NA



# check replicates 
pdf('batch_correlation_ProCOC_20180827.pdf', height=20, width=25)
	par(mfrow=c(4,5), mai=c(0.7, 0.8, 0.7, 0.8), cex.axis=2.5, las=1, mgp = c(4.2, 1.5, 0), cex.lab=2.5)
	plot(ProCOC_log[,1], ProCOC_log[,2], xlim=c(8,19), ylim=c(8,19), xlab=colnames(ProCOC_log)[1], ylab=colnames(ProCOC_log)[2])
		text(17,9.5, labels=paste("R=",round(cor(ProCOC_log[,1], ProCOC_log[,2], use='pairwise.complete.obs'),digits=4), sep=""), cex=2.5)
	plot(ProCOC_log[,3], ProCOC_log[,4], xlim=c(8,19), ylim=c(8,19), xlab=colnames(ProCOC_log)[3], ylab=colnames(ProCOC_log)[4])
		text(17,9.5, labels=paste("R=",round(cor(ProCOC_log[,3], ProCOC_log[,4], use='pairwise.complete.obs'),digits=4), sep=""), cex=2.5)
	plot(ProCOC_log[,8], ProCOC_log[,9], xlim=c(8,19), ylim=c(8,19), xlab=colnames(ProCOC_log)[8], ylab=colnames(ProCOC_log)[9])
		text(17,9.5, labels=paste("R=",round(cor(ProCOC_log[,8], ProCOC_log[,9], use='pairwise.complete.obs'),digits=4), sep=""), cex=2.5)
	plot(ProCOC_log[,11], ProCOC_log[,12], xlim=c(8,19), ylim=c(8,19), xlab=colnames(ProCOC_log)[11], ylab=colnames(ProCOC_log)[12])
		text(17,9.5, labels=paste("R=",round(cor(ProCOC_log[,11], ProCOC_log[,12], use='pairwise.complete.obs'),digits=4), sep=""), cex=2.5)
	plot(ProCOC_log[,13], ProCOC_log[,14], xlim=c(8,19), ylim=c(8,19), xlab=colnames(ProCOC_log)[13], ylab=colnames(ProCOC_log)[14])
		text(17,9.5, labels=paste("R=",round(cor(ProCOC_log[,13], ProCOC_log[,14], use='pairwise.complete.obs'),digits=4), sep=""), cex=2.5)
	plot(ProCOC_log[,15], ProCOC_log[,16], xlim=c(8,19), ylim=c(8,19), xlab=colnames(ProCOC_log)[15], ylab=colnames(ProCOC_log)[16])
		text(17,9.5, labels=paste("R=",round(cor(ProCOC_log[,15], ProCOC_log[,16], use='pairwise.complete.obs'),digits=4), sep=""), cex=2.5)
	plot(ProCOC_log[,22], ProCOC_log[,23], xlim=c(8,19), ylim=c(8,19), xlab=colnames(ProCOC_log)[22], ylab=colnames(ProCOC_log)[23])
		text(17,9.5, labels=paste("R=",round(cor(ProCOC_log[,22], ProCOC_log[,23], use='pairwise.complete.obs'),digits=4), sep=""), cex=2.5)
	plot(ProCOC_log[,27], ProCOC_log[,28], xlim=c(8,19), ylim=c(8,19), xlab=colnames(ProCOC_log)[27], ylab=colnames(ProCOC_log)[28])
		text(17,9.5, labels=paste("R=",round(cor(ProCOC_log[,27], ProCOC_log[,28], use='pairwise.complete.obs'),digits=4), sep=""), cex=2.5)
	plot(ProCOC_log[,30], ProCOC_log[,31], xlim=c(8,19), ylim=c(8,19), xlab=colnames(ProCOC_log)[30], ylab=colnames(ProCOC_log)[31])
		text(17,9.5, labels=paste("R=",round(cor(ProCOC_log[,30], ProCOC_log[,31], use='pairwise.complete.obs'),digits=4), sep=""), cex=2.5)
	plot(ProCOC_log[,35], ProCOC_log[,36], xlim=c(8,19), ylim=c(8,19), xlab=colnames(ProCOC_log)[35], ylab=colnames(ProCOC_log)[36])
		text(17,9.5, labels=paste("R=",round(cor(ProCOC_log[,35], ProCOC_log[,36], use='pairwise.complete.obs'),digits=4), sep=""), cex=2.5)
	plot(ProCOC_log[,37], ProCOC_log[,38], xlim=c(8,19), ylim=c(8,19), xlab=colnames(ProCOC_log)[37], ylab=colnames(ProCOC_log)[38])
		text(17,9.5, labels=paste("R=",round(cor(ProCOC_log[,37], ProCOC_log[,38], use='pairwise.complete.obs'),digits=4), sep=""), cex=2.5)
	plot(ProCOC_log[,39], ProCOC_log[,40], xlim=c(8,19), ylim=c(8,19), xlab=colnames(ProCOC_log)[39], ylab=colnames(ProCOC_log)[40])
		text(17,9.5, labels=paste("R=",round(cor(ProCOC_log[,39], ProCOC_log[,40], use='pairwise.complete.obs'),digits=4), sep=""), cex=2.5)
	plot(ProCOC_log[,42], ProCOC_log[,43], xlim=c(8,19), ylim=c(8,19), xlab=colnames(ProCOC_log)[42], ylab=colnames(ProCOC_log)[43])
		text(17,9.5, labels=paste("R=",round(cor(ProCOC_log[,42], ProCOC_log[,43], use='pairwise.complete.obs'),digits=4), sep=""), cex=2.5)
	plot(ProCOC_log[,44], ProCOC_log[,45], xlim=c(8,19), ylim=c(8,19), xlab=colnames(ProCOC_log)[44], ylab=colnames(ProCOC_log)[45])
		text(17,9.5, labels=paste("R=",round(cor(ProCOC_log[,44], ProCOC_log[,45], use='pairwise.complete.obs'),digits=4), sep=""), cex=2.5)
	plot(ProCOC_log[,46], ProCOC_log[,47], xlim=c(8,19), ylim=c(8,19), xlab=colnames(ProCOC_log)[46], ylab=colnames(ProCOC_log)[47])
		text(17,9.5, labels=paste("R=",round(cor(ProCOC_log[,46], ProCOC_log[,47], use='pairwise.complete.obs'),digits=4), sep=""), cex=2.5)
	plot(ProCOC_log[,48], ProCOC_log[,49], xlim=c(8,19), ylim=c(8,19), xlab=colnames(ProCOC_log)[48], ylab=colnames(ProCOC_log)[49])
		text(17,9.5, labels=paste("R=",round(cor(ProCOC_log[,48], ProCOC_log[,49], use='pairwise.complete.obs'),digits=4), sep=""), cex=2.5)
	plot(ProCOC_log[,50], ProCOC_log[,51], xlim=c(8,19), ylim=c(8,19), xlab=colnames(ProCOC_log)[50], ylab=colnames(ProCOC_log)[51])
		text(17,9.5, labels=paste("R=",round(cor(ProCOC_log[,50], ProCOC_log[,51], use='pairwise.complete.obs'),digits=4), sep=""), cex=2.5)
	plot(ProCOC_log[,52], ProCOC_log[,53], xlim=c(8,19), ylim=c(8,19), xlab=colnames(ProCOC_log)[52], ylab=colnames(ProCOC_log)[53])
		text(17,9.5, labels=paste("R=",round(cor(ProCOC_log[,52], ProCOC_log[,53], use='pairwise.complete.obs'),digits=4), sep=""), cex=2.5)
	plot(ProCOC_log[,54], ProCOC_log[,55], xlim=c(8,19), ylim=c(8,19), xlab=colnames(ProCOC_log)[54], ylab=colnames(ProCOC_log)[55])
		text(17,9.5, labels=paste("R=",round(cor(ProCOC_log[,54], ProCOC_log[,55], use='pairwise.complete.obs'),digits=4), sep=""), cex=2.5)
	plot(ProCOC_log[,56], ProCOC_log[,57], xlim=c(8,19), ylim=c(8,19), xlab=colnames(ProCOC_log)[56], ylab=colnames(ProCOC_log)[57])
		text(17,9.5, labels=paste("R=",round(cor(ProCOC_log[,56], ProCOC_log[,57], use='pairwise.complete.obs'),digits=4), sep=""), cex=2.5)
dev.off()	


pdf('replicate_correlation_ProCOC_20180827.pdf', height=4.5, width=10)
	par(mfrow=c(1,2),mai=c(1,1,0.1,0.5),cex.axis=1.5, las=1, mgp = c(3, 1.5, 0), cex.lab=1.5)
	plot(ProCOC_log[,5], ProCOC_log[,58], xlim=c(8,19), ylim=c(8,19), xlab=colnames(ProCOC_log)[5], ylab=colnames(ProCOC_log)[58])
		text(17,9.5, labels=paste("R=",round(cor(ProCOC_log[,5], ProCOC_log[,58], use='pairwise.complete.obs'),digits=4), sep=""), cex=1.5)
	plot(ProCOC_log[,50], ProCOC_log[,59], xlim=c(8,19), ylim=c(8,19), xlab=colnames(ProCOC_log)[50], ylab=colnames(ProCOC_log)[59])
		text(17,9.5, labels=paste("R=",round(cor(ProCOC_log[,50], ProCOC_log[,59], use='pairwise.complete.obs'),digits=4), sep=""), cex=1.5)
dev.off()

#check replicates in clustering 
batch <- c(1,2,1,2,1,1,2,1,2,1,1,2,1,2,1,2,2,1,1,1,1,1,2,1,1,1,1,2,1,1,2,1,2,2,1,2,1,2,1,2,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,1)
length(batch)
[1] 59
replicate <- c(0,0,0,0,1,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0 ,0,0,0, 0,0,0, 0,0,0, 0,0,0,  0,0,0, 0,0,0,  0,2,0, 0,0,0, 0,0,0, 1,2)
anno_col <- data.frame(batch, replicate)
rownames(anno_col)<- colnames(ProCOC_log)
library(pheatmap)
pheatmap(ProCOC_log, show_rownames=F, show_colnames=F, filename='heatmap_protein_data_20180827.pdf',
border_color=NA, cluster_rows=F, cluster_cols=T, annotation_col=anno_col)


## Step 4: Internal fold change calculation

ProCOC_mean  <- apply(ProCOC_log[,1:57], 2, mean, na.rm=T)
length(ProCOC_mean)
[1] 57

ProCOC_log2FC <- matrix(0, nrow(ProCOC_log),57)

for (i in 1:57){
	for (j in 1:nrow(ProCOC_log)){
		tmp <- ProCOC_log[j,i] - ProCOC_mean[i]

		ProCOC_log2FC[j, i] <- tmp
	}
}

rownames(ProCOC_log2FC) <- rownames(ProCOC_log)
colnames(ProCOC_log2FC) <- colnames(ProCOC_log)[1:57]

####       Using all samples seperately
#gleason <- c(9, 9, 9, 9, 8, 8, 9, 9, 9, 8, 9, 9, 9, 9, 9, 9, 9, 7, 7, 7.5, 7, 7.5, 7.5, 7, 7, 
#7, 7, 7, 7, 7, 7, 7, 6, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6)

#prot_gleason <- sapply(seq_len(nrow(ProCOC_log2FC)), function(i)
#cor.test(ProCOC_log2FC[i,],gleason, method="pearson"))
#padjust.prot <- p.adjust(prot_gleason[3,], method='fdr')
#prot_gleason <- rbind(prot_gleason, padjust.prot)




#### Step 5: Combining replicates by averaging

ProCOC_log2FC_reduced <- matrix(0, 194, 37)
ProCOC_log2FC_reduced[,1] <- apply(ProCOC_log2FC[,1:2], 1, mean, na.rm=T)
ProCOC_log2FC_reduced[,2] <- apply(ProCOC_log2FC[,3:4], 1, mean, na.rm=T)
ProCOC_log2FC_reduced[,3:5] <- ProCOC_log2FC[,5:7]
ProCOC_log2FC_reduced[,6] <- apply(ProCOC_log2FC[,8:9], 1, mean, na.rm=T)
ProCOC_log2FC_reduced[,7] <- ProCOC_log2FC[,10]
ProCOC_log2FC_reduced[,8] <- apply(ProCOC_log2FC[,11:12], 1, mean, na.rm=T)
ProCOC_log2FC_reduced[,9] <- apply(ProCOC_log2FC[,13:14], 1, mean, na.rm=T)
ProCOC_log2FC_reduced[,10] <- apply(ProCOC_log2FC[,15:16], 1, mean, na.rm=T)
ProCOC_log2FC_reduced[,11:15] <- ProCOC_log2FC[,17:21]
ProCOC_log2FC_reduced[,16] <- apply(ProCOC_log2FC[,22:23], 1, mean, na.rm=T)
ProCOC_log2FC_reduced[,17:19] <- ProCOC_log2FC[,24:26]
ProCOC_log2FC_reduced[,20] <- apply(ProCOC_log2FC[,27:28], 1, mean, na.rm=T)
ProCOC_log2FC_reduced[,21] <- ProCOC_log2FC[,29]
ProCOC_log2FC_reduced[,22] <- apply(ProCOC_log2FC[,30:31], 1, mean, na.rm=T)
ProCOC_log2FC_reduced[,23:25] <- apply(ProCOC_log2FC[,32:34], 1, mean, na.rm=T)
ProCOC_log2FC_reduced[,26] <- apply(ProCOC_log2FC[,35:36], 1, mean, na.rm=T)
ProCOC_log2FC_reduced[,27] <- apply(ProCOC_log2FC[,37:38], 1, mean, na.rm=T)
ProCOC_log2FC_reduced[,28] <- apply(ProCOC_log2FC[,39:40], 1, mean, na.rm=T)
ProCOC_log2FC_reduced[,29] <- ProCOC_log2FC[,41]
ProCOC_log2FC_reduced[,30] <- apply(ProCOC_log2FC[,42:43], 1, mean, na.rm=T)
ProCOC_log2FC_reduced[,31] <- apply(ProCOC_log2FC[,44:45], 1, mean, na.rm=T)
ProCOC_log2FC_reduced[,32] <- apply(ProCOC_log2FC[,46:47], 1, mean, na.rm=T)
ProCOC_log2FC_reduced[,33] <- apply(ProCOC_log2FC[,48:49], 1, mean, na.rm=T)
ProCOC_log2FC_reduced[,34] <- apply(ProCOC_log2FC[,50:51], 1, mean, na.rm=T)
ProCOC_log2FC_reduced[,35] <- apply(ProCOC_log2FC[,52:53], 1, mean, na.rm=T)
ProCOC_log2FC_reduced[,36] <- apply(ProCOC_log2FC[,54:55], 1, mean, na.rm=T)
ProCOC_log2FC_reduced[,37] <- apply(ProCOC_log2FC[,56:57], 1, mean, na.rm=T)

rownames(ProCOC_log2FC_reduced)<- rownames(ProCOC_log2FC)


colnames(ProCOC_log2FC_reduced) <- c("high_100", "high_201","high_21","high_218", "high_283", "high_285", "high_289", "high_291", "high_299", "high_339", "high_55.2", "intermediate_22", "intermediate_25", "intermediate_26", "intermediate_296", "intermediate_31",   "intermediate_319", "intermediate_320", "intermediate_49", "intermediate_57", "intermediate_79", "intermediate_83", "intermediate_87", "low_12", "low_125", "low_126", "low_137", "low_148", "low_149", "low_163", "low_222",  "low_276", "low_317", "low_355", "low_46", "low_53", "low_92")       

## Step 6: Correlation analysis

gleason <- c(9, 9, 8, 8, 9, 9, 8, 9, 9, 9, 9, 7, 7, 7.5, 7, 7.5, 7, 7, 7, 7, 7, 7, 7, 6, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6)

prot_gleason <- sapply(seq_len(nrow(ProCOC_log2FC_reduced)), function(i)
	cor.test(ProCOC_log2FC_reduced[i,],gleason, method="pearson"))
padjust.prot <- p.adjust(prot_gleason[3,], method='fdr')
prot_gleason <- rbind(prot_gleason, padjust.prot)

colnames(prot_gleason)<- rownames(ProCOC_log2FC_reduced)

#extract pearson phi
a <- as.numeric(prot_gleason[4,])
names(a)<- colnames(prot_gleason)

#get top 10 positive and negative correlating proteins
head(sort(a, decreasing=T), 10)
   P58166    Q4VNC0    P07339    Q9UBG0    Q86VB7    P13473    P03951    P23142 
0.5532087 0.5459124 0.5290736 0.4885483 0.4640167 0.4556688 0.4535990 0.4506214 
   O15204    P01033 
0.4192939 0.4115924
tail(sort(a, decreasing=T), 10)
    Q14624     P35555     P35858     P22792     P01009     P04003     P07333 
-0.3616691 -0.3633026 -0.3650601 -0.3682092 -0.3757177 -0.3763515 -0.4005627 
    P29622     P02748     P01011 
-0.4307700 -0.4421052 -0.5478691 
   
> names(head(sort(a, decreasing=T), 10))
 [1] "P58166" "Q4VNC0" "P07339" "Q9UBG0" "Q86VB7" "P13473" "P03951" "P23142"
 [9] "O15204" "P01033"
> names(tail(sort(a, decreasing=T), 10))
 [1] "Q14624" "P35555" "P35858" "P22792" "P01009" "P04003" "P07333" "P29622"
 [9] "P02748" "P01011"

prot_gleason[c(3,4,10),names(head(sort(a, decreasing=T), 10))]
             P58166       Q4VNC0       P07339      Q9UBG0      Q86VB7     
p.value      0.0005672377 0.0004742355 0.001083898 0.002155396 0.004363388
estimate     0.5532087    0.5459124    0.5290736   0.4885483   0.4640167  
padjust.prot 0.03668137   0.03668137   0.05256904  0.08362936  0.1105441  
             P13473     P03951      P23142      O15204     P01033    
p.value      0.00459492 0.004807546 0.005128333 0.01091235 0.01402734
estimate     0.4556688  0.453599    0.4506214   0.4192939  0.4115924 
padjust.prot 0.1105441  0.1105441   0.1105441   0.1764163  0.1814202 
> prot_gleason[c(3,4,10),names(tail(sort(a, decreasing=T), 10))]
             Q14624     P35555     P35858     P22792     P01009     P04003    
p.value      0.02784261 0.03193891 0.02630616 0.02494268 0.02192632 0.02367795
estimate     -0.3616691 -0.3633026 -0.3650601 -0.3682092 -0.3757177 -0.3763515
padjust.prot 0.2077487  0.2117192  0.205441   0.205441   0.205441   0.205441  
             P07333     P29622      P02748     P01011      
p.value      0.0255473  0.007777966 0.0061499  0.0004481817
estimate     -0.4005627 -0.43077    -0.4421052 -0.5478691  
padjust.prot 0.205441   0.137175    0.1193081  0.03668137  



## Step 7: Check if correlation is better than random 

mean(abs(a))
[1] 0.1880399


# mix gleason scores randomly 100-times
mean_distribution <- NULL
for(i in 1:100){
	tmp_gleason <- sample(gleason)
	cor_tmp <- sapply(seq_len(nrow(ProCOC_log2FC_reduced)), function(i)
	cor.test(ProCOC_log2FC_reduced[i,],tmp_gleason, method="pearson"))
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
t = -34.145, df = 99, p-value < 2.2e-16
alternative hypothesis: true mean is not equal to 0.1880399
95 percent confidence interval:
 0.1292266 0.1356867
sample estimates:
mean of x 
0.1324567 

## yes random corrleations are significantly lower than real gleason assignments 

## Step 8: Save data

#save protein inference

save(data_protein, ProCOC_log, file='ProCOC_protein_level_20180827.Rdata')
save(data_top3, data_filtered, data_peptide, file='ProCOC_protein_inference_workflow_20180827.Rdata')

save(a, names_prot1, ProCOC_log2FC, ProCOC_log2FC_reduced, mean_distribution, prot_gleason, file='ProCOC1_correlation_analysis_20180827.Rdata')

P58166
Q4VNC0
P07339
Q9UBG0
Q86VB7
P13473
P03951
P23142
O15204
P01033
Q14624
P35555
P35858
P22792
P01009
P04003
P07333
P29622
P02748
P01011


