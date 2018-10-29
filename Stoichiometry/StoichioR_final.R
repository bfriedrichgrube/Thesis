on crick

	
> setwd('Stoichio_20180912/')

library(StochioR)
int_pairs = read.table ("int_in_prot_compl_w_annot", sep = "\t", stringsAsFactors = FALSE, header = T)
load('Stoichio_input_20180912.Rdata')


data <- StochioR(read_expression_data = TRUE, data = quant_data_all,
        cols_with_reference_data = c(59:103),
        log_transformed = TRUE,
        ProteinPairs = int_pairs, Reference = "reference", FractionThr = 0.1,
        score_threshold = 3, n_samples_fraction = 0.2, n_samples_threshold = 12,
        perc_thr = 0.2, modif = 1)


read_expression_data
read_protein_pairs
variability_protein_ratios
perturbed_protein_ratios
The number of samples in which the ratio has to been above the threshold is: 11.6
[1] "Top 0.1 and 1% upper and lower values are: 6.45, 3.86, -6.46 and -3.6."
[1] "Top 1% of the absolute values for the modified z-scores is 4.52."
representative_pairs
Modif has been set to 2.

summary(data)
                           Length Class         Mode
QuantData                  1      ExpressionSet S4  
ProteinRatios              2      -none-        list
variability_protein_ratios 6      data.frame    list
perturbed_protein_ratios   4      -none-        list
representative_pairs       5      data.frame    list

biogrid <- read.delim2('coreComplexes.txt', sep='\t', header=T)
> a <- unlist(strsplit(data$representative_pairs[,1], '-')
+ )
length(a)
[1] 52
head(a)
[1] "P62841" "P42677" "P62917" "Q15287" "P08238" "Q15185"
b <- a[seq(1,length(a), by=2)]
head(b)
[1] "P62841" "P62917" "P08238" "Q9Y6M9" "Q15287" "O75533"
c <- a[seq(2,length(a), by=2)]
head(c)
[1] "P42677" "Q15287" "Q15185" "P28331" "O75533" "P41223"

d <- cbind(b,c)
head(d)
     b        c       
[1,] "P62841" "P42677"
[2,] "P62917" "Q15287"
[3,] "P08238" "Q15185"
[4,] "Q9Y6M9" "P28331"
[5,] "Q15287" "O75533"
[6,] "O75533" "P41223"

complex_known <- list()
for (i in 1:length(b)){

	tmp1 <- biogrid[grep(b[i], biogrid$subunits.UniProt.IDs.),1]
	tmp2 <- biogrid[grep(c[i], biogrid$subunits.UniProt.IDs.),1]
	tmp3 <- intersect(tmp1, tmp2)

	if (length(tmp3)==0) {
		complex_known[[i]] <- 'no complex known'
	}
	else {
		tmp4 <- biogrid[which(biogrid$ComplexID %in% tmp3),c(1:3,9:10)]
		complex_known[[i]] <- tmp4
	}
names(complex_known[[i]]) <- paste(b[i], c[i], collapse=' ')

}

tmp1 <- int_pairs[,1]
class(int_pairs)
[1] "data.frame"
tmp2 <- int_pairs[,2]
tmp <- c(tmp1, tmp2)
prot_mapped <- intersect(rownames(quant_data_all), tmp)
length(prot_mapped)
[1] 953
head(prot_mapped)
[1] "P42285" "Q8IX12" "Q7Z434" "Q9NTI5" "O76094" "Q7L2H7"
write.table(prot_mapped, file='prot_mapped.txt', quote=F, row.names=F)


check random pairing

set.seed(0)
> ProtA <- sample(int_pairs[,1])
> dim(int_pairs)
[1] 99329     3
> length(ProtA)
[1] 99329
> head(ProtA)
[1] "Q9UBU7" "Q9NTW7" "P42679" "P04843" "P52298" "P52815"
> int_pairs_random <- cbind(ProtA, int_pairs[,2:3])
> head(int_pairs_random)
   ProtA  ProtB          Complex_reported
1 Q9UBU7 Q9HBY0                  Reactome
2 Q9NTW7 Q86UR1                  Reactome
3 P42679 P63000                  Reactome
4 P04843 P19878 Int3D_Structure; Reactome
5 P52298 P15153                  Reactome
6 P52815 Q8IZP0       Int3D_Dom_dom_model


StoichioR does not run through with a random network -> no pairs identified. 


#check if siginificant samples are the same in similar complexes

data$representative_pairs[c(2,5,6,7,10,13,16,20,21,22),c(1,3)]
    Protein_pair
6  P62917-Q15287
4  Q15287-O75533
21 O75533-P41223
2  Q6P2Q9-Q53GS9
26 Q15287-O75643
27 Q15287-Q05048
16 Q15029-Q15287
24 Q15393-Q9BWJ5
8  Q6P2Q9-Q9UG63
10 P08865-Q15287
                                                                                                                          Significant_samples
6                       HRD_1,HRD_102,HRD_103,HRD_29,HRD_34,HRD_43,HRD_5,HRD_59,HRD_60,HRD_61,HRD_65,HRD_68,HRD_77,HRD_83,HRD_85,HRD_87,HRD_9
4                              HRD_1,HRD_102,HRD_103,HRD_29,HRD_34,HRD_43,HRD_5,HRD_59,HRD_60,HRD_61,HRD_65,HRD_77,HRD_83,HRD_85,HRD_87,HRD_9
21 HRD_1,HRD_102,HRD_103,HRD_23,HRD_29,HRD_34,HRD_43,HRD_44,HRD_5,HRD_61,HRD_63,HRD_68,HRD_69,HRD_70,HRD_72,HRD_82,HRD_83,HRD_84,HRD_87,HRD_9
2                            HRD_1,HRD_102,HRD_103,HRD_14,HRD_23,HRD_29,HRD_34,HRD_36,HRD_59,HRD_60,HRD_65,HRD_71,HRD_75,HRD_85,HRD_93,HRD_98
26                             HRD_1,HRD_102,HRD_103,HRD_29,HRD_34,HRD_43,HRD_5,HRD_59,HRD_60,HRD_61,HRD_65,HRD_77,HRD_83,HRD_85,HRD_87,HRD_9
27                                             HRD_1,HRD_29,HRD_34,HRD_43,HRD_5,HRD_59,HRD_60,HRD_61,HRD_65,HRD_77,HRD_83,HRD_85,HRD_87,HRD_9
16                                           HRD_1,HRD_102,HRD_29,HRD_34,HRD_43,HRD_59,HRD_60,HRD_61,HRD_65,HRD_77,HRD_83,HRD_85,HRD_87,HRD_9
24                                                   HRD_1,HRD_14,HRD_23,HRD_43,HRD_6,HRD_61,HRD_65,HRD_68,HRD_73,HRD_75,HRD_77,HRD_88,HRD_97
8                                                    HRD_1,HRD_103,HRD_24,HRD_34,HRD_36,HRD_43,HRD_5,HRD_61,HRD_65,HRD_82,HRD_83,HRD_87,HRD_9
10                                                   HRD_1,HRD_102,HRD_29,HRD_34,HRD_43,HRD_5,HRD_59,HRD_60,HRD_61,HRD_65,HRD_77,HRD_83,HRD_9


        HRD_1 10x
HRD_102 7x
HRD_103 6x 
HRD_14 2x
HRD_23 3x
HRD_24 1x
        HRD_29 8x
        HRD_34 9x
        HRD_43 9x
HRD_44 1x
HRD_5 7x 
HRD_59 7x
HRD_60 7x
HRD_61 9x
HRD_63 1x
        HRD_65 9x
HRD_68 3x 
HRD_69 1x
HRD_70 1x
HRD_71 1x
HRD_72 1x
HRD_77 7x
HRD_82 2x
        HRD_83 8x
HRD_84 1x
HRD_85 6x
        HRD_87 7x
HRD_88 1x
        HRD_9 8x
HRD_93 1x
HRD_97 1x
HRD_98 1x


#### Analysis of results - Plot ratios and distributions ####

## plot distribution boxplots of top spliceosomal interaction pairs

## RNPS1, SNRNP200, RPL8, SF3B1, BUD31, PRPF8, USP39

tmp1 <- quant_data_all['Q15287',] - quant_data_all['P62917',]
tmp2 <- quant_data_all['Q15287',] - quant_data_all['O75643',]
tmp3 <- quant_data_all['Q15287',] - quant_data_all['O75533',]
tmp4 <- quant_data_all['O75533',] - quant_data_all['P41223',]
tmp5 <- quant_data_all['Q6P2Q9',] - quant_data_all['Q53GS9',]


pdf('boxplot_distributions_20180924.pdf', width=20, height=15)
par(cex.axis=2.5, las=1, mar=c(25,10,1,1))
boxplot(tmp1[1:58],tmp1[59:103], tmp2[1:58],tmp2[59:103], tmp3[1:58],tmp3[59:103], tmp4[1:58],tmp4[59:103],tmp5[1:58],tmp5[59:103], notch=T, col=c("#E495A5","#E495A5","#BDAB66","#BDAB66","#65BC8C","#65BWarning message:#55B8D0","#C29DDE","#C29DDE"), border=rep(c('black', 'red'), times=5), boxlwd=4, xaxt='n')
text(c('RNPS1-RPL8 (HRD)', 'RNPS1-RPL8 (non-HRD)', 'RNPS1-SNRNP200 (HRD)','RNPS1-SNRNP200 (non-HRD)','RNPS1-SF3B1 (HRD)','RNPS1-SF3B1 (non-HRD)','SF3B1-BUD31 (HRD)','SF3B1-BUD31 (non-HRD)','PRPF8-USP39 > RD)','PRPF8-USP39 (non-HRD)'),cex=2.5, x=c(1:10), y=-4.5, srt=45, xpd=T, adj=1)
dev.off()

pdf('RNPS1_RPL8.pdf', width=10, height=5)
par(mfrow=c(1,2), las=1)
plot(quant_data_all['Q15287',1:58], xlab='', ylab='', type='l', lwd='3', col='palevioletred1', ylim=c(8,22))
lines(quant_data_all['P62917',1:58], type='l', lwd='3', col='lightgoldenrod')
plot(quant_data_all['Q15287',59:103], xlab='', ylab='', type='l', lwd='3', col='palevioletred1', ylim=c(8, 22))
lines(quant_data_all['P62917',59:103], type='l', lwd='3', col='lightgoldenrod')
dev.off()


pdf('RNPS1_SNRNP200.pdf', width=10, height=5)
par(mfrow=c(1,2), las=1)
plot(quant_data_all['Q15287',1:58], xlab='', ylab='', type='l', lwd='3', col='palevioletred1', ylim=c(8,22))
lines(quant_data_all['O75643',1:58], type='l', lwd='3', col='lightgoldenrod')
plot(quant_data_all['Q15287',59:103], xlab='', ylab='', type='l', lwd='3', col='palevioletred1', ylim=c(8, 22))
lines(quant_data_all['O75643',59:103], type='l', lwd='3', col='lightgoldenrod')
dev.off()

pdf('RNPS1_SF3B1.pdf', width=10, height=5)
par(mfrow=c(1,2), las=1)
plot(quant_data_all['Q15287',1:58], xlab='', ylab='', type='l', lwd='3', col='palevioletred1', ylim=c(8,22))
lines(quant_data_all['O75533',1:58], type='l', lwd='3', col='lightgoldenrod')
plot(quant_data_all['Q15287',59:103], xlab='', ylab='', type='l', lwd='3', col='palevioletred1', ylim=c(8, 22))
lines(quant_data_all['O75533',59:103], type='l', lwd='3', col='lightgoldenrod')
dev.off()

pdf('SF3B1_BUD31.pdf', width=10, height=5)
par(mfrow=c(1,2), las=1)
plot(quant_data_all['O75533',1:58], xlab='', ylab='', type='l', lwd='3', col='palevioletred1', ylim=c(8,22))
lines(quant_data_all['P41223',1:58], type='l', lwd='3', col='lightgoldenrod')
plot(quant_data_all['O75533',59:103], xlab='', ylab='', type='l', lwd='3', col='palevioletred1', ylim=c(8, 22))
lines(quant_data_all['P41223',59:103], type='l', lwd='3', col='lightgoldenrod')
dev.off()

pdf('PRPF8_USP39.pdf', width=10, height=5)
par(mfrow=c(1,2), las=1)
plot(quant_data_all['Q6P2Q9',1:58], xlab='', ylab='', type='l', lwd='3', col='palevioletred1', ylim=c(8,22))
lines(quant_data_all['Q53GS9',1:58], type='l', lwd='3', col='lightgoldenrod')
plot(quant_data_all['Q6P2Q9',59:103], xlab='', ylab='', type='l', lwd='3', col='palevioletred1', ylim=c(8, 22))
lines(quant_data_all['Q53GS9',59:103], type='l', lwd='3', col='lightgoldenrod')
dev.off()


## picked samples vs. overall reference distribution



# compare picked ratios to reference distribution
a <- tmp1[ c('HRD_1','HRD_102','HRD_103','HRD_29','HRD_34','HRD_43','HRD_5','HRD_59','HRD_60','HRD_61','HRD_65','HRD_68','HRD_77','HRD_83','HRD_85','HRD_87','HRD_9')]

pdf("ratio_distributions_RNPS1_RPL8.pdf", width=10, height=10)
par(cex.axis=2, las=1)
plot(density(tmp1[59:103]), lwd='3', col='grey', xlab='', ylab='', main='', ylim=c(-0.1, 0.5))
#lines(density(tmp1[1:58]), lwd='3', col='green', xlab='', ylab='', main='')
points(x=a, y=jitter(rep(0,length(a)),1),lwd='3', col='red')
dev.off()

##
a <- tmp2[ c('HRD_1','HRD_102','HRD_103','HRD_29','HRD_34','HRD_43','HRD_5','HRD_59','HRD_60','HRD_61','HRD_65','HRD_77','HRD_83','HRD_85','HRD_87','HRD_9')]

pdf("ratio_distributions_RNPS1_SNRNP200.pdf", width=10, height=10)
par(cex.axis=2, las=1)
plot(density(tmp2[59:103]), lwd='3', col='grey', xlab='', ylab='', main='', ylim=c(-0.1, 0.5))
points(x=a, y=jitter(rep(0,length(a)),1),lwd='3', col='red')
dev.off()


##
a <- tmp3[ c('HRD_1','HRD_102','HRD_103','HRD_29','HRD_34','HRD_43','HRD_5','HRD_59','HRD_60','HRD_61','HRD_65','HRD_77','HRD_83','HRD_85','HRD_87','HRD_9')]

pdf("ratio_distributions_RNPS1_SF3B1.pdf", width=10, height=10)
par(cex.axis=2, las=1)
plot(density(tmp3[59:103]), lwd='3', col='grey', xlab='', ylab='', main='', ylim=c(-0.1, 0.5))
points(x=a, y=jitter(rep(0,length(a)),1),lwd='3', col='red')
dev.off()

##
a <- tmp4[c('HRD_1','HRD_102','HRD_103','HRD_23','HRD_29','HRD_34','HRD_43','HRD_44','HRD_5','HRD_61','HRD_63','HRD_68','HRD_69','HRD_70','HRD_72','HRD_82','HRD_83','HRD_84','HRD_87','HRD_9')]

pdf("ratio_distributions_SF3B1_BUD31.pdf", width=10, height=10)
par(cex.axis=2, las=1)
plot(density(tmp4[59:103]), lwd='3', col='grey', xlab='', ylab='', main='', ylim=c(-0.1, 0.5))
points(x=a, y=jitter(rep(0,length(a)),1),lwd='3', col='red')
dev.off()

##
a <- tmp5[ c('HRD_1','HRD_102','HRD_103','HRD_14','HRD_23','HRD_29','HRD_34','HRD_36','HRD_59','HRD_60','HRD_65','HRD_71','HRD_75','HRD_85','HRD_93','HRD_98')]

pdf("ratio_distributions_PRPF8_USP39.pdf", width=10, height=10)
par(cex.axis=2, las=1)
plot(density(tmp5[59:103]), lwd='3', col='grey', xlab='', ylab='', main='', ylim=c(-0.1, 0.5))
points(x=a, y=jitter(rep(0,length(a)),1),lwd='3', col='red')
dev.off()




