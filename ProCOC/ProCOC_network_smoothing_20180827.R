## network_smoothing on partially filtered network


setwd('Documents/Biomarker_Cancer/ProCOC/ProCOC_final/ProCOC1_correlation_analysis/')
load('ProCOC1_correlation_analysis_20180827.Rdata')
ls()
[1] "a"                     "mean_distribution"     "names_prot1"          
[4] "ProCOC_log2FC"         "ProCOC_log2FC_reduced" "prot_gleason" 

network <- read.table('/home/betty/Downloads/9606.protein.links.detailed.v10.5.txt.gz', quote='\"', comment.char="", stringsAsFactors=F, header=T)
network[,1] <- as.character(substr(network[,1],6,20))
network[,2] <- as.character(substr(network[,2],6,20))
head(network[,1:2])
         protein1        protein2
1 ENSP00000000233 ENSP00000263431
2 ENSP00000000233 ENSP00000353863
3 ENSP00000000233 ENSP00000342026
4 ENSP00000000233 ENSP00000240874
5 ENSP00000000233 ENSP00000379847
6 ENSP00000000233 ENSP00000400088

dim(network)
[1] 11353056       10

network <- network[network$combined_score>=800,]
dim(network)
[1] 632286     10
head(a)
     O00391      O00533      O14786      O75636      O75882      O95445 
 0.10157997  0.22653287  0.14357709 -0.31863613  0.25414949  0.03160904 
names_prot <- names(a)

mapping_table <- read.delim('M2018', sep='\t', dec='.', header=T, as.is=T)
head(mapping_table)

mapping_table <- mapping_table[!duplicated(mapping_table$From),]
dim(mapping_table)
[1] 193   2 

a <- a[grep('P01860',names(a), invert=T)]
length(a)
[1] 192
a <- a[grep('sp|Q3SZR3|A1AG_BOVIN',names(a), invert=T)]
a <- a[grep('sp|Q58D62|FETUB_BOVIN',names(a), invert=T)]
length(a)
[1] 191

mapping_table <- mapping_table[1:191,]
dim(mapping_table)
[1] 191   2

expr_mat <- data.frame(as.character(mapping_table[,2]), as.numeric(a), as.numeric(a))
dim(expr_mat)
[1] 191   3

head(expr_mat)
  as.character.mapping_table...2.. as.numeric.a. as.numeric.a..1
1                  ENSP00000356572    0.10157997      0.10157997
2                  ENSP00000256509    0.22653287      0.22653287
3                  ENSP00000265371    0.14357709      0.14357709
4                  ENSP00000270879   -0.31863613     -0.31863613
5                  ENSP00000262919    0.25414949      0.25414949
6                  ENSP00000365081    0.03160904      0.03160904


colnames(expr_mat)<- c('gene_id', 'correlation_1', 'correlation_2')
library(BioNetSmooth)
Loading required package: igraph

Attaching package: ‘igraph’

The following objects are masked from ‘package:stats’:

    decompose, spectrum

The following object is masked from ‘package:base’:

    union

Loading required package: Matrix

Attaching package: ‘BioNetSmooth’

The following objects are masked _by_ ‘.GlobalEnv’:

    expr_mat, network

Warning message:
package ‘BioNetSmooth’ was built under R version 3.3.0 
>

length(intersect(network[,1],expr_mat[,1]))
[1] 128
netmap <- network_mapping(network, expr_mat, type = "laplacian", merge.by = "gene_id", global = TRUE)
Network with: 13273 nodes

netmap$smoothmat <- network_smoothing(net = netmap$G, mat_intensities = netmap$mat_intensities, conditions = colnames (netmap$mat_intensities), iter = 20, alpha = 0.5, network_type = "laplacian")
summary(netmap)
                Length    Class     Mode     
G               176172529 dgCMatrix S4       
mat_intensities     26546 -none-    numeric  
gene_names          13273 -none-    character
smoothmat           26546 -none-    numeric  


b <- netmap$smoothmat[,1]
names(b)<- netmap$gene_names
d<- b[order(b)]
length(b)
[1] 13273
x <- head(d, 133)
d <-b[order(b, decreasing=T)]
y <- head(d, 133)

write.table(names(y), file='positive_prots_20180827.txt', quote=F, row.names=F)
write.table(names(x), file='negative_prots_20180827.txt', quote=F, row.names=F)
positive_prots <- read.delim('uniprot-', sep='\t', dec='.', header=T, as.is=T)

write.table(positive_prots[,2], file='positive_prots_20180827.txt', quote=F, row.names=F)

negative_prots <- read.delim('uniprot-yourlist%3AM20180827A7434721E10EE6586998A056CCD0537EA5A24DZ.tab', sep='\t', dec='.', header=T, as.is=T)

write.table(negative_prots[,2], file='negative_prots_20180827.txt', quote=F, row.names=F)
