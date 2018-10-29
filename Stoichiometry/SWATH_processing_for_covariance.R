## Data processing for Stoichio workflow  ## 
## Less missing values -> less imputations ##

## Step 1: Prepare R session ##

setwd('Documents/Biomarker_Cancer/CPTAC_Covariance/')

## Step 2: Load in data ##

data <- read.delim('/home/betty/Documents/Biomarker_Cancer/Ovarian_JohnHopkins/2017_10_Paper/CPTAC_103_wocomp_newlib_jumbo_min2max3_quantnorm_batchcorr/fragments_for_protein_quantification.txt', sep='\t', dec='.', as.is=T, header=T)
dim(data)
[1] 150408    106

## Step 3: Filter data and select top3 fragments

data_filtered <- data[1:90421,]
dim(data_filtered)
[1] 90421   106
tail(data_filtered[,1])
[1] "Q9Y6Y8" "Q9Y6Y8" "Q9Y6Y8" "Q9Y6Y8" "Q9Y6Y8" "Q9Y6Y8"

data_filtered$fragment_sum <- apply(as.matrix(data_filtered[,4:106]), 1, sum, na.rm=T)
 
l <- data_filtered[order(data_filtered$fragment_sum, decreasing=T),]
d <- by(l, l$Peptide, head, n=3)
data_top3 <- Reduce(rbind, d)

## Step 4: Aggregate fragments to peptides by summing the top3 fragments

length(unique(data_top3$Peptide))
[1] 15684

data_peptide <- aggregate(data_top3[,4:106], by=list(data_top3$Protein, data_top3$Peptide), FUN=sum)

dim(data_peptide)
[1] 15684   105

length(unique(data_peptide$Group.1))
[1] 2914

rownames(data_protein) <- data_protein[,1]
data_protein[,1]<- NULL
data_protein <- as.matrix(data_protein)
dim(data_protein)
[1] 2914  103

data_log <- log2(data_protein)
data_log[which(is.infinite(data_log))] <- NA
dim(data_log)
[1] 2914  103


## Step 5: Check protein abundance and missing values

mean_abund <- apply(data_protein, 1, mean, na.rm=T)
length(mean_abund)
[1] 2914
head(mean_abund)
  A0AV96   A0AVT1   A0FGR8   A0MZ66   A1L0T0   A2RRP1 
8779.571 7296.150 8139.078 2931.655 5623.489 2063.641 
mean_abund <- sort(mean_abund)
length(mean_abund)
[1] 2910

head(mean_abund)
  O60443   O60287   P21757   Q03013   Q5T2E6   Q96HW7 
1182.027 1276.589 1288.127 1292.845 1313.430 1330.859 


data_protein_sorted <- data_log[names(mean_abund),]
data_protein_sorted[1:3, 1:3]

SWATH_missing <- NULL
for (i in 1:nrow(data_protein_sorted)) {
tmp <- length(which(is.nan(as.matrix(data_protein_sorted[i,]))))
SWATH_missing <- c(SWATH_missing, tmp)
}

length(SWATH_missing)
[1] 2910
names(SWATH_missing) <- rownames(data_protein_sorted)
head(SWATH_missing)
O60443 O60287 P21757 Q03013 Q5T2E6 Q96HW7 
    87     90     89     94    100     68 


quantile(mean_abund, c(0, 0.1, 0.2, 0.3, 0.4, 0.5))
      0%      10%      20%      30%      40%      50% 
1182.027 3861.542 5191.340 6165.489 7265.060 8630.774 


        0%        10%        20%        30%        40%        50%        60% 
  1182.027   3861.542   5191.340   6165.489   7265.060   8630.774  10199.196 
       70%        80%        90%       100% 
 12494.333  15950.821  24879.027 828064.649 

length(which(mean_abund<=5191.340))
[1] 582

length(which(mean_abund<=7265.060))
[1] 1164

length(which(mean_abund<=10199.196))
[1] 1746

length(which(mean_abund<=15950.821))
[1] 2328




a <- names(which(SWATH_missing[1:582]<11))
b <- names(which(SWATH_missing[583:1164]<21))
c <- names(which(SWATH_missing[1165:1746]<31))
d <- names(which(SWATH_missing[1747:2328]<42))
e <- names(which(SWATH_missing[2329:2910]<52))


SWATH_filtered <- data_protein_sorted[c(a,b,c,d,e),]
dim(SWATH_filtered)
[1] 1475  103
length(which(is.na(SWATH_filtered)))
13487

13487/(103*1475)
[1] 0.08877407

length(which(is.na(data_protein)))
[1] 101689
> dim(data_protein)
[1] 2914  103

101689/(103*2914)
[1] 0.338803

from 34% missing values down to 9% missing values

SWATH_filtered[1:4, 1:4]
count_missing <- length(which(is.na(SWATH_filtered)))

new_mean <- mean(SWATH_filtered, na.rm=T)*0.8
new_sd <- sd(SWATH_filtered, na.rm=T)*0.4
set.seed(0)
background <- rnorm(count_missing, mean=new_mean, sd=new_sd)

#plot in order to check where random distribution lays
data_hist <- hist(SWATH_filtered)
back_hist <- hist(background)
plot(data_hist, col='lightblue', main='Distribution of data and distribution of background', xlab='log2 intensities')
plot(back_hist, col='red', main='Distribution of data and distribution of background', xlab='log2 intensities', add=T)
SWATH_imputed <- SWATH_filtered
length(which(is.na(SWATH_imputed)))
[1] 21914
SWATH_imputed[which(is.na(SWATH_imputed))]<-background
length(which(is.na(SWATH_imputed)))
[1] 0

save(SWATH_imputed, SWATH_filtered, data_protein, data_peptide, data_top3, data_log, data_filtered, file='CPTAC_covariance_SWATH_processing_20180912.Rdata')

quant_data_all <- SWATH_imputed

save(quant_data_all, file='Stoichio_input_20180912.Rdata')

