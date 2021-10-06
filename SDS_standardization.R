# script to standardize raw SDS values

# Read raw SDS file
SDS <- read.table("chr22_sds_test.txt", header=TRUE)

# Bin the SDS dataframe into bins based on 1% DAF intervals
bin_index <- findInterval(SDS$DAF, seq(.05, .95, .01) )

# check if there is any bin with only one value 
# sum(table(bin_index)==1)>0
# FALSE

# initialize standardized SDS column
SDS$sSDS <- rep( NA, dim(SDS)[1] )

# perform the standardization
# daf_interval becomes a different 1% DAF interval on each iteration
for( daf_interval in unique(sort(bin_index)) ){
  interval_sds <- SDS[ bin_index==daf_interval , "rSDS"] # extract the sds for that bin interval
  standardize_sds <- (interval_sds - mean(interval_sds)) / sd(interval_sds)
  SDS$sSDS[bin_index==daf_interval] <- standardize_sds
}

# "P values are two-sided tail probabilities of standard normal." - Field, et al. 2016 
# "Two-tailed p-values were converted by whole genome-wide standardized SDS z-scores." - Peikuan Cong, et al. (Westlake)
# to calculate 2-sided pvals the formula is pvalue2sided = 2*pnorm(-abs(z)); where z are the standardize Z-scores
SDS$pval <- 2*pnorm( -abs(SDS$sSDS) )

# plot
library(wordcloud)

plot(SDS$POS, -log10(SDS$pval),pch=16, ylab="-log10 Pval", xlab="chr22",las=2, col="grey")
abline(h=60, lty = 2, col="grey") # threshold
points(SDS$POS[-log10(SDS$pval)>60], -log10(SDS$pval)[-log10(SDS$pval)>60],pch=16) # snps that pass thr on a diff color
pos_text <- SDS$ID[-log10(SDS$pval)>60] # extract text
# define wordlayout
nc <- wordlayout( SDS$POS[-log10(SDS$pval)>60], -log10(SDS$pval)[-log10(SDS$pval)>60], words=pos_text)
nc[,1] <- nc[,1]-1000000
text(nc[,1],nc[,2], label=pos_text, cex=1)


