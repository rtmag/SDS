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

