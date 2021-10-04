library(utils)
library(stats)

#-----------------------------------------------------------
# Global constants/parameters
#-----------------------------------------------------------

PRECISION      = 4
E_GRID_NUM_POINTS    = 50
E_GRID_SCALE_FACTOR  = 20
OPTIM_NUM_ITERATIONS    = 5

# when test-SNPs are close to a boundary (e.g., end-of-chromosome or centromere) individuals may don't have any signleton
# between the test-SNP and the nearby boundary. In such cases we truncate the observations (see below).
# However, if this happens to more than 5% of the individuals we don't make any prediction and skip over the test-SNP.

skip_boundary_missing_singletons_fraction_threshold = 0.05


#-----------------------------------------------------------
# Process the arguments, load the files, etc.
#-----------------------------------------------------------

args <- commandArgs(TRUE)

if (length(args)==1 && (args[1]=="-h" | args[1]=="--h" | args[1]=="-help" | args[1]=="--help")) { 
   write(HELP, stdout())
   quit()
}

if (length(args) < 6) {
   write("\nExpecting more arguments... see details with the -h flag.\n", stdout())
   quit()
}

# s_file

if (args[1]=="-") {
    sin_snp_file_name <- "stdin"
} else {
    sin_snp_file_name <- args[1]
}

sin_snp_fh <- file(sin_snp_file_name, 'r')

# t_file

if (args[2]=="-") {
    test_snp_file_name <- "stdin"
} else {
    test_snp_file_name <- args[2]
}

test_snp_fh <- file(test_snp_file_name, 'r')

# o_file

if (args[3]=="-") {
    sin_observability_file_name <- "stdin"
} else {
    sin_observability_file_name <- args[3]
}

# b_file

if (args[4]=="-") {
    boundaries_file_name <- "stdin"
} else {
    boundaries_file_name <- args[4]
}

# g_file

if (args[5]=="-") {
    gamma_shape_file_name <- "stdin"
} else {
    gamma_shape_file_name <- args[5]
}

# "initial guess" param is a starting point for optimization,
# and the center of a grid from which other starting points are
# randomly selected

E_grid_center <- as.numeric(args[6])

# an upper bound on the number of singletons per individual
# -- so that R's read.table works to read the singletons file allowing different number of singletons per individual
# -- default to 10000, which is a fairly conservative bound. but if you know your input files, use a tighter bound

MAX_SINGLETONS_PER_INDV = 10000

if (length(args)>6) {
   MAX_SINGLETONS_PER_INDV <- as.numeric(args[7])
}



# Read the singletons file and store in memory

singletons_df<-read.table(sin_snp_fh, sep="", header=FALSE, row.names=NULL, col.names = paste0("V",seq_len(MAX_SINGLETONS_PER_INDV)), fill=TRUE)
non_na_columns <- sum(colSums(singletons_df,na.rm=TRUE)>0)
singletons_df <- singletons_df[,-((non_na_columns+1):MAX_SINGLETONS_PER_INDV)]
singletons_nrow=nrow(singletons_df)
singletons<-as.matrix(singletons_df, nrow=singletons_nrow)
rm(singletons_df)
close(sin_snp_fh)

singletons_current_indices <- 1:singletons_nrow*0+1

sin_observability_fh <- file(sin_observability_file_name, 'r')
line<-unlist(strsplit(readLines(sin_observability_fh,n=1),"\t|\ "))

if (length(line) > 0) {
   sin_observability <- as.numeric(c(line[1:length(line)]))
} else {
   sin_observability <- 1:singletons_nrow*0+1
}

# scaling sin_observability does not affect the reported statistic.
# however, the initial point is sensitive to this scaling, so we always
# scale to mean observability one, as in the full observability case
sin_observability = sin_observability/mean(sin_observability)

# Read the boundaries file

boundaries_fh <- file(boundaries_file_name, 'r')
boundaries_df <- read.table(boundaries_fh, sep="", header=FALSE, row.names=NULL)
boundaries_nrow = nrow(boundaries_df)
boundaries <- as.matrix(boundaries_df, nrow=boundaries_nrow)
rm(boundaries_df)
close(boundaries_fh)
boundaries_current_index <- 1

# Read the gamma-shape parameters file

gamma_shape_fh <- file(gamma_shape_file_name, 'r')
gamma_shape_df <- read.table(gamma_shape_fh, sep="", header=FALSE, row.names=NULL)
gamma_shape_nrow = nrow(gamma_shape_df)
gamma_shape <- as.matrix(gamma_shape_df, nrow=gamma_shape_nrow)
gamma_shape_freq  <- gamma_shape[,1]
gamma_shape_shape <- gamma_shape[,2]
rm(gamma_shape_df)
close(gamma_shape_fh)

#-----------------------------------------------------------
# Internal functions
#-----------------------------------------------------------

# interpolate the gamma shape paramter for a given allele frequency (using the gamma shape input file)

# note: in theory the shape parameter should not depend only on frequency but also on
# whether the allele is ancestral or derived. however the differences are small
# and these shape parameters are estimated from simple demographic models that are
# almost surely mispecified... which makes these further small differences negligible.
# any such potential frequency-dependent biases should be eliminated when the raw SDS scores
# are standardized by derived allele frequency bins

get_gamma_shape <- function(freq) {
   index <- findInterval(freq,gamma_shape_freq)
   if (index == 0) {
      index = 1
   }
   if (index == length(gamma_shape_freq)) {
      index = length(gamma_shape_freq) - 1
   }
   #freq is between index & index+1
   y1 = gamma_shape_shape[index]
   y2 = gamma_shape_shape[index+1]
   x1 = gamma_shape_freq[index]
   x2 = gamma_shape_freq[index+1]
   shape_interpolated = y1 + (y2-y1)*(freq-x1)/(x2-x1)
   return(shape_interpolated)
}

# Compute log(a+b) from log(a) and log(b) 

#logsum <- function(log_a,log_b) {
#  # why not use log_a (-Inf) == -Inf
#  if (log_a <= -Inf | log_b <= -Inf) {
#     if (log_a <= -Inf) {res = log_b}
#     else {res = log_a}
#  }
#  else {
#     if (log_a > log_b) {res = log_a + log(1 + exp(log_b - log_a))}
#     else {res = log_b + log(1 + exp(log_a - log_b))}
#  }
#  return(res)
#}

logsum <- function(log_a,log_b) {
   # log_a will be a vector
   # log_b is a single value

   # init res as the vector log_a
   res <- log_a

   # if (log_a <= -Inf):
   res[which(log_a <= -Inf)] <- log_b[which(log_a <= -Inf)]

   # else
   if( !(log_b <= -Inf) ){
        res[ which( log_a > log_b ) ] <- (log_a + log(1 + exp(log_b - log_a)))[which( log_a > log_b )]
     res[ which( log_a <= log_b ) ] <- (log_b + log(1 + exp(log_a - log_b)))[which( log_a <= log_b )]
   }
   return(res)
}

# Log-Likelihood (returns the *minus* log-likelihood because 'nlm' is *minimizing* the objective function)

f_minus_log_likelihood <- function(x) {
  res_LL <- f_minus_log_likelihood0( c(x[1],x[2]) )
  return(res_LL)
}

f_minus_log_likelihood0 <- function(x) {
  res_LL <- f_minus_log_likelihood1(x)
  if (is.nan(res_LL)) {
     #return(-Inf)
     return(-1e15)
  }
  return(res_LL)
}

f_minus_log_likelihood1 <- function(x) {

  # E1,E2 are the model parameters: the mean tip lengths of the ancestral and derived alleles (in mutation rate units)
  logE1 = x[1]
  logE2 = x[2]

  # A1,A2 - the gamma shape parameters
  # B1,B2 - the gamma rate parameters
  logB1 = logA1 - logE1
  logB2 = logA2 - logE2

  # LL0,LL1,LL2 - log-likelihod components for the three genotype groups 
  LL0 = LL1 = LL2 = 0

  # n0,n1,n2 - number of individuals in each genotype group
  if (n0>0) {
     LL0 = 2.0*A1*(logB1 - mean(logsum(log(dat0),logB1))) + mean(log(dat0)) + log(2) + logA1 - 2.0*mean(logsum(log(dat0),logB1)) + mean(logsum(log(2)+logA1,0))
  }
  if (n2>0) {
     LL2 = 2.0*A2*(logB2 - mean(logsum(log(dat2),logB2))) + mean(log(dat2)) + log(2) + logA2 - 2.0*mean(logsum(log(dat2),logB2)) + mean(logsum(log(2)+logA2,0))
  }
  if (n1>0) {
     LL1 = A1*(logB1 - mean(logsum(log(dat1),logB1))) + A2*(logB2 - mean(logsum(log(dat1),logB2))) + mean(log(dat1)) + mean(logsum(logsum(-2.0*logsum(log(dat1),logB1)+logA1+logsum(logA1,0), -2.0*logsum(log(dat1),logB2)+logA2+logsum(logA2,0)), log(2)+logA1+logA2-logsum(log(dat1),logB1)-logsum(log(dat1),logB2)))
  }
  LL = LL0*n0 + LL1*n1 + LL2*n2

  return(-LL)
}


#-----------------------------------------------------------
# Main loop over test-SNPs
#-----------------------------------------------------------

# Write header line

my_header = paste( "ID", "AA", "DA", "POS", "DAF", "nG0", "nG1", "nG2", "rSDS", "SuggestedInitPoint", sep="\t" )
write(my_header,stdout())
#write(my_header,"chr22_sds_test.txt")

# Read the test SNPs and process them one by one (altough we read them in chunks of 10000 snps)

options(warn=-1)

while(length(lines<-readLines(test_snp_fh, n<-10000))>0) { 
    for (i in 1:length(lines)) {

   line<-unlist(strsplit(lines[i],"\\s+"))

   test_snp_id <- line[1]

   #allele1 will correspond to the ancestral allele (allele2 for the derived)

   test_snp_allele1 <- line[2]
   test_snp_allele2 <- line[3]
   test_snp_location <- as.numeric(line[4])
   test_snp_genotypes <- as.numeric(c(line[5:length(line)]))
   test_snp_Nminus1 = length(test_snp_genotypes)-1

   tmp_singletons_upstream_intervals <- 1:length(test_snp_genotypes)*NA
   tmp_singletons_downstream_intervals <- 1:length(test_snp_genotypes)*NA

   # find the boundaries for the current test-SNP
    # first condition: have we looked beyond the final chromosome boundary
    # second condition: is the snp coordinate higher than the end of the current boundary
    # if both true, check the next boundary
   while (boundaries_current_index <= boundaries_nrow && boundaries[boundaries_current_index,2] < test_snp_location) {
      boundaries_current_index = boundaries_current_index+1
   }

   # first condition: have we looked beyond the final chromosome boundary
   # second condition: is the snp coordinate upstream of the start of the boundary
   if (boundaries_current_index <= boundaries_nrow && boundaries[boundaries_current_index,1] <= test_snp_location) {

      # compute for each individual the distance to the nearest upstream singleton and the distance to the nearest downstream singleton

      # save the current boundaries to a tmp variable
      tmp_boundary_upstream = boundaries[boundaries_current_index,1]
      tmp_boundary_downstream = boundaries[boundaries_current_index,2]

      # for each sample (rows); ind_i keeps track of the sample
      for (ind_i in 1:singletons_nrow) {
           # while loop on every singleton (columns of the singleton matrix)
           # first condition: checks if we have not looped through all the singletons
           # second condition: checks that the singleton we are currently checking is not NA
         while(singletons_current_indices[ind_i]<=dim(singletons)[2] && !is.na(singletons[ind_i,singletons_current_indices[ind_i]]) ) {
             # if the current singleton coordinate is higher than the test snp, breaks the while
            if (as.numeric(singletons[ind_i,singletons_current_indices[ind_i]]) >= test_snp_location) {
          break
          }
           # else (if the current singleton coordinate is not higher or equal thatn the test snp)
           # for the current sample, change to the next singleton 
            singletons_current_indices[ind_i] = singletons_current_indices[ind_i]+1
           }
         # if the current singleton index in the sample ind_i, is not the first singleton
         if (singletons_current_indices[ind_i] >= 2) {
             # save the previos singleton in the tmp variable
             tmp_singleton_location = as.numeric(singletons[ind_i,singletons_current_indices[ind_i]-1])
             # check if the previous singleton is after of the boundary start; get the left distance of the singleton to the test snp
             if (tmp_singleton_location >= tmp_boundary_upstream) {
               tmp_singletons_upstream_intervals[ind_i] = test_snp_location - tmp_singleton_location
             }
         }
          # first condition: the singleton current index is not out of bounds
          # second condition: the current singleton is not NA
         if (singletons_current_indices[ind_i]<=dim(singletons)[2] && !is.na(singletons[ind_i,singletons_current_indices[ind_i]])) {
             # if the current singleton coordinate is after or equal the coordinate of the test snp
             if (as.numeric(singletons[ind_i,singletons_current_indices[ind_i]]) >= test_snp_location) {
            tmp_singleton_location = as.numeric(singletons[ind_i,singletons_current_indices[ind_i]])
           # if the singleton coordinate is before the end of the boundary; calulated the right distance between the singleton and the test snp
           if (tmp_singleton_location <= tmp_boundary_downstream) {
             tmp_singletons_downstream_intervals[ind_i] = tmp_singleton_location - test_snp_location
           }
          }
         }
     }

      # For snps close to boundary we are still left with NA values for nearest singletons.
      # If more than 100*skip_boundary_missing_singletons_fraction_threshold % of the individuals are NA we skip this snp;
      # Otherwise, we set the NA entries to the largest distance observed.

      # tmp_singletons_upstream_intervals contains the left distance from the closest singletons (per sample) to test snp
      # tmp_singletons_downstream_intervals contains the right distance from the closest singletons (per sample) to test snp
      # this if checks that there is singleton information for 95% of the samples to the left and right of the test snp
      if (mean(is.na(tmp_singletons_upstream_intervals)) <= skip_boundary_missing_singletons_fraction_threshold && mean(is.na(tmp_singletons_downstream_intervals)) <= skip_boundary_missing_singletons_fraction_threshold) {
         #  print("passed the 95% missing singleton threshold")
         # replace NAs with the maximum distance value
         tmp_singletons_upstream_intervals[is.na(tmp_singletons_upstream_intervals)] = max(tmp_singletons_upstream_intervals,na.rm=TRUE)
         tmp_singletons_downstream_intervals[is.na(tmp_singletons_downstream_intervals)] = max(tmp_singletons_downstream_intervals,na.rm=TRUE)

         tmp_singleton_intervals = tmp_singletons_upstream_intervals + tmp_singletons_downstream_intervals

         # Correct for singleton "observability" (e.g., the prob. that a singleton is detected due to low sequencing depth)

         tmp_singleton_intervals = tmp_singleton_intervals * sin_observability
         # there was a bug here: if there are NA in the genotypes list the AF is NA
         test_snp_af = mean(test_snp_genotypes, na.rm=TRUE)/2

         #
         # Compute rSDS
         #

         dat0 <- tmp_singleton_intervals[which(test_snp_genotypes==0)]
         dat1 <- tmp_singleton_intervals[which(test_snp_genotypes==1)]
         dat2 <- tmp_singleton_intervals[which(test_snp_genotypes==2)]
         dat012 <- tmp_singleton_intervals
         n0 = length(dat0)
         n1 = length(dat1)
         n2 = length(dat2)


         # gamma shape parameters by the allele frequency (the ancestral/derived annotation is ignored)

         A1 = get_gamma_shape(1-test_snp_af)
         A2 = get_gamma_shape(test_snp_af)
         logA1 = log(A1)
         logA2 = log(A2)

         # a grid of possible initial points around the input guess
 
         logE = seq(from=log(E_grid_center)-log(E_GRID_SCALE_FACTOR),to=log(E_grid_center)+log(E_GRID_SCALE_FACTOR),length=E_GRID_NUM_POINTS)

         best_MLE_nlm_iteration  = 0
         best_MLE_log_likelihood = -Inf
         best_MLE_param_estimate = -Inf
         best_MLE_nlm_iteration  = 0

         # we make run the optimization from several random starting points on the grid,
         # as well as from the grid center which is the input guess

         iter = 1
         while (iter <= OPTIM_NUM_ITERATIONS) {
            tmp_init_point = c(sample(logE,1),sample(logE,1))
             tmp_nlm_fit <- nlm(f_minus_log_likelihood, tmp_init_point)
             tmp_MLE_est <- tmp_nlm_fit$estimate
             tmp_MLE_ll <- -1.0*tmp_nlm_fit$minimum
             if (tmp_MLE_ll > best_MLE_log_likelihood) {
                best_MLE_log_likelihood = tmp_MLE_ll
                best_MLE_param_estimate = tmp_MLE_est
                best_MLE_nlm_iteration = iter
            }
             iter = iter+1
         }

         # test the E_grid_center as starting point

         tmp_init_point = c(log(E_grid_center),log(E_grid_center))
         tmp_nlm_fit <- nlm(f_minus_log_likelihood, tmp_init_point)
         tmp_MLE_est <- tmp_nlm_fit$estimate
         tmp_MLE_ll <- -1.0*tmp_nlm_fit$minimum
         if (tmp_MLE_ll > best_MLE_log_likelihood) {
               best_MLE_log_likelihood = tmp_MLE_ll
                best_MLE_param_estimate = tmp_MLE_est
                best_MLE_nlm_iteration = iter
         }


         # that's it - write the results for the best (over the different initial guesses) maximum likelihood estimate

         # raw SDS is the log ratio of inferred mean tip length
         # (the mean tip lengths here were parameterized here in log space, so rSDS is the difference betweeen the estimated parameters) 

         test_snp_rSDS <- best_MLE_param_estimate[1] - best_MLE_param_estimate[2]

         suggested_init_point_param = paste("1e",format(mean(best_MLE_param_estimate)/log(10),digits=1),sep="")
         my_str = paste( test_snp_id, test_snp_allele1, test_snp_allele2, test_snp_location,
                     format(test_snp_af,digits=PRECISION),
               n0,
               n1,
               n2,
                     format(test_snp_rSDS,digits=PRECISION),
               suggested_init_point_param,
               sep="\t" )
         write(my_str,stdout())
         #write(my_str,"chr22_sds_test.txt")
         #print(my_str)
           }
   }
    }
}
