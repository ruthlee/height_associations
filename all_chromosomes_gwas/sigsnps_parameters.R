setwd("~/R.Projects/height_associations/all_chromosomes_gwas")

raw.snps <- read.table("all_snps_GENO")

colnames(raw.snps) <- c("#Chr", "Position", "ID", "Ref", "Alt", "Test", "Obs_ct", "Beta", "Se", "T-stat", "P")

# we want to only keep rows with significant pvalues in raw.snps
sigsnps <- raw.snps[ raw.snps$P <= 10^-5 & !is.na(raw.snps$P), ]

# reading in joe's blocks file
joes.blocks <- read.table ( "fourier_ls-all.bed", stringsAsFactors = FALSE , header = TRUE )

# chromosome ID to match domdev frame chromosome ID's
joes.blocks$chr <- gsub ( 'chr' , '' , joes.blocks$chr )

# row ID creation
joes.blocks <- cbind (joes.blocks, "I" = 1:nrow (joes.blocks) )

  # progress bar for assigning joe's blocks numbers to each row in significant snps frame (longest part of program)
pb <- txtProgressBar( min = 0, max = nrow ( sigsnps ) , style= 3 )
for ( i in 1:nrow ( sigsnps ) ) {
  for ( j in 1:nrow ( joes.blocks ) ) {
    if ( sigsnps [ i, "Position" ] >= joes.blocks [ j, "start" ] &
         sigsnps [ i, "Position" ] <= joes.blocks [ j, "stop" ] &
         sigsnps [ i, "#Chr" ] == joes.blocks [ j, "chr" ] ) {
      sigsnps [ i, "Block" ] <- joes.blocks [ j, "ID" ]
    }
  }
  setTxtProgressBar( pb, i )
}

# this took forever, but the csv sigsnps_blocks is the result

write.table(sigsnps, "sigsnps_blocks.csv", sep=",")

sigsnps <- read.table("sigsnps_blocks.csv", sep=",", stringsAsFactors = FALSE )

colnames(sigsnps) <- c("#Chr", "Position", "ID", "Ref", "Alt", "Test", "Obs_ct", "Beta", "Se", "T-stat", "P", "Block")

# Pulling the most significant snps out of each block according to the dominance deviation or average effect pvalue. The function returns warnings for some blocks but will still output a dataframe of significant snps.
most_sig_snps <- function(dataframe, avgeff = FALSE, domdev = FALSE) {
    if (avgeff) {
        blocks <- unique(dataframe$Block)
        sigsnps_df <- data.frame(matrix(ncol = ncol(dataframe)))
        colnames(sigsnps_df) <- c("#Chr", "Position", "ID", "Ref", "Alt", "Test", "Obs_ct", "Beta", "Se", "T-stat", "P", "Block")
        for(i in 1:length(blocks)) {
            new.df <- dataframe[dataframe$Block == blocks[i] & dataframe$Test == "ADD",]
            if(nrow(new.df) != 0){
                sigsnps_df[i, ] <- new.df[new.df$P == min(new.df$P), ]
            }
        }
        sigsnps_df <- sigsnps_df[ !is.na(sigsnps_df$P), ]
    }

    if (domdev) {
        blocks <- unique(dataframe$Block)
        sigsnps_df <- data.frame(matrix(ncol = ncol(dataframe)))
        colnames(sigsnps_df) <- c("#Chr", "Position", "ID", "Ref", "Alt", "Test", "Obs_ct", "Beta", "Se", "T-stat", "P", "Block")
        for(i in 1:length(blocks)) {
            new.df <- dataframe[dataframe$Block == blocks[i] & dataframe$Test == "DOMDEV",]
            if(nrow(new.df) != 0){
                sigsnps_df[i, ] <- new.df[new.df$P == min(new.df$P), ]
            }
        }
        sigsnps_df <- sigsnps_df[ !is.na(sigsnps_df$P), ]
    }

    return(sigsnps_df)
}

# pulling the most significant snps for dominance deviation p-value
most_sigsnps <- most_sig_snps(sigsnps, avgeff = TRUE )

# extracting the list of significant snps
snp_list <- most_sigsnps$ID

# finding the corresponding additive effect sizes for the significant (via domdev) snps
avgeff_list <- rep(0,length(snp_list))

for (i in 1:length(snp_list)){
    if (nrow(sigsnps[sigsnps$ID == snp_list[i] & sigsnps$Test == "ADD", ]) != 0 ) {
        avgeff_list[i] <- sigsnps[sigsnps$ID == snp_list[i] & sigsnps$Test == "ADD", "Beta"]
    }
}

# finally, finding the mean additive effect sizes and dominance deviation values (I think these were "Beta" in the table) and printing.

mean_domdev <- mean(most_sigsnps$Beta)
mean_avgeff <- mean(avgeff_list)
print(paste("Mean D: ", mean_domdev))
print(paste("Mean A: ", mean_avgeff))

# Mean D: 1.03387707369682
# Mean A: 1.16120469982025


# Checking whether the additive effect sizes are significantly different from zero.

# Getting rid of NA rows
x <- raw.snps[!is.na(raw.snps$Beta), ]
add.df <- x[x$Test == "ADD", ]

write.table(add.df, "add_ttest_df.csv", sep = ",")

add.df <- read.table("add_ttest_df.csv", sep=",", stringsAsFactors = FALSE )

mean(add.df$Beta) # mean = 0.09516468

t.test(add.df$Beta, mu=0)

# RESULTS OF T-TEST
#	One Sample t-test
#
# data:  add.df$Beta
# t = 143.79, df = 638060, p-value < 2.2e-16
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#  0.09386749 0.09646187
# sample estimates:
#  mean of x
# 0.09516468
#
                                        #





