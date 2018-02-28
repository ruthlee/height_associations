setwd("~/R.Projects/height_associations/all_chromosomes_gwas")

# Read in all significant snps
sigsnps <- read.table("sigsnps_blocks.csv", sep=",", stringsAsFactors = FALSE )
colnames(sigsnps) <- c("Chr", "Position", "ID", "Min", "Maj", "Test", "Obs_ct", "Beta", "Se", "T-stat", "P", "Block")

# Read in dominance deviation significant snps
domdev_sigsnps <- read.table("domdev_sigsnps.csv", stringsAsFactors=FALSE)
colnames(domdev_sigsnps) <- c("Chr", "Position", "ID", "Min", "Maj", "Test", "Obs_ct", "Beta", "Se", "T-stat", "P", "Block")

# Read in frequency data
chr1_freqs <- read.table("allele_freqs/chr1_allele_freq_GBR.frq", row.names=NULL, stringsAsFactors = FALSE)
colnames(chr1_freqs) <- c("Chr", "Pos", "x", "y", "MajAlleleFreq", "MinAlleleFreq")

# Configure allele frequency columns in raw frequency table
for (i in 1:nrow(chr1_freqs)) {
    chr1_freqs[i, "MajAlleleFreq"]  <- gsub("[AGCT:]", "", chr1_freqs[i, "MajAlleleFreq"])
    chr1_freqs[i, "MinAlleleFreq"]  <- gsub("[AGCT:]", "", chr1_freqs[i, "MinAlleleFreq"])
}
chr1_freqs$MajAlleleFreq <- as.numeric(chr1_freqs$MajAlleleFreq)
chr1_freqs$MinAlleleFreq <- as.numeric(chr1_freqs$MinAlleleFreq)

# Get additive effect sizes by finding same SNPs as domdev_sigsnps
add_sigsnps <- sigsnps[sigsnps$Test == "ADD" & sigsnps$ID == domdev_sigsnps$ID[1], ]
for (i in 2:nrow(domdev_sigsnps)) {
    add_sigsnps <- rbind(add_sigsnps, sigsnps[sigsnps$Test == "ADD" & sigsnps$ID == domdev_sigsnps$ID[i], ])
}

# at this point we have chr1_freq df configured, add_sigsnps, domdev_sigsnps.

# make domdev_sigsnps and add_sigsnps match up in SNPs
domdev_sigsnps <- domdev_sigsnps[ domdev_sigsnps$ID %in% add_sigsnps$ID, ]

# get parameter list
A <- add_sigsnps[add_sigsnps$Chr==1, "Beta"]
D <- domdev_sigsnps[domdev_sigsnps$Chr==1, "Beta"]

# match positions with frequencies and add to domdev_sigsnps data frame
chr1_freqs <- chr1_freqs[ chr1_freqs$Pos %in% domdev_sigsnps$Position, ]

# now we plug these parameters into the summation equation for Qx, using the frequency data

get_Qx <- function(A, D, p) {
    mat <- matrix(nrow=length(p), ncol=length(p))
    epsilon <- mean(p)
    for( i in 1:length(p)){
         mat[ i , i ] <- ((1/2) * A[i] * (p[i] - epsilon))^2
                      + (A[i] * D[i] * (1 - 2 * p[i]) * (p[i] - epsilon)^2)
                      + (D[i] * (1-2 * p[i]) * (p[i] - epsilon))^2
         for (j in 1:length(p)) {
              if ( i != j ) {
                  mat[i, j] <- ((1/4) * A[i] * A[j] * (p[i] - epsilon) *  (p[j] - epsilon)) + (1/2) * A[i] * D[j] * (1 - 2 * p[j]) * (p[i] - epsilon) * (p[j] - epsilon) + (1/2) *  D[i] * A[j] * (1 - 2 * p[i]) * (p[i] - epsilon) * (p[j] - epsilon) +  D[i] * (1 - 2 * p[i]) * (p[i] - epsilon) * D[j] * (1 - 2 * p[j]) * (p[j] - epsilon)
              }
         }
    }

    Qx <- sum(diag(mat)) + sum(mat[row(mat) != col(mat)])
    return(Qx)
}

get_Qx(A, D, p = chr1_freqs$MinAlleleFreq)

# Qx = 21.5324 (MajAlleleFreq)
# Qx = 1608.5324 (MinAlleleFreq)

# now to iterate over all chromosome frequencies

get_freqs <- function(population) {
    chr_freqs <- data.frame(matrix(ncol=6, nrow=0))
    colnames(chr_freqs) <- c("Chr", "Pos", "x", "y", "MajAlleleFreq", "MinAlleleFreq")

    for (i in c(1:8,10:11,13:14,16:18,20:22)){
        chr_i_freqs <- read.table(paste("allele_freqs/chr", i, "_allele_freq_", population, ".frq", sep = ""), row.names = NULL, stringsAsFactors = FALSE)
        colnames(chr_i_freqs) <- c("Chr", "Pos", "x", "y", "MajAlleleFreq", "MinAlleleFreq")
        chr_freqs <- rbind(chr_freqs, chr_i_freqs)
    }
    return(chr_freqs)
}

# getting frequency list for all chromosomes in GBR and TSI populations
gbr.freqs <- get_freqs("GBR")
tsi.freqs <- get_freqs("TSI")

# pruning sigsnps lists to only snps in our frequency lists
domdev_sigsnps <- domdev_sigsnps[domdev_sigsnps$Position %in% gbr.freqs$Pos, ]
add_sigsnps <- add_sigsnps[add_sigsnps$Position %in% gbr.freqs$Pos, ]

gbr.freqs <- gbr.freqs[ gbr.freqs$Pos %in% add_sigsnps$Position, ]
tsi.freqs <- tsi.freqs[ tsi.freqs$Pos %in% add_sigsnps$Position, ]

# Configure allele frequency columns in raw frequency table
for (i in 1:nrow(gbr.freqs)) {
    gbr.freqs[i, "MajAlleleFreq"]  <- gsub("[AGCT:]", "", gbr.freqs[i, "MajAlleleFreq"])
    gbr.freqs[i, "MinAlleleFreq"]  <- gsub("[AGCT:]", "", gbr.freqs[i, "MinAlleleFreq"])
    tsi.freqs[i, "MajAlleleFreq"]  <- gsub("[AGCT:]", "", tsi.freqs[i, "MajAlleleFreq"])
    tsi.freqs[i, "MinAlleleFreq"]  <- gsub("[AGCT:]", "", tsi.freqs[i, "MinAlleleFreq"])
}
gbr.freqs$MajAlleleFreq <- as.numeric(gbr.freqs$MajAlleleFreq)
gbr.freqs$MinAlleleFreq <- as.numeric(gbr.freqs$MinAlleleFreq)
tsi.freqs$MajAlleleFreq <- as.numeric(tsi.freqs$MajAlleleFreq)
tsi.freqs$MinAlleleFreq <- as.numeric(tsi.freqs$MinAlleleFreq)

A <- add_sigsnps$Beta
D <- domdev_sigsnps$Beta
p <- gbr.freqs$MinAlleleFreq

gbr_Qx <- get_Qx(A, D, p)
# 284394.9

A <- add_sigsnps$Beta
D <- domdev_sigsnps$Beta
p <- tsi.freqs$MinAlleleFreq

tsi_Qx <- get_Qx(A, D, p)
# 277076.7

