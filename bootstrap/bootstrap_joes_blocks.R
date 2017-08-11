bootstrap.joes.blocks <- function ( avgeff.sig.cutoff, replicates, avgeff.thinning = FALSE, domdev.thinning = FALSE ) {

  setwd( "~/R.Projects/height_associations" )
  source( "height_association_functions.R" )
  
  # creating empty domdev frame
  domdev.frame <- data.frame( stringsAsFactors = FALSE )
  
  # Filling domdev frame with all SNPs from chromosomes 1-22 
  for ( i in 1:22 ) {
    avgeff.sig.frame <- make.cutoff.frame( i, avgeff.sig.cutoff = avgeff.sig.cutoff )
    domdev.frame.i <- cbind ( i, avgeff.sig.frame$position, avgeff.sig.frame$dominance.deviation, avgeff.sig.frame$domdev.pvalue, avgeff.sig.frame$pvalue )
    domdev.frame <- rbind.data.frame( domdev.frame, domdev.frame.i )
  }
  
  # getting rid of NA's in column 3 (dominance deviation values) and adding row ID to each row 
  domdev.frame <- domdev.frame [ !is.na( domdev.frame$V3 ), ]
  domdev.frame <- cbind ( domdev.frame, "ID" = numeric( length = nrow( domdev.frame ) ) )
  
  setwd ( "~/R.Projects/height_associations/joes.blocks/EUR" )
  
  # reading in joe's blocks file
  joes.blocks <- read.table ( "fourier_ls-all.bed", stringsAsFactors = FALSE , header = TRUE )
  
  # chromosome ID to match domdev frame chromosome ID's 
  joes.blocks$chr <- gsub ( 'chr' , '' , joes.blocks$chr )
  
  # row ID creation
  joes.blocks <- cbind (joes.blocks, "ID" = 1:nrow (joes.blocks) )
  
  # progress bar for assigning joe's blocks numbers to each row in domdev frame (longest part of program) 
  pb <- txtProgressBar( min = 0, max = nrow ( domdev.frame ) , style= 3 )
  
  for ( i in 1:nrow ( domdev.frame ) ) {
    for ( j in 1:nrow ( joes.blocks ) ) { 
      if ( domdev.frame [ i, 2 ] >= joes.blocks [ j, 2 ] &
           domdev.frame [ i, 2 ] <= joes.blocks [ j, 3 ] & 
           domdev.frame [ i, 1 ] == joes.blocks [ j, 1 ] ) {
        domdev.frame [ i, 6 ] <- joes.blocks [ j, 4 ]
      }
    }
    setTxtProgressBar( pb, i )
  }
  
  # getting list of joe's blocks in the domdev frame 
  domdev.blocks <- unique( domdev.frame$ID )
  
  # without thinning, bootstrap and get list of mean replicates that you can histogram. 
  if ( avgeff.thinning == FALSE & domdev.thinning == FALSE ) { 
  
    mean.replicates <- numeric()
    
    for ( i in 1:replicates ) {
      resampled.blocks <- sample ( domdev.blocks , size = length ( domdev.blocks ) , replace = TRUE )
      boot.replicate <- domdev.frame [ domdev.frame$ID %in% resampled.blocks, 3 ]
      mean.replicates [ i ] <- mean( boot.replicate )
    }
    
    return ( mean.replicates )
    
  }
  
  # with domdev thinning, do same as above but first choose snps with lowest domdev pvalue (V4) in each block.
  # Most sig domdevs in most.sig.snps
  if ( domdev.thinning ) {
    
    most.sig.snps <- numeric() 
    mean.replicates <- numeric()
    
    for ( i in 1:length( domdev.blocks ) ) {
      snps <- domdev.frame [ domdev.frame$ID %in% domdev.blocks[ i ] , ]
      most.sig.snps[ i ] <- snps [ snps$V4 == min ( snps$V4 ), 3 ]
    }
    
    for ( i in 1:replicates ) {
      snp.sample <- sample ( most.sig.snps , size = length ( most.sig.snps ) , replace = TRUE )
      mean.replicates [ i ] <- mean ( snp.sample )
    }
    
    return ( mean.replicates )
    
  }
  
  # with avgeff thinning, do same as above but first choose snps with lowest avgeff pvalue (V5) in each block.
  # Most sig domdevs in most.sig.snps
  if ( avgeff.thinning ) {
    
    most.sig.snps <- numeric() 
    mean.replicates <- numeric()
    
    for ( i in 1:length( domdev.blocks ) ) {
      snps <- domdev.frame [ domdev.frame$ID %in% domdev.blocks[ i ] , ]
      most.sig.snps[ i ] <- snps [ snps$V5 == min ( snps$V5 ), 3 ]
    }
    
    for ( i in 1:replicates ) {
      snp.sample <- sample ( most.sig.snps , size = length ( most.sig.snps ) , replace = TRUE )
      mean.replicates [ i ] <- mean ( snp.sample )
    }
    
    return ( mean.replicates )     
  }
  
}

replicates <- bootstrap.joes.blocks ( avgeff.sig.cutoff = 10^-8,
                        replicates = 10000,
                        avgeff.thinning = TRUE ) 
  

hist( replicates, breaks = 100, main = paste ( "Most sig avgeff pvalue bootstrap, avgeff cutoff = 10^-8, reps = 10000, mean = ", mean (replicates) ) )


