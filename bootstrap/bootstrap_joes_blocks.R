bootstrap.joes.blocks <- function ( avgeff.sig.cutoff, replicates, avgeff.thinning = FALSE, domdev.thinning = FALSE ) {

  setwd( "~/R.Projects/height_associations" )
  source( "height_association_functions.R" )
  
  domdev.frame <- data.frame( stringsAsFactors = FALSE )
  
  for ( i in 1:22 ) {
    avgeff.sig.frame <- make.cutoff.frame( i, avgeff.sig.cutoff = avgeff.sig.cutoff )
    domdev.frame.i <- cbind ( i, avgeff.sig.frame$position, avgeff.sig.frame$dominance.deviation, avgeff.sig.frame$domdev.pvalue, avgeff.sig.frame$pvalue )
    domdev.frame <- rbind.data.frame( domdev.frame, domdev.frame.i )
  }
  
  domdev.frame <- domdev.frame [ !is.na( domdev.frame$V3 ), ]
  domdev.frame <- cbind ( domdev.frame, "ID" = numeric( length = nrow( domdev.frame ) ) )
  
  setwd ( "~/R.Projects/height_associations/joes.blocks/EUR" )
  
  joes.blocks <- read.table ( "fourier_ls-all.bed", stringsAsFactors = FALSE , header = TRUE )
  
  joes.blocks$chr <- gsub ( 'chr' , '' , joes.blocks$chr )
  
  joes.blocks <- cbind (joes.blocks, "ID" = 1:nrow (joes.blocks) )
  
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
  
  domdev.blocks <- unique( domdev.frame$ID )
  
  if ( avgeff.thinning == FALSE & domdev.thinning == FALSE ) { 
  
    mean.replicates <- numeric()
    
    for ( i in 1:replicates ) {
      resampled.blocks <- sample ( domdev.blocks , size = length ( domdev.blocks ) , replace = TRUE )
      boot.replicate <- domdev.frame [ domdev.frame$ID %in% resampled.blocks, 3 ]
      mean.replicates [ i ] <- mean( boot.replicate )
    }
    
    return ( mean.replicates )
    
  }
  
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


