make.i.frame <- function ( chromosome ) {
  setwd( "~/R.Projects/height_associations/ukb_height_gwas" )
  
  i <- chromosome 
  freq.frame <- paste ( "chr", i, ".frq", sep = "" )
  add.frame <- paste ( "height_add_chr", i, ".assoc.linear", sep = "" )
  geno.frame <- paste ( "height_geno_chr", i, ".assoc.linear", sep = "" )
  
  freq.frame <- read.table ( freq.frame , stringsAsFactors = FALSE , header = TRUE )
  add.frame <- read.table ( add.frame , stringsAsFactors = FALSE , header = TRUE )
  geno.frame <- read.table ( geno.frame , stringsAsFactors = FALSE , header = TRUE )
  
  i.frame <- data.frame ( "position" = add.frame [ , c ( 3 ) ], 
                          "frequency" = freq.frame [ , c ( 5 ) ],
                          "avg.effect" = add.frame [ , c( 7 ) ],
                          "pvalue" = add.frame [ , c( 9 ) ],
                          "homozygous.effect" = geno.frame [ geno.frame$TEST == "ADD", c( 7 ) ], 
                          "homoeff.pvalue" = geno.frame [ geno.frame$TEST == "ADD", c( 9 ) ],
                          "dominance.deviation " = geno.frame [ geno.frame$TEST == "DOMDEV", c( 7 ) ],
                          "domdev.pvalue" = geno.frame [ geno.frame$TEST == "DOMDEV", c( 9 ) ],
                          "genotypic.test" = geno.frame [ geno.frame$TEST == "GENO_2DF", c( 7 ) ],
                          "genotypic.pvalue" = geno.frame [ geno.frame$TEST == "GENO_2DF", c( 9 ) ] )
  
  
  return ( i.frame )
  
}
# Creates a data frame with position, frequency, avg.effect, avg effect pvalue,
# homozygous effect, homoeff pvalue, domdev, domdev pvalue, genotypic effect
# size, and genotypic p values as columns with all SNPS from chromosome i. 

make.cutoff.frame <- function ( chromosome, avgeff.sig.cutoff = FALSE, domdev.sig.cutoff = FALSE ) {
  
  i.frame <- make.i.frame( chromosome )
  
  if ( avgeff.sig.cutoff ) {
    avgeff.sig <- i.frame [ i.frame$pvalue < avgeff.sig.cutoff, ]
    return ( avgeff.sig )
  } else if ( domdev.sig.cutoff ) {
    domdev.sig <- i.frame [ i.frame$domdev.pvalue < domdev.sig.cutoff, ]
    return ( domdev.sig )
  }
  
}
# Creates a data frame from i.frame from one chromosome based on cutoffs for 
# average effect size or dominance deviation. 

append.cutoff.frame <- function ( start.chr, end.chr, avgeff.sig.cutoff = FALSE, domdev.sig.cutoff = FALSE ) {
  avgeff.sig.frame <- data.frame( stringsAsFactors = FALSE ) 
  if ( avgeff.sig.cutoff ) {
    for ( i in start.chr:end.chr ) {
      avgeff.sig.frame.i <- make.cutoff.frame( chromosome = i, avgeff.sig.cutoff = 10^-8 )
      avgeff.sig.frame <- rbind.data.frame( avgeff.sig.frame, avgeff.sig.frame.i )
    }
    return ( avgeff.sig.frame )
    
  } else if ( domdev.sig.cutoff ) {
    domdev.sig.frame.i <- make.cutoff.frame( chromosome = i, domdev.sig.cutoff = 10^-8 )
    domdev.sig.frame <- rbind.data.frame( domdev.sig.frame, domdev.sig.frame.i )
    return ( domdev.sig.frame )
  }
  
}
# Creates a data frame from make.cutoff.frame from a specified range of chromosomes
# based on cutoffs for average effect size or dominance deviation. 

qq.plot <- function ( start.chr, end.chr, avgeff.sig.cutoff = FALSE , genotypic.cutoff = FALSE , avgeff.plot = FALSE, homoeff.plot = FALSE, domdev.plot = FALSE ) {
  
  #recover()
  setwd( "~/R.Projects/height_associations/ukb_height_gwas" )
  
  i.frame <- append.cutoff.frame( start.chr, end.chr, avgeff.sig.cutoff, genotypic.cutoff ) 
  
  if ( avgeff.sig.cutoff ) {
    avgeff.sig <- i.frame 
    
    logexpvalues <- - log ( ppoints ( nrow ( avgeff.sig ) ), 10 )
    
    if ( avgeff.plot ) {
      sorted.frame <- avgeff.sig [ order( avgeff.sig$pvalue ), ]
      logpvalues <- - log ( sorted.frame$pvalue, 10 )
      plot ( logexpvalues, logpvalues, main = "average effect qq plot", xlab = "-log expected pvalues", ylab = "-log avg eff pvalues", 
             pch = 20, col = rgb ( 0, 0, 0, 0.2 ) )
      abline ( 0, 1 )
    }
    
    if ( homoeff.plot ) {
      sorted.frame <- avgeff.sig [ order( avgeff.sig$homoeff.pvalue ), ]
      logpvalues <- - log ( sorted.frame$homoeff.pvalue, 10 )
      plot ( logexpvalues, logpvalues, main = "homoeff qq plot", xlab = "-log expected pvalues", ylab = "-log homoeff pvalues",
             pch = 20, col = rgb ( 0, 0, 0, 0.2 ) )
      abline ( 0, 1 )
    }
    
    if ( domdev.plot ) {
      sorted.frame <- avgeff.sig [ order( avgeff.sig$domdev.pvalue ), ]
      logpvalues <- - log ( sorted.frame$domdev.pvalue, 10 )
      plot ( logexpvalues, logpvalues, main = "domdev qq plot", xlab = "-log expected pvalues", ylab = "-log domdev pvalues", 
             pch = 20, col = rgb ( 0, 0, 0, 0.2 ) )
      abline ( 0, 1 )
    }
    
  }
  
  if ( genotypic.cutoff ) {
    
    genotypic.sig <- i.frame 
    
    logexpvalues <- - log ( ppoints ( nrow ( genotypic.sig ) ), 10 )
    
    if ( avgeff.plot ) {
      sorted.frame <- genotypic.sig [ order( genotypic.sig$pvalue ), ]
      logpvalues <- - log ( sorted.frame$pvalue, 10 )
      plot ( logexpvalues, logpvalues, main = "average effect qq plot", xlab = "-log expected pvalues", ylab = "-log avg eff pvalues", 
             pch = 20, col = rgb ( 0, 0, 0, 0.2 ) )
      abline ( 0, 1 )
    }
    
    if ( homoeff.plot ) {
      sorted.frame <- genotypic.sig [ order( genotypic.sig$homoeff.pvalue ), ]
      pdf ( file = paste ( "homoeff qq plot, chromosome ", chromosome, ", genotypic.sig cutoff = ", genotypic.cutoff ) )
      plot ( logexpvalues, logpvalues, main = "homoeff qq plot", xlab = "-log expected pvalues", ylab = "-log homoeff pvalues",
             pch = 20, col = rgb ( 0, 0, 0, 0.2 ) )
      abline ( 0, 1 )
    }
    
    if ( domdev.plot ) {
      sorted.frame <- genotypic.sig [ order( genotypic.sig$domdev.pvalue ), ]
      logpvalues <- - log ( sorted.frame$domdev.pvalue, 10 )
      plot ( logexpvalues, logpvalues, main = "domdev qq plot", xlab = "-log expected pvalues", ylab = "-log domdev pvalues", 
             pch = 20, col = rgb ( 0, 0, 0, 0.2 ) )
      abline ( 0, 1 )
    }
    
  }
  
  return ( i.frame )
  
}
# makes qq plots of SNPs from a chromosome against avg eff pvalues, domdev pvalues
# and homoeff pvalues. 


joes.blocks.thinning <- function ( avgeff.sig.cutoff, avgeff.thinning = FALSE, domdev.thinning = FALSE ) {
  
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
  
  if ( domdev.thinning ) {
    
    most.sig.snps <- numeric() 
    
    for ( i in 1:length( domdev.blocks ) ) {
      snps <- domdev.frame [ domdev.frame$ID %in% domdev.blocks[ i ] , ]
      most.sig.snps[ i ] <- snps [ snps$V4 == min ( snps$V4 ), 3 ]
    }
    
    return ( most.sig.snps )
    
  }
  
  if ( avgeff.thinning ) {
    
    most.sig.snps <- numeric() 
    
    for ( i in 1:length( domdev.blocks ) ) {
      snps <- domdev.frame [ domdev.frame$ID %in% domdev.blocks[ i ] , ]
      most.sig.snps[ i ] <- snps [ snps$V5 == min ( snps$V5 ), 3 ]
    }
    
    return ( most.sig.snps )
    
  }
  
}

# Thins based on average effect pvalue or dominance deviation pvalue, selecting
# the lowest pvalue in each of joe's blocks with a significant SNP based on an
# average effect size pvalue. 




