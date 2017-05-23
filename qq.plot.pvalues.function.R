qq.plot <- function ( chromosome, avgeff.sig.cutoff = FALSE , genotypic.cutoff = FALSE , avgeff.plot = FALSE, homoeff.plot = FALSE, domdev.plot = FALSE ) {
  
  #recover()
  setwd( "~/R.Projects/height_associations/ukb_height_gwas" )
  
  i.frame <- make.i.frame( chromosome ) 
  
  if ( avgeff.sig.cutoff ) {
    avgeff.sig <- i.frame [ i.frame$pvalue < avgeff.sig.cutoff, ]
    
    logexpvalues <- - log ( ppoints ( nrow ( avgeff.sig ) ), 10 )
    
    if ( avgeff.plot ) {
      sorted.frame <- avgeff.sig [ order( avgeff.sig$pvalue ), ]
      logpvalues <- - log ( sorted.frame$pvalue, 10 )
      plot ( logexpvalues, logpvalues, main = "average effect qq plot", xlab = "-log expected pvalues", ylab = "-log avg eff pvalues", 
             pch = 20, col = rgb ( 0, 0, 0, 0.2 ), xlim = c( 0, 10 ), ylim = c( 0, 80 ) )
      abline ( 0, 1 )
    }
    
    if ( homoeff.plot ) {
      sorted.frame <- avgeff.sig [ order( avgeff.sig$homoeff.pvalue ), ]
      logpvalues <- - log ( sorted.frame$homoeff.pvalue, 10 )
      plot ( logexpvalues, logpvalues, main = "homoeff qq plot", xlab = "-log expected pvalues", ylab = "-log homoeff pvalues",
             pch = 20, col = rgb ( 0, 0, 0, 0.2 ), xlim = c( 0, 10 ), ylim = c( 0, 80 ) )
      abline ( 0, 1 )
    }
    
    if ( domdev.plot ) {
      sorted.frame <- avgeff.sig [ order( avgeff.sig$domdev.pvalue ), ]
      logpvalues <- - log ( sorted.frame$domdev.pvalue, 10 )
      plot ( logexpvalues, logpvalues, main = "domdev qq plot", xlab = "-log expected pvalues", ylab = "-log domdev pvalues", 
             pch = 20, col = rgb ( 0, 0, 0, 0.2 ), xlim = c( 0, 10 ), ylim = c( 0, 10 ) )
      abline ( 0, 1 )
    }
    
    sorted.frame <- avgeff.sig [ order( avgeff.sig$domdev.pvalue ), ]
    logpvalues <- - log ( sorted.frame$domdev.pvalue, 10 )
    return ( sorted.frame )
    
  }
  
  if ( genotypic.cutoff ) {
    
    genotypic.sig <- i.frame [ i.frame$genotypic.pvalue < genotypic.cutoff, ]
    
    logexpvalues <- - log ( ppoints ( nrow ( genotypic.sig ) ), 10 )
    
    if ( avgeff.plot ) {
      sorted.frame <- genotypic.sig [ order( genotypic.sig$pvalue ), ]
      logpvalues <- - log ( sorted.frame$pvalue, 10 )
      plot ( logexpvalues, logpvalues, main = "average effect qq plot", xlab = "-log expected pvalues", ylab = "-log avg eff pvalues", 
             pch = 20, col = rgb ( 0, 0, 0, 0.2 ), xlim = c( 0, 10 ), ylim = c( 0, 80 ) )
      abline ( 0, 1 )
    }
    
    if ( homoeff.plot ) {
      sorted.frame <- genotypic.sig [ order( genotypic.sig$homoeff.pvalue ), ]
      pdf ( file = paste ( "homoeff qq plot, chromosome ", chromosome, ", genotypic.sig cutoff = ", genotypic.cutoff ) )
      plot ( logexpvalues, logpvalues, main = "homoeff qq plot", xlab = "-log expected pvalues", ylab = "-log homoeff pvalues",
             pch = 20, col = rgb ( 0, 0, 0, 0.2 ), xlim = c( 0, 10 ), ylim = c( 0, 80 ) )
      abline ( 0, 1 )
    }
    
    if ( domdev.plot ) {
      sorted.frame <- genotypic.sig [ order( genotypic.sig$domdev.pvalue ), ]
      logpvalues <- - log ( sorted.frame$domdev.pvalue, 10 )
      plot ( logexpvalues, logpvalues, main = "domdev qq plot", xlab = "-log expected pvalues", ylab = "-log domdev pvalues", 
             pch = 20, col = rgb ( 0, 0, 0, 0.2 ), xlim = c( 0, 10 ), ylim = c( 0, 10 ) )
      abline ( 0, 1 )
    }
    
  }
  
}

qq.plot ( chromosome = 1,
          avgeff.sig.cutoff = 10^-8,
          genotypic.cutoff = FALSE,
          avgeff.plot = FALSE,
          homoeff.plot = FALSE,
          domdev.plot = TRUE )




