setwd ( "R.Projects" )
setwd ( "height_associations" )
setwd ( "ukb_height_gwas" )

chr.homoeff.analysis <- function ( chromosome, homoeff.sig.cutoff, avgeffplot = FALSE, homo.vs.domdev = TRUE ) { 
  
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
                          "domdev.pvalue" = geno.frame [ geno.frame$TEST == "DOMDEV", c( 9 ) ] ) 
  
  
  homoeff.sig <- i.frame [ i.frame$homoeff.pvalue < homoeff.sig.cutoff, ]
  
  if ( avgeffplot ) {
    plot <- plot( NULL , xlim = c ( -0.15, 0.15 ), ylim = c ( -0.15, 0.15 ), ylab="Homozygous/Dominance Effect", xlab="Avg Effect" 
                  , main = paste ( "Red= homozygous effect, blue = domdev, homoeff cutoff =", homoeff.sig.cutoff, ", chromosome =", chromosome ) )
    points ( homoeff.sig$avg.effect, homoeff.sig$homozygous.effect, pch = 20 , col = rgb ( 0.9 , 0.2 , 0 , 0.5 ) )
    points ( homoeff.sig$avg.effect, homoeff.sig$dominance.deviation, pch = 20, col = rgb ( 0 , 0.2 , 1 , 0.5  ) )
  }
  
  if ( homo.vs.domdev ) {
    plot ( homoeff.sig$homozygous.effect, homoeff.sig$dominance.deviation, main = paste ( "homozygous effect vs dominance deviation, homoeff cutoff =", homoeff.sig.cutoff, ", chromosome=", chromosome ), xlab = "homozygous effect size", ylab = "dominance deviation eff size", pch = 20 , col = rgb ( 0.9 , 0.2 , 0 , 0.5 ) )
  }
  
}


chr.homoeff.analysis( chromosome = 3, 
                     homoeff.sig.cutoff = 1,
                     avgeffplot = FALSE, 
                     homo.vs.domdev = TRUE )


