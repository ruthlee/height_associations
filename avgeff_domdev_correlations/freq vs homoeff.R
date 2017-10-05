setwd ( "R.Projects" )
setwd ( "height_associations" )
setwd ( "ukb_height_gwas" )


chr.freq.vs.domdev <- function ( chromosome, homoeff.sig.cutoff, freq.vs.domdev.plot = TRUE ) { 

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
  
  first.set <- homoeff.sig [ homoeff.sig$dominance.deviation > -0.5 & homoeff.sig$dominance.deviation < -0.1 , ]
  second.set <- homoeff.sig [ homoeff.sig$dominance.deviation > -0.1 & homoeff.sig$dominance.deviation < 0 , ]
  third.set <- homoeff.sig [ homoeff.sig$dominance.deviation > 0 & homoeff.sig$dominance.deviation < 0.1 , ]
  fourth.set <- homoeff.sig [homoeff.sig$dominance.deviation > 0.1 & homoeff.sig$dominance.deviation < 0.5 , ] 

  if ( freq.vs.domdev.plot ) {
    plot <- plot( NULL , xlim = c ( 0, 1 ), ylim = c ( -0.5, 0.5 ), xlab="frequency", ylab="homozygous effect",  
                  main = paste ( "Red/orange= negative, blue/green = positive, homoeff cutoff =", homoeff.sig.cutoff, ", chromosome =", chromosome ) )
    points ( first.set$frequency, first.set$homozygous.effect, pch = 20 , col = rgb ( 1 , 0 , 0 , 0.4 ) )
    points ( second.set$frequency, second.set$homozygous.effect, pch = 20 , col = rgb ( 1 , 0.5 , 0 , 0.4 ) )
    points ( third.set$frequency, third.set$homozygous.effect, pch = 20 , col = rgb ( 0 , 1 , 0 , 0.4 ) )
    points ( fourth.set$frequency, fourth.set$homozygous.effect, pch = 20 , col = rgb ( 0 , 0 , 1 , 0.4 ) )
  }
}


chr.freq.vs.domdev ( chromosome = 1,
                     homoeff.sig.cutoff = 1,  # All SNPs
                     freq.vs.domdev.plot = TRUE )




