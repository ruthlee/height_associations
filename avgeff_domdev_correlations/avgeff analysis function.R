setwd ( "R.Projects" )
setwd ( "height_associations" )
setwd ( "ukb_height_gwas" )

chr.avgeff.analysis <- function ( chromosome, avgeff.sig.cutoff, avgeffplot = FALSE, homo.vs.domdev = TRUE ) {

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
                     "homozygous effect" = geno.frame [ geno.frame$TEST == "ADD", c( 7 ) ],
                     "homoeff pvalue" = geno.frame [ geno.frame$TEST == "ADD", c( 9 ) ],
                     "dominance deviation " = geno.frame [ geno.frame$TEST == "DOMDEV", c( 7 ) ],
                     "domdev pvalue" = geno.frame [ geno.frame$TEST == "DOMDEV", c( 9 ) ] )


  avgeff.sig <- i.frame [ i.frame$pvalue < avgeff.sig.cutoff, ]

  if ( avgeffplot ) {
    plot <- plot( NULL , xlim = c ( -0.15, 0.15 ), ylim = c ( -0.15, 0.15 ), ylab="Homozygous/Dominance Effect", xlab="Avg Effect"
                  , main = paste ( "Red= homozygous effect, blue = domdev, avgeff cutoff =", avgeff.sig.cutoff, ", chromosome =", chromosome ) )
    points ( avgeff.sig$avg.effect, avgeff.sig$homozygous.effect, pch = 20 , col = rgb ( 0.9 , 0.2 , 0 , 0.5 ) )
    points ( avgeff.sig$avg.effect, avgeff.sig$dominance.deviation, pch = 20, col = rgb ( 0 , 0.2 , 1 , 0.5  ) )
  }

  if ( homo.vs.domdev ) {
    plot ( avgeff.sig$homozygous.effect, avgeff.sig$dominance.deviation, main = paste ( "homozygous effect vs dominance deviation, avgeff cutoff =", avgeff.sig.cutoff, ", chromosome=", chromosome ), xlab = "homozygous effect size", ylab = "dominance deviation eff size", pch = 20 , col = rgb ( 0.9 , 0.2 , 0 , 0.5 ) )
  }

}


chr.avgeff.analysis( chromosome = 3,
                     avgeff.sig.cutoff = 1,
                     avgeffplot = TRUE,
                     homo.vs.domdev = TRUE )


