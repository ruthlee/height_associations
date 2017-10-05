setwd( "~/R.Projects/height_associations" )
source( "qq.plot.pvalues.function.R" )

frame <- qq.plot ( chromosome = 1,
                   avgeff.sig.cutoff = 10^-8,
                   genotypic.cutoff = FALSE,
                   avgeff.plot = FALSE,
                   homoeff.plot = FALSE,
                   domdev.plot = TRUE )


# Created 2 plots for every chromosome 

for ( i in c ( 10^-8, 1 ) ) {
  for ( j in 1:22 )
    qq.plot ( chromosome = j, 
              avgeff.sig.cutoff = i,
              avgeff.plot = FALSE,
              homoeff.plot = FALSE,
              domdev.plot = TRUE )
}



# All chromosomes on one plot. 

plot ( NULL, main = "domdev qq plot, all chromosomes", xlab = "-log expected pvalues", ylab = "-log domdev pvalues", 
       xlim = c( 0, 5 ), ylim = c( 0, 5 ) )

colors <- rainbow ( 22, alpha = 0.3 )

for ( i in 1:22 ) {
  pvalues <- qq.plot ( chromosome = i,
                       avgeff.sig.cutoff = 10^-8,
                       avgeff.plot = FALSE,
                       homoeff.plot = FALSE,
                       domdev.plot = FALSE )
  logexpvalues <- - log ( ppoints ( length ( pvalues ) ), 10 )
  points ( logexpvalues, pvalues, pch = 20, col = colors [ i ] )
  abline ( 0, 1 )
}



# All SNPs from one chromosome in ONE qq-plot

setwd ( "/Users/Kyelee/R.Projects/height_associations")
source ( "height_association_functions.R" )

sig.frame <- qq.plot ( 1, 22, avgeff.sig.cutoff = 10^-8, domdev.plot = TRUE )

# finding mean dominance deviation effect size

mean( sig.frame$dominance.deviation, na.rm = TRUE )

# mean for significance of 10^-8: 0.003135951
# mean for significance of 1: 0.0006821619

hist ( domdev, breaks = 100, main = "Dominance Deviations Effect Sizes, all chromosomes, avg eff size sig = 10^-8" )

