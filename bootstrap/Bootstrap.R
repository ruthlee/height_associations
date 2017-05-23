# mean for significance of 10^-8: 0.003135951
# mean for significance of 1: 0.0006821619

setwd( "~/R.Projects/height_associations/boot/R" )

source ( "bootfuns.q" )

domdev <- (append.cutoff.frame( 1, 22, 10^-8 ))$dominance.deviation

domdev <- domdev [ !is.na(domdev) ]

samplemean <- function( domdev, d) {
  
  return ( mean ( domdev[d] ) )
  
}

boot.data <- boot ( data = domdev,
       statistic = samplemean,
       R = 1000 )

print(boot.data)
plot(boot.data)
boot.ci( boot.data,
         conf = 0.95,
         type= "norm" )

mean (domdev)
