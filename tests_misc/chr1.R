setwd ( "R.Projects" )
setwd ( "height_associations" )
setwd ( "ukb_height_gwas" )
chr1.frq <- read.table ( "chr1.frq" , stringsAsFactors = FALSE , header = TRUE )
height.chr1 <- read.table ( "height_add_chr1.assoc.linear" , stringsAsFactors = FALSE , header = TRUE )
geno.chr1 <- read.table ( "height_geno_chr1.assoc.linear" , stringsAsFactors = FALSE , header = TRUE )


# c-binding all the parts of the tables that we need into one data frame -- NOT a 
# matrix because we need to do create a logical vector of TRUE's and FALSE's for when
# p is less than 10^-8. The dollar sign calls the particular column in the data frame.
# Once we have the vector, we can make it the row of the data frame using the bracket 
# notation. The logical vector will make all the rows in which FALSE appears to disappear.
# This is because we are assigning a new data frame in which the only rows called are the
# ones that are TRUE. 



homozygous.effect <- geno.chr1 [ 1:nrow ( geno.chr1 ) %% 3 == 1 | 1:nrow ( geno.chr1 ) %/% 3 == 1/3, 7 ]
dominance.deviation <- geno.chr1 [ 1:nrow ( geno.chr1 ) %% 3 == 2 | 1:nrow ( geno.chr1 ) %/% 3 == 2/3, 7 ]
p.homo.eff <- geno.chr1 [ 1:nrow ( geno.chr1 ) %% 3 == 1 | 1:nrow ( geno.chr1 ) %/% 3 == 1/3, 9 ]
p.dom.dev <- geno.chr1 [ 1:nrow ( geno.chr1 ) %% 3 == 1 | 1:nrow ( geno.chr1 ) %/% 3 == 1/3, 9 ]

# Alternatively, could have used logical vector geno.chr1$TEST == "Add" and called data
# frame geno.chr1 [ geno.chr1$TEST == 'ADD', ]  

chr1 <- data.frame ( "position" = height.chr1 [ , c ( 3 ) ], 
                "frequency" = chr1.frq [ , c ( 5 ) ],
                "avg.effect" = height.chr1 [ , c( 7 ) ],
                "pvalue" = height.chr1 [ , c( 9 ) ],
                "homozygous effect" = homozygous.effect, 
                "homoeff pvalue" = p.homo.eff,
                "dominance deviation " = dominance.deviation,
                "domdev pvalue" = p.dom.dev ) 

avgeff.sig.chr1 <- chr1 [ chr1$pvalue < 10^-10, ]

plot <- plot( NULL , xlim = c ( -0.15, 0.15 ), ylim = c ( -0.15, 0.15 ), ylab="Homozygous/Dominance Effect", xlab="Avg Effect" 
              , main = "Red= homozygous effect, blue = domdev" )

points ( avgeff.sig.chr1$avg.effect, avgeff.sig.chr1$homozygous.effect, pch = 20 , col = rgb ( 0.9 , 0.2 , 0 , 0.5 ) )
points ( avgeff.sig.chr1$avg.effect, avgeff.sig.chr1$dominance.deviation, pch = 20, col = rgb ( 0 , 0.2 , 1 , 0.5  ) )

plot ( avgeff.sig.chr1$homozygous.effect, avgeff.sig.chr1$dominance.deviation, main = "homozygous effect vs dominance deviation", xlab = "homozygous effect size", ylab = "dominance deviation eff size", pch = 20 , col = rgb ( 0.9 , 0.2 , 0 , 0.5 ) )




# looking at frequency of SNPs below each p-value cutoff. 

snpcounts <- numeric()
for ( i in 0:10 ) {
  snpcounts [ i + 1 ] <- nrow ( chr1 [ chr1$pvalue < 10^(-i), ] )
}
  
plot ( c(0:10), snpcounts, xlab = "-log(pthreshold)", ylab = "Number of SNPs", pch = 20, main = "SNP Frequency" )

print(snpcounts)



sorted.frame <- chr1 [ order ( chr1$pvalue ), ]
head(sorted.frame)


