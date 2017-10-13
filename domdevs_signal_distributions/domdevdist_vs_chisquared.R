cutoff <- qchisq(0.95, 1)
sum ( x > cutoff ) / length(x)


# From the list of Qx statistic functions we got from the domdevdist function, we can see what proportion of those values are greater than a certain chi-squared statistic with a user-specified p-value. This is to see how changing domdev values for the Qx statistic changes the chi-squared distribution-- we expect a greater proportion of Qx statistic values to be greater than the expected proportion as domdev increases.

prop_over_exp.p <- function ( Qx = numeric(), pvalue, exp.df ) {
    cutoff <- qchisq( 1-pvalue, exp.df )
    sum ( Qx > cutoff ) / length ( x )
}

# Now we can make the plot of domdev values vs proportion over expected pvalue cutoff.

domdev_vs_exp.p <- function ( domdev = numeric(), pvalue, exp.df, size = 100, population = 1, reps = 1000, generations = 10, epsilonrange1 = 0.4, epsilonrange2 =0.6, avgeff = 0.5 ) {

    prop <- numeric()

    pb <- txtProgressBar ( min = 0, max = length(domdev) )

    for ( i in 1:length(domdev) ) {

        setTxtProgressBar ( pb, i )

        Qx <- domdev_distributions ( size,
                              population,
                              reps,
                              generations,
                              epsilonrange1,
                              epsilonrange2,
                              avgeff,
                              domdev = domdev [ i ],
                              loud = FALSE )

        prop [ i ] <- prop_over_exp.p ( Qx, pvalue, exp.df )
    }

    plot ( domdev, prop, main = paste( "Domdev of distribution vs proportion of distributions over p = ", pvalue, " threshold"), xlab = "domdev", ylab = "proportion", pch = 16 )

    return( prop )

}

domdev <- seq ( 0, 1, 0.1 )

# Testing the function with population = 1

prop <- domdev_vs_exp.p ( domdev, 0.05, 1 )


pb <- txtProgressBar( min = 0 , max = length(domdev))

for ( i in 1:length(domdev) ) {
    setTxtProgressBar(pb, i )
    x <- domdev_distributions ( size = 100,
                       population = 1,
                       reps = 1000,
                       generations = 10,
                       epsilonrange1 = 0.4,
                       epsilonrange2 = 0.6,
                       avgeff = 0.5,
                       domdev = domdev[i],
                       loud = FALSE)
    pdf ( paste ( "plot", i, ".pdf", sep = "") )
    plot(ecdf(x))
    plot(ecdf(rchisq(1000, df=1)), add= TRUE, col= "red")
    dev.off()
}







