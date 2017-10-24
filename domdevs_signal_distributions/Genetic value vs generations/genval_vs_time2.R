genval_vs_time <- function ( size, replicates, generations, epsilonrange1, epsilonrange2, avgeff, domdev , genval = TRUE ) {

    Fst <- 1/size

    if ( genval == FALSE ) {
                                        # just the allele frequencies over generations

        all_freqs <- list ()

        plot( 0, main = paste ( "Allele frequency vs time, replicates = ", replicates, ", avgeff = ", avgeff, " , domdev = ", domdev ),  xlab = "Generations", ylab = "Allele frequency", xlim = c( 0, generations+1), ylim = c(0, 1), type = "l")

        for ( i in 1:replicates ) {
            afterdrift <- numeric()
            epsilon <- runif ( 1, epsilonrange1, epsilonrange2 )
            afterdrift [ 1 ]  <- epsilon

            for ( j in 1:(generations) ) {
                var <- Fst * afterdrift [ j ] * ( 1 - afterdrift [ j ] )
                afterdrift [ j + 1 ] <-  rnorm( 1, afterdrift [ j ], sqrt ( var ) )

                afterdrift  <- replace ( afterdrift, is.nan(afterdrift) | afterdrift <= 0, 0 )
                afterdrift <- replace ( afterdrift, afterdrift >= 1, 1 )

            }

            lines ( seq_len ( generations+1 ), afterdrift, col = rgb(190, 190, 190, alpha = 100, maxColorValue = 255 ) )

            all_freqs [[ i ]] <- afterdrift

        }

        all_freqs_df <- as.data.frame( all_freqs )
        all_freqs_df$avg = rowMeans( all_freqs_df )

        mean_freqs <- all_freqs_df$avg

        lines ( seq_len ( generations + 1 ), mean_freqs, col = "red", lwd = 2 )

        return( all_freqs_df )
    }

    ####

    if ( genval == TRUE ) {

        plot( 0, main =  paste ( "Genetic value vs time, replicates = ", replicates, ", avgeff = ", avgeff, " , domdev = ", domdev ),  xlab = "Generations", ylab = "Genetic value", xlim = c( 0, generations+1), ylim = c(-1, 1), type = "l")

        all_genvals <- list()

        for ( i in 1:replicates ) {
            afterdrift <- numeric()
            genval <- numeric()
            epsilon <- runif ( 1, epsilonrange1, epsilonrange2 )
            afterdrift [ 1 ]  <- epsilon

            for ( j in 1:(generations) ) {
                var <- Fst * afterdrift [ j ] * ( 1 - afterdrift [ j ] )
                alpha <- ( 1 / 2 ) * avgeff + domdev * ( 1 - 2 * afterdrift [ j ] )

                afterdrift  <- replace ( afterdrift, is.nan(afterdrift) | afterdrift <= 0, 0 )

                genval [ j ] <- ( afterdrift [ j ] - epsilon ) * alpha

                afterdrift [ j + 1 ] <-  rnorm( 1, afterdrift [ j ], sqrt ( var ) )

                afterdrift  <- replace ( afterdrift, is.nan(afterdrift) | afterdrift <= 0, 0 )
                afterdrift <- replace ( afterdrift, afterdrift >= 1, 1 )

            }

            lines ( seq_len ( generations ), genval, col = rgb(190, 190, 190, alpha = 100, maxColorValue = 255 ) )

            all_genvals [[ i ]] <- genval

        }

        all_genvals_df <- as.data.frame ( all_genvals )
        all_genvals_df$avg = rowMeans ( all_genvals_df )

        mean_genvals <- all_genvals_df$avg

        lines ( seq_len ( generations ), mean_genvals, col = "red", lwd = 2 )

        return( all_genvals_df )

    }

}



genval_vs_time (size = 100,
                replicates = 100,
                generations = 5000,
                epsilonrange1 = 0.5,
                epsilonrange2 = 0.5,
                avgeff = 0,
                domdev = 1,
                genval = TRUE )




domdev <- -seq( 0.2, 1, 0.2 )

for ( i in 1:length(domdev)) {
    pdf ( paste ( "domdev=", domdev[i], sep = "" ) )
    x <- genval_vs_time (size = 10000,
                         replicates = 100,
                         generations = 5000,
                         epsilonrange1 = 0.5,
                         epsilonrange2 = 0.5,
                         avgeff = 0.5,
                         domdev = domdev [ i ],
                         genval = TRUE )
    dev.off()
}


