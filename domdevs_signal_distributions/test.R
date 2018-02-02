domdev_distributions <- function ( size, population, reps, generations, epsilonrange1, epsilonrange2, avgeff, domdev, loud = TRUE ) {

  epsilon <- runif ( size, epsilonrange1, epsilonrange2 ) # size refers to the number of loci to calculate the statistic over.

    Fst <-(generations/size) # approximation for drift

   # if ( size <= generations ) {
   #     stop ( "Population size should be greater than generations for accurate results to use normal approximation for drift.")
   # }


    var <- Fst * epsilon * ( 1 - epsilon )
    alpha <- (1 / 2) * avgeff + domdev * ( 1 - epsilon )
    afterdrift <- list ()
    Qx <- numeric()

    if ( loud ) {
        pb <- txtProgressBar( min = 0, max = reps )
    }

    for ( h in 1:reps) {

        if ( loud ) {
            setTxtProgressBar ( pb, h )
        }

        sum_pops <- numeric()

                                        # obtaining a list (length = populations) of vectors for the frequencies after drift for each loci given the Normal approximation to drift

        for ( i in 1:population ) {
            norm_drift <- numeric ()
            for ( j in 1:size ) {
                norm_drift [ j ] <- rnorm ( 1 , epsilon [ j ], sqrt ( var[ j ] ) )
            }
            afterdrift [[ i ]] <- norm_drift
        }

    # now we calculate one Qx statistic, summed across the given number of populations. We use matrix method developed in "equation_verification.R"

        for ( i in 1:population ) {

            p <- afterdrift [[ i ]]
            A <- rep ( avgeff, size )
            D <- rep ( domdev, size )

            sum_matrix <- matrix_maker ( size = size, A = A, D = D, p = p, epsilon = epsilon )

            sum1 <- sum ( diag ( sum_matrix ) )
            sum2 <- sum ( sum_matrix [ row ( sum_matrix ) != col ( sum_matrix ) ] )

            sum_pops [ i ] <- sum1 + sum2

        }

                                        # "correct" additive variance uses epsilon values, but in practice we can only take allele frequencies from present day populations, so choose arbitrary afterdrift population to calculate vA with (which should increase the variance?

        vA <- sum ( alpha^2  * afterdrift [[ 1 ]] * ( 1 - afterdrift[[ 1 ]] ) )

        Qx [ h ]  <-  ( 1 / ( Fst * vA ) ) * sum ( sum_pops )

    }

    if ( loud ) {
        hist( Qx, main = paste ( "domdev = ", domdev, ", avgeff = ", avgeff, ", populations = ", population, ", reps = ", reps ))
    }

    return ( sum_matrix)

}

z <- domdev_distributions ( size = 100,
                       population = 1,
                       reps = 1000,
                       generations = 10,
                       epsilonrange1 = 0.5,
                       epsilonrange2 = 0.5,
                       avgeff = 0.5,
                       domdev = -0.2 )

mean (z)
var (z)
