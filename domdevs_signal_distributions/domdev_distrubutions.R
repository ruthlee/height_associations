# We want to directly simulate the Qx statistic using the terms derived mathematically. We can do this by simulating a random distribution of ancestral allele frequencies over a given number of populations.

domdev_distributions <- function ( size, population, reps, generations, epsilonrange1, epsilonrange2, avgeff, domdev ) {

  epsilon <- runif ( size, epsilonrange1, epsilonrange2 ) # size refers to the number of loci to calculate the statistic over.
  drift <- generations / size

  if ( size <= generations ) {
    stop ( "Population size should be greater than generations for accurate results to use normal approximation for drift.")
  }


  var <- drift * epsilon * ( 1 - epsilon )
  afterdrift <- list ()
  Qx <- numeric()

  pb <- txtProgressBar( min = 0, max = reps )

  for ( h in 1:reps) {
      setTxtProgressBar ( pb, h )
      sum_pops <- numeric()

                                        # obtaining a list (length = populations) of vectors for the frequencies after drift for each loci given the Normal approximation to drift

    for ( i in 1:population ) {
        drift <- numeric ()
        for ( j in 1:size ) {
            drift [ j ] <- rnorm ( 1 , epsilon [ j ], sqrt ( var[ j ] ) )
        }
        afterdrift [[ i ]] <- drift
    }

    # now we calculate one Qx statistic, summed across the given number of populations. We use matrix method developed in "equation_verification.R"

    for ( i in 1:population ) {

      p <- afterdrift [[ i ]]
      A <- rep ( avgeff, size )
      D <- rep ( domdev, size )

      sum_matrix <- matrix ( ncol = size, nrow = size )

      for ( j in 1:size) {

        sum_matrix [ j , j ] <- ( (1/2) * A[ j ] * ( p [ j ] - epsilon [ j ] ) ) ^ 2 + ( A[ j ] * D [ j ] * ( 1 - 2 * p [ j ] ) * ( p [ j ] - epsilon [ j ] )^2 ) + ( D [ j ] * ( 1 - 2 * p [ j ] ) * ( p [ j ] - epsilon [ j ] ) ^ 2 )

      }

      sum1 <- sum ( diag ( sum_matrix ) )

      for ( j in 1:size ) {
        for ( k in 1:size ) {

          if ( j != k ) {
            sum_matrix [ j, k ] <- ( ( 1 / 4 ) * A [ j ] * A [ k ] * ( p [ j ] - epsilon[ j ] ) *  ( p [ k ] - epsilon [ k ] ) ) + (1/2) * A [ j ] * D [ k ] * ( 1 - 2 * p [ k ] ) * ( p [ j ] - epsilon [ j ] ) * ( p [ k ] - epsilon [ k ] ) + (1/2) *  D [ j ] * A [ k ] * ( 1 - 2 * p [ j ] ) * ( p [ j ] - epsilon [ j ] ) * ( p [ k ] - epsilon [ k ] ) + D [ j ] * ( 1 - 2 * p [ j ] ) * ( p [ j ] - epsilon [ j ] ) * D [ k ] * ( 1 - 2 * p [ k ] ) * ( p [ k ] - epsilon [ k ] )
          }
        }
      }

      sum2 <- sum ( sum_matrix [ row ( sum_matrix ) != col ( sum_matrix ) ] )

      sum_pops [ i ] <- sum1 + sum2

    }

    Qx [ h ]  <- sum ( sum_pops )

  }

  hist( Qx, main = paste ( "domdev = ", domdev, ", avgeff = ", avgeff, ", populations = ", population, ", reps = ", reps ))

  return ( Qx )

}



x <- domdev_distributions ( size = 100,
                       population = 5,
                       reps = 1000,
                       generations = 10,
                       epsilonrange1 = 0.4,
                       epsilonrange2 = 0.6,
                       avgeff = 0.5,
                       domdev = 0.3 )

mean (x)

hist( x, main = "domdev = 0.3, avgeff = 0.5, population = 1, size = 100, reps = 1000" )
