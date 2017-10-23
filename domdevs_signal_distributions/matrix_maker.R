# Function to split domdev signal expansion into terms, with a toggle to switch terms on and off to see how each term contributes to the overall effect. The plan is to integrate this function into domdev_distributions.

matrix_maker <- function ( term1 = TRUE, term2 = TRUE, term3 = TRUE, term4 = TRUE, term5 = TRUE, term6 = TRUE, size, A, D, p, epsilon ) {

    sum_matrix <- matrix ( rep(0, size^2 ), nrow= size, ncol=size )

    for ( j in 1:size ) {
        t1 <- ( (1/2) * A[ j ] * ( p [ j ] - epsilon [ j ] ) ) ^ 2
        t3 <- ( A[ j ] * D [ j ] * ( 1 - 2 * p [ j ] ) * ( p [ j ] - epsilon [ j ] )^2 )
        t5 <- ( D [ j ] * ( 1 - 2 * p [ j ] ) * ( p [ j ] - epsilon [ j ] ) )^2

            if ( term1 ) {
                sum_matrix [ j, j ] <- sum_matrix [ j, j ] + t1
            }

             if ( term3 ) {
                sum_matrix [ j, j ] <- sum_matrix [ j, j ] + t3
            }

             if ( term5 ) {
                sum_matrix [ j, j ] <- sum_matrix [ j, j ] + t5
            }
    }

    for ( j in 1:size) {
        for ( k in 1:size ) {
            t2 <- ( ( 1 / 4 ) * A [ j ] * A [ k ] * ( p [ j ] - epsilon[ j ] ) *  ( p [ k ] - epsilon [ k ] ) )
            t4 <- (1/2) * A [ j ] * D [ k ] * ( 1 - 2 * p [ k ] ) * ( p [ j ] - epsilon [ j ] ) * ( p [ k ] - epsilon [ k ] ) + (1/2) *  D [ j ] * A [ k ] * ( 1 - 2 * p [ j ] ) * ( p [ j ] - epsilon [ j ] ) * ( p [ k ] - epsilon [ k ] )
            t6 <- D [ j ] * ( 1 - 2 * p [ j ] ) * ( p [ j ] - epsilon [ j ] ) * D [ k ] * ( 1 - 2 * p [ k ] ) * ( p [ k ] - epsilon [ k ] )

            sum_matrix [ j, k ] <- t2 + t4 + t6

            if ( term2 ) {
                sum_matrix [ j, k ] <- sum_matrix [ j, k ] + t2
            }

             if ( term4 ) {
                sum_matrix [ j, k ] <- sum_matrix [ j, k ] + t4
            }

             if ( term6 ) {
                sum_matrix [ j, k ] <- sum_matrix [ j, k ] + t6
            }
        }
    }

    return ( sum_matrix )


}

matrix_maker ( size = 100,
              A = rep(0.5,100),
              D = rep(0,100),
              p = rep(0.5, 100),
              epsilon = rep(0.4, 100))


