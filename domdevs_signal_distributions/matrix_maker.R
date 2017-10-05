# Function to split domdev signal expansion into terms, with a toggle to switch terms on and off to see how each term contributes to the overall effect. The plan is to integrate this function into domdev_distributions.

matrix_maker <- function ( term1 = TRUE, term2 = TRUE, term3 = TRUE, term4 = TRUE, term5 = TRUE, term6 = TRUE, size, A, D, p, epsilon ) {

    sum_matrix <- matrix (  ncol = size, nrow = size )

    for ( j in 1:size ) {
        term1 <- ( (1/2) * A[ j ] * ( p [ j ] - epsilon [ j ] ) ) ^ 2
        term3 <- ( A[ j ] * D [ j ] * ( 1 - 2 * p [ j ] ) * ( p [ j ] - epsilon [ j ] )^2 )
        term5 <- ( D [ j ] * ( 1 - 2 * p [ j ] ) * ( p [ j ] - epsilon [ j ] ) ^ 2 )

        sum_matrix [ j, j ] <- term1 + term3 + term5

            if ( term1 == FALSE ) {
                sum_matrix [ j, k ] <- sum_matrix [ j, k ] - term2
            }

             if ( term3 == FALSE ) {
                sum_matrix [ j, k ] <- sum_matrix [ j, k ] - term4
            }

             if ( term5 == FALSE ) {
                sum_matrix [ j, k ] <- sum_matrix [ j, k ] - term6
            }
    }

    for ( j in 1:size) {
        for ( k in 1:size ) {
            term2 <- ( ( 1 / 4 ) * A [ j ] * A [ k ] * ( p [ j ] - epsilon[ j ] ) *  ( p [ k ] - epsilon [ k ] ) )
            term4 <- (1/2) * A [ j ] * D [ k ] * ( 1 - 2 * p [ k ] ) * ( p [ j ] - epsilon [ j ] ) * ( p [ k ] - epsilon [ k ] ) + (1/2) *  D [ j ] * A [ k ] * ( 1 - 2 * p [ j ] ) * ( p [ j ] - epsilon [ j ] ) * ( p [ k ] - epsilon [ k ] )
            term6 <- D [ j ] * ( 1 - 2 * p [ j ] ) * ( p [ j ] - epsilon [ j ] ) * D [ k ] * ( 1 - 2 * p [ k ] ) * ( p [ k ] - epsilon [ k ] )

            sum_matrix [ j, k ] <- term2 + term4 + term6

            if ( term2 == FALSE ) {
                sum_matrix [ j, k ] <- sum_matrix [ j, k ] - term2
            }

             if ( term4 == FALSE ) {
                sum_matrix [ j, k ] <- sum_matrix [ j, k ] - term4
            }

             if ( term6 == FALSE ) {
                sum_matrix [ j, k ] <- sum_matrix [ j, k ] - term6
            }
        }
    }

    return ( sum_matrix )


}


