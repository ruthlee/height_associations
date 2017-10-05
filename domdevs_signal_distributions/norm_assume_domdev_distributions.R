# This function will make use of the MultiVarChiSquared function (which uses the normal approximation of drift) to evaluate what happens to the average effect when it depends on the homozygous effect and the dominance deviation.

equ22_Test <- function( size, population, reps, generations, epsilonrange1, epsilonrange2, avgeff, domdev ) {
                                        # First figure out what the ancestral mean and variance would be, given a normal distribution of the ancestral frequency.

    epsilon <- runif ( size, epsilonrange1, epsilonrange2 )
    alpha <- numeric()
    drift <- generations / ( size )
    afterdrift <- list ()
    chisquared <- numeric ()

    if ( size <= generations ) {
        stop ( "Population size should be greater than generations for accurate results.")
     }

    mean <- epsilon
    var <- drift * epsilon * ( 1 - epsilon )

    for ( h in 1:reps) {

        for ( i in 1:population ) {
            afterdrift [[ i ]] <- rnorm ( size , mean , sqrt ( var ) )
            alpha <- ( 1/2 ) * avgeff + domdev * ( 1 - 2 * afterdrift [[ i ]] )
        }

        # given the list of frequencies after drift (each element of the list is a vector of size "size", we can now calculate a vector of alphas from the average effects and dominance deviations specified above.

        z <- numeric ()
        estepsilon <- numeric ()

        for ( j in 1:population ) {
            z [ j ] <- sum ( alpha * afterdrift [[ j ]] )
            }
                                        # put in a switch -- if i specify only one population, automatically just use ancestral value.

        multimean <- mean ( z )

        for ( k in 1:length( afterdrift ) ) {
            estepsilon [ k ] <- mean ( afterdrift [[ k ]] )
            }

        multivar <- ( sum ( alpha^2 * estepsilon * ( 1 - estepsilon ) ) )

        chisquared [ h ] <- sum ( ( z - multimean ) ^ 2 / ( multivar * drift ) )
    }

    hist ( chisquared, breaks = 50, main = print ( paste ( "Eq 22 Histogram, avgeff = ", avgeff, ", domdev = ", domdev ) ) )
    return (chisquared)

}

x <- equ22_Test ( size = 1000,
             population = 2,
             reps = 10000,
             generations = 100,
             epsilonrange1 = 0.4,
             epsilonrange2 = 0.6,
             avgeff = 0.5,
             domdev = 0 )








