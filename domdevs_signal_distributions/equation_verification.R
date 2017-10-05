# creating parameter vectors

epsilon <- runif ( 100, 0.4, 0.6 )
drift <- 10/100
var <- drift * epsilon * ( 1 - epsilon )
size <- 100
afterdrift <- numeric()

for ( j in 1:size ) {
    afterdrift [ j ] <- rnorm ( 1 , epsilon [ j ], sqrt ( var[ j ] ) )
}

p <- afterdrift
A <- rep ( 0.5, 100 )
D <- rep ( 0.01, 100 )

# CHECKING EQUATION 6

alpha <- (1/2) * A + D * ( 1 - 2 * p )
sum_equ6 <- sum ( alpha * ( p - epsilon ) ) * sum ( alpha * ( p - epsilon ) )
print( sum_equ6 )


# CHECKING EQUATION 10


matrix2 <- matrix ( nrow = 100, ncol = 100 )

for ( i in 1:100 ) {
    matrix2 [ i , i ] <- ( (1/2) * A[ i ] * ( p [ i ] - epsilon [ i ] ) ) ^ 2
                      + ( A[ i ] * D [ i ] * ( 1 - 2 * p [ i ] ) * ( p [ i ] - epsilon [ i ] )^2 )
                      + ( D [ i ] * ( 1 - 2 * p [ i ] ) * ( p [ i ] - epsilon [ i ] ) )^ 2
}

sum1 <- sum ( diag ( matrix2 ) )

for ( i in 1:100 ) {
  for ( j in 1:100 ) {
    if ( i != j ) {
        matrix2 [ i , j ] <- ( ( 1 / 4 ) * A [ i ] * A [ j ] * ( p [ i ] - epsilon[ i ] ) *  ( p [ j ] - epsilon [ j ] ) )
        + (1/2) * A [ i ] * D [ j ] * ( 1 - 2 * p [ j ] ) * ( p [ i ] - epsilon [ i ] ) * ( p [ j ] - epsilon [ j ] ) + (1/2) *  D [ i ] * A [ j ] * ( 1 - 2 * p [ i ] ) * ( p [ i ] - epsilon [ i ] ) * ( p [ j ] - epsilon [ j ] )
        +  D [ i ] * ( 1 - 2 * p [ i ] ) * ( p [ i ] - epsilon [ i ] ) * D [ j ] * ( 1 - 2 * p [ j ] ) * ( p [ j ] - epsilon [ j ] )
    }
  }
}

sum2 <- sum ( matrix2 [ row ( matrix2 ) != col ( matrix2 ) ] )

sum_eq10 <- sum1 + sum2
print(sum_eq10)











