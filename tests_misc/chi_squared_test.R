# Now we are simulating a normal distribution taking frequencies and effect sizes from
# the data as "ancestral" frequncies and simulating a bunch of random normally distributed
# numbers to plot.

epsilon <- 0.5
alpha <- 0.1

mean <- sum ( alpha * epsilon )
variance <- sum ( alpha^2 * epsilon * ( 1 - epsilon ) )

x <- rnorm ( 1000000, mean, variance )

hist( x, 100 )

x2 <- ( ( x - mean ) / ( variance ) ) ^ 2

hist ( x2, 100 , main = "chi-squared with 1 deg freedom" )
