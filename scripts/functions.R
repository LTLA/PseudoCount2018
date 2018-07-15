#' Simulates NB-distributed counts and computes the difference in means of log-values.
SIMULATE_ERR <- function(s1, s2, mu, pseudo=1, size=1000, N=10000) {
    if (is.infinite(size)) {
        FUN <- function(mu) rpois(N, lambda=mu)
    } else {
        FUN <- function(mu) rnbinom(N, mu=mu, size=size)
    }

    output <- numeric(length(mu))
    for (i in seq_along(mu)) {
        A <- FUN(mu[i] * s1)
        B <- FUN(mu[i] * s2)
        output[i] <- mean(log(A/s1 + pseudo)) - mean(log(B/s2 + pseudo))
    }
    return(output)
}

dir.create("pics", showWarnings=FALSE)
