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

#' Creating a plot of the maximum errors.
PLOT_MAX_ERR <- function(sf, mat, col, legend=TRUE) {
    plot(sf, numeric(length(sf)), ylim=range(mat, na.rm=TRUE), type="n", log="x",
        xlab="First size factor", ylab=expression(Log[2]-"fold change"),
        main=sprintf("Dispersion of %s", disp))

    for (i in seq_len(nrow(mat))) {
        current <- seq_len(i)
        lines(sf[current], mat[i,current], lwd=2)
        lines(sf[current], mat[i,current], lwd=1.5, col=col[i])
        points(sf[current], mat[i,current], pch=21, bg=col[i], cex=1.5)
    }

	if (legend) {
        legend("topright", fill=col, legend=sf, title="Second size factor") 
    }
    return(NULL)
}

dir.create("pics", showWarnings=FALSE)
