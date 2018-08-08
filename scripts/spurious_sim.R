#########################################
# Defining some common functions.

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

#' Finds the maximum error, which should occur around the pseudo-count.
FIND_MAX_ERR <- function(s1, s2, pseudo=1, size=1000, N=10000) {
    mus <- 10^seq(-3, 3, length.out=100) * pseudo
    output <- SIMULATE_ERR(s1, s2, mu=mus, pseudo=pseudo, size=size, N=N)
    max(abs(output))
}

dir.create("pics", showWarnings=FALSE)
all.sf <- c(0.01, 0.1, 1, 10, 100)
col <- viridis::viridis(length(all.sf))
dispersions <- c(0, 0.1, 1, 10)

#########################################
# Determines the maximum error of the log-transformation under a variety of conditions.

set.seed(1000)
pdf("pics/max_effect.pdf")
for (disp in dispersions) { 

    collected <- vector("list", 10)
    for (it in seq_along(collected)) {
        output <- matrix(0, length(all.sf), length(all.sf))
        for (i1 in seq_along(all.sf)) {
            for (i2 in seq_len(i1)) {
                output[i1, i2] <- FIND_MAX_ERR(s1=all.sf[i1], s2=all.sf[i2], size=1/disp)
            }
        }
        collected[[it]] <- output
    }

    collected <- do.call(abind::abind, c(collected, list(along=3)))
    output.mean <- output.se <- matrix(0, length(all.sf), length(all.sf))
    for (i in seq_along(all.sf)) {
        for (j in seq_len(i)) {
            current <- collected[i,j,]
            output.mean[i,j] <- mean(current)
            output.se[i,j] <- sd(current)/sqrt(length(current))
        }
    }

    output.mean <- output.mean/log(2)
    output.se <- output.se/log(2)
    upper <- output.mean + output.se

    plot(all.sf, numeric(length(all.sf)), ylim=range(upper), type="n", log="x",
        xlab="First size factor", ylab=expression("Maximum difference in"~log[2]*"-expression"),
        main=sprintf("Dispersion of %s", disp))

    for (i in seq_len(nrow(output))) {
        current <- seq_len(i)
        cur.mean <- output.mean[i,current]
        lines(all.sf[current], cur.mean, lwd=2)
        lines(all.sf[current], cur.mean, lwd=1.5, col=col[i])
        points(all.sf[current], cur.mean, pch=21, bg=col[i], cex=1.5)
    }
    legend("topright", fill=col, legend=all.sf, title="Second size factor") 
}
dev.off()

#########################################
# Tests the theoretical error bounds.

set.seed(1000)
pdf("pics/error_bound.pdf")
for (disp in dispersions) { 
    output <- matrix(0, length(all.sf), length(all.sf))

    for (i1 in seq_along(all.sf)) {
        for (i2 in seq_len(i1-1)) {
            pseudo <- abs(1/all.sf[i1] - 1/all.sf[i2])
            output[i1, i2] <- FIND_MAX_ERR(s1=all.sf[i1], s2=all.sf[i2], pseudo=pseudo, size=1/disp)
        }
    }

    output <- output/log(2)
    output <- Matrix::forceSymmetric(output, "L")
    diag(output) <- NA
    plot(all.sf, numeric(length(all.sf)), ylim=c(0, 0.3), type="n", log="x",
        xlab="First size factor", ylab=expression("Maximum difference in"~log[2]*"-expression"),
        main=sprintf("Dispersion of %s", disp))

    for (i in seq_len(nrow(output))) {
        current <- output[i,]
        lines(all.sf, current, lwd=2)
        lines(all.sf, current, lwd=1.5, col=col[i])

        if (i!=1L && i!=length(all.sf)) {
            gap <- c(i-1, i+1)
            lines(all.sf[gap], current[gap], col=col[i], lty=2)
        }

        points(all.sf, current, pch=21, bg=col[i], cex=1.5)
    }

    abline(h=1/8/log(2), col="red")
    legend("topright", fill=col, legend=all.sf, title="Second size factor") 
}
dev.off()
