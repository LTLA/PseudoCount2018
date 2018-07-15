# Simulates the effect of greater shrinkage on the effect size estimates.

all.dispersions <- c(0, 0.1, 1, 10) 
all.mu <- c(0.01, 0.1, 1, 10, 100)
N <- 1000

pdf("pics/power.pdf", width=12, height=10)
par(mfrow=c(length(all.dispersions), length(all.mu)),
    mar=c(4.1, 4.1, 3.1, 1.1))
legended <- FALSE    

for (disp in all.dispersions) { 
    if (disp==0) {
        FUN <- function(mu) rpois(N, lambda=mu)
    } else {
        FUN <- function(mu) rnbinom(N, mu=mu, size=1/disp)
    }

    for (mu in all.mu) {
        all.fc <- c(1, 1.5, 2, 4, 6, 8, 10)
        output.ref <- output.obs <- output.100 <- numeric(length(all.fc))
        for (i in seq_along(all.fc)) {
            A <- FUN(mu)
            B <- FUN(mu * all.fc[i])
            
            t1 <- t.test(log2(B + 1), log2(A + 1))
            t10 <- t.test(log2(B + 10), log2(A + 10))
            t100 <- t.test(log2(B + 100), log2(A + 100))
            output.ref[i] <- t1$statistic 
            output.obs[i] <- t10$statistic
            output.100[i] <- t100$statistic
        }

        ylim <- range(c(output.ref, output.obs, output.100))
        plot(all.fc, output.ref, col="black", ylim=ylim, xlab="Fold change", ylab="Welch t-statistic", 
            main=sprintf("Dispersion = %s, mean = %s", disp, mu), log="x", pch=4)
        lines(all.fc, output.ref, col="black")
        points(all.fc, output.obs, col="red", pch=16)
        lines(all.fc, output.obs, col="red")
        points(all.fc, output.100, col="blue", pch=15)
        lines(all.fc, output.100, col="blue")

        if (!legended) {
            legend("topleft", pch=c(4, 16, 15), lwd=1, col=c("black", "red", "blue"), legend=c(1, 10, 100), title="Pseudo-count")
            legended <- TRUE
        }
    }
}

dev.off()
