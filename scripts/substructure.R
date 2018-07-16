
#' Computes a larger pseudo-count based on the size factors.
INFLATE_PSEUDO <- function(sf, p=0.05) {
    max(1, abs(diff(1/quantile(sf, c(p, 1-p)))))
}

#' Makes a PCA plot.
PLOT_PCA <- function(mat, ...) {
    out <- prcomp(mat)
    total.var <- sum(out$sdev^2)
    plot(out$x[,1], out$x[,2], ...,
        xlab=sprintf("PC1 (%.1f%%)", out$sdev[1]^2/total.var * 100),
        ylab=sprintf("PC2 (%.1f%%)", out$sdev[2]^2/total.var * 100)
    )
}

# Simulations with two separate clusters of observations.

set.seed(99999)
ncells <- 500
lambda <- rep(c(0.1, 10), each=ncells/2)
ngenes <- 1000

X <- matrix(rpois(ngenes*ncells, lambda=lambda), ncol=ncells, byrow=TRUE)
size.fac <- lambda/mean(lambda)    
Z <- log2(t(X)/size.fac + 1)

pdf("pics/clusters.pdf")
PLOT_PCA(Z, bg=ifelse(lambda > 1, "salmon", "dodgerblue"), pch=21)
legend("topright", pt.bg=c("salmon", "dodgerblue"), pch=21, legend=c("Large", "Small"))
dev.off()

# Simulations with a log-transformation-induced trajectory.

set.seed(99999)
ncells <- 500
lambda <- runif(ncells, 0.1, 5)
ngenes <- 1000

X <- matrix(rpois(ngenes*ncells, lambda=lambda), ncol=ncells, byrow=TRUE)
size.fac <- lambda/mean(lambda)    
Z <- log2(t(X)/size.fac + 1)

all.col <- viridis::magma(20) 
col <- all.col[cut(log(size.fac), 20)]

pdf("pics/trajectory.pdf")    
layout(cbind(1,2), width=c(6, 1))
par(mar=c(5.1, 4.1, 4.1,0.1))
PLOT_PCA(Z, bg=col, pch=21)

par(mar=c(0.1, 0.1, 0.1,0.1))
plot(0,0, type="n", axes=FALSE, ylab="", xlab="", ylim=c(-1, 1), xlim=c(-1, 0.5))
start.loc <- seq(-0.5, 0.5, length.out=length(all.col))
interval <- diff(start.loc)[1]
rect(-0.3, start.loc, 0.3, start.loc+interval, col=all.col, border=NA)
text(0, -0.5, pos=1, "Low", cex=1.5)
text(0, 0.5+interval,  pos=3, "High", cex=1.5)
dev.off()
