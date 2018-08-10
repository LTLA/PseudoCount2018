# This script evaluates alternative transformations.

set.seed(99999)
ncells <- 500
ngenes <- 10000
abundance <- 2^runif(ngenes, -2, 10)

# Specifying the size factor differences between cells.
MULT <- rep(c(0.1, 10), each=ncells/2)
first.group <- MULT==0.1
second.group <- !first.group

size.fac <- MULT 
size.fac <- size.fac/mean(size.fac)

# Central plotting function.
par(cex.lab=1.4, cex.main=1.5)
PLOT_FUN <- function(g1, g2, ...) {
    Y <- g1-g2
    X <- log10(abundance)
    smoothScatter(X, Y, xlab=expression(Log[10]~"mean abundance"), ylab="Difference in group means", ...)
    fit <- loess(Y ~ X, span=0.1)
    curve(predict(fit, data.frame(X=x)), add=TRUE, col="red", lwd=2)
}

###################################
# Creating pictures of common transformations.

for (scenario in 1:2) {
    if (scenario==1L) {
        lambda <- outer(abundance, size.fac, "*")
        X <- matrix(rpois(ngenes*ncells, lambda=lambda), ncol=ncells)
        ylim <- c(-3, 0.5)
        pdf("pics/alternative_pois.pdf", width=4, height=6)
    } else {
        lambda <- outer(abundance, size.fac, "*")
        X <- matrix(rnbinom(ngenes*ncells, mu=lambda, size=1), ncol=ncells)
        ylim <- c(-4, 0.5)
        pdf("pics/alternative_nb.pdf", width=4, height=6)
    }

    # Square rooting.
    rooted <- sqrt( t(t(X)/size.fac) )
    g1 <- rowMeans(rooted[,first.group])
    g2 <- rowMeans(rooted[,second.group])
    PLOT_FUN(g1, g2, main="Square root", ylim=ylim)
        
    # Log-transform 
    log1p <- log2( t(t(X)/size.fac) + 1 )
    g1 <- rowMeans(log1p[,first.group])
    g2 <- rowMeans(log1p[,second.group])
    PLOT_FUN(g1, g2, main="Log", ylim=ylim)
    
    # VST.
    library(DESeq2)
    dds <- DESeqDataSetFromMatrix(X, colData=DataFrame(row.names=seq_len(ncells)), ~1)
    sizeFactors(dds) <- size.fac
    vst.out <- vst(dds, fitType="mean")
    
    g1 <- rowMeans(assay(vst.out)[,first.group])
    g2 <- rowMeans(assay(vst.out)[,second.group])
    PLOT_FUN(g1, g2, main="DESeq2 VST", ylim=ylim)

    # Our proposed approach.
    pseudo <- log2( t(t(X)/size.fac) + abs(diff(1/range(size.fac))) )
    g1 <- rowMeans(pseudo[,first.group])
    g2 <- rowMeans(pseudo[,second.group])
    PLOT_FUN(g1, g2, main="Increased pseudo-count", ylim=ylim)

    dev.off()
}

###################################
# Examining SCnorm's behaviour.

library(SCnorm)
ncells <- 500
sf <- rep(c(0.1, 10), each=ncells/2)
ngenes <- 10000
abundances <- runif(ngenes, 0, 10)

X <- matrix(rpois(ngenes*ncells, lambda=outer(abundances, sf)), ncol=ncells)
rownames(X) <- paste0("GENE_", seq_len(ngenes))
colnames(X) <- paste0("CELL_", seq_len(ncells))

normed <- SCnorm(X, Conditions=rep(1, ncells))
g1 <- rowMeans(log2(metadata(normed)$NormalizedData[,sf==1] + 1))
g2 <- rowMeans(log2(metadata(normed)$NormalizedData[,sf==100] + 1))

pdf("scnorm.pdf")
plot(abundances, g1-g2, xlab="True abundance", ylab="Log2-fold change")
dev.off()
