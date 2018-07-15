# Evaluating alternative transformations.

set.seed(99999)
ncells <- 500
lambda <- rep(c(1, 100), each=ncells/2)
ngenes <- 10000
X <- matrix(rpois(ngenes*ncells, lambda=lambda), ncol=ncells, byrow=TRUE)

size.fac <- lambda/mean(lambda)
PLOT_FUN <- function(g1, g2, main) {
    boxplot(list(Smaller=g1, Larger=g2), xlab="", ylab="Mean expression", main=main, col="grey80")
}

pdf("pics/alternatives.pdf")

# Square rooting.
rooted <- sqrt( t(t(X)/size.fac) )
g1 <- rowMeans(rooted[,lambda==1])
g2 <- rowMeans(rooted[,lambda==100])
PLOT_FUN(g1, g2, main="Square root")
    
# Log-transform 
arcsinh <- log2( t(t(X)/size.fac) + 1)
g1 <- rowMeans(arcsinh[,lambda==1])
g2 <- rowMeans(arcsinh[,lambda==100])
PLOT_FUN(g1, g2, main="Log2(X + 1)")

# VST.
library(DESeq2)
dds <- DESeqDataSetFromMatrix(X, colData=DataFrame(row.names=seq_len(ncells)), ~1)
sizeFactors(dds) <- size.fac
vst.out <- vst(dds, fitType="mean")

g1 <- rowMeans(assay(vst.out)[,lambda==1])
g2 <- rowMeans(assay(vst.out)[,lambda==100])
PLOT_FUN(g1, g2, main="VST")

# Our proposed approach.
pseudo <- log2( t(t(X)/size.fac) + abs(diff(1/range(size.fac))) )
g1 <- rowMeans(pseudo[,lambda==1])
g2 <- rowMeans(pseudo[,lambda==100])
PLOT_FUN(g1, g2, "Increased pseudo-count")
dev.off()


