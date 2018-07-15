# Evaluating alternative transformations.

set.seed(99999)
ncells <- 500
lambda <- rep(c(1, 100), each=ncells/2)
ngenes <- 10000

size.fac <- lambda/mean(lambda)
PLOT_FUN <- function(g1, g2, main) {
    boxplot(list(Smaller=g1, Larger=g2), xlab="", ylab="Mean expression", main=main, col="grey80")
}

###################################
# Creating pictures of common transformations.

pdf("pics/alternatives.pdf", width=10, height=4)
par(mfrow=c(1,4), mar=c(4.1, 4.1, 4.1, 1.1), cex.lab=1.2)

for (scenario in 1:2) {
    if (scenario==1L) {
        X <- matrix(rpois(ngenes*ncells, lambda=lambda), ncol=ncells, byrow=TRUE)
    } else {
        X <- matrix(rnbinom(ngenes*ncells, mu=lambda, size=1), ncol=ncells, byrow=TRUE)
    }

    # Square rooting.
    rooted <- sqrt( t(t(X)/size.fac) )
    g1 <- rowMeans(rooted[,lambda==1])
    g2 <- rowMeans(rooted[,lambda==100])
    PLOT_FUN(g1, g2, main="Square root")
        
    # Log-transform 
    log1p <- log2( t(t(X)/size.fac) + 1)
    g1 <- rowMeans(log1p[,lambda==1])
    g2 <- rowMeans(log1p[,lambda==100])
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
}

dev.off()

###################################
# Examining SCnorm's behaviour.

library(SCnorm)
ncells <- 500
sf <- rep(c(1, 100), each=ncells/2)
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
