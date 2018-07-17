# Checks for log-fold change distortions in public ERCC data.

library(BiocFileCache)
bfc <- BiocFileCache(ask=FALSE)
fname <- bfcrpath(bfc, "http://cf.10xgenomics.com/samples/cell-exp/1.1.0/ercc/ercc_raw_gene_bc_matrices.tar.gz")

tempdir <- tempfile()
dir.create(tempdir)
untar(fname, exdir=tempdir)

library(DropletUtils)
sce <- read10xCounts(file.path(tempdir, "matrices_mex/ercc92"))

library(Matrix)
lib.sizes <- colSums(counts(sce))
keep <- lib.sizes > 100

sce <- sce[,keep]
lib.sizes <- lib.sizes[keep]
size.facs <- lib.sizes/mean(lib.sizes)
normed <- t(t(counts(sce))/size.facs)

# Computing log-fold changes between the lowest and highest subsets.
all.prop <- 0.2
lower <- lib.sizes < quantile(lib.sizes, all.prop)
upper <- lib.sizes > quantile(lib.sizes, 1-all.prop)

left <- normed[,lower]
right <- normed[,upper]

all.lfcs <- rowMeans(log2(right+1)) - rowMeans(log2(left+1))
all.ref <- log2((rowMeans(right)+1)/(rowMeans(left)+1))

big.pseudo <- abs(1/median(size.facs[lower]) - 1/median(size.facs[upper]))
all.big <- rowMeans(log2(right+big.pseudo)) - rowMeans(log2(left+big.pseudo))

ave.count <- rowMeans(normed)

# Creating a plot.
ylim <- range(c(all.big, all.ref, all.lfcs))

pdf("pics/ercc.pdf")
plot(ave.count, all.lfcs, xlab="Average count", log="x", ylab="Difference in log-values", pch=16, ylim=ylim)
plot(ave.count, all.ref, xlab="Average count", log="x", ylab="Log-fold change in mean counts", pch=16, ylim=ylim)
plot(ave.count, all.big, xlab="Average count", log="x", ylab="Difference in log-values", pch=16, ylim=ylim)
dev.off()

