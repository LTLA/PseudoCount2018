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

# Creating a plot.
all.extra <- all.lfcs - all.ref
all.extra.big <- all.big - all.ref
combined <- c(all.extra, all.extra.big)
breaks <- seq(min(combined), max(combined), length.out=20)

pdf("pics/ercc.pdf")
hist(all.extra, xlab="Error", breaks=breaks, col="grey80", main="Error with pseudo-count = 1")
hist(all.extra.big, xlab="Error", breaks=breaks, col="grey80", main="Error with large pseudo-count")

all.lfcs <- c(all.ref, all.lfcs, all.big)
fc.breaks <- seq(min(all.lfcs), max(all.lfcs), length.out=20)
hist(all.ref, xlab="Difference in log-means", breaks=fc.breaks, col="grey80", main="")
hist(all.lfcs, xlab="Difference in mean-logs (+ 1)", breaks=fc.breaks, col="grey80", main="")
hist(all.big, xlab="Difference in mean-logs (+ large)", breaks=fc.breaks, col="grey80", main="")
dev.off()

