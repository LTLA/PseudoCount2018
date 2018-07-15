# Tests the theoretical error bounds.

source("functions.R")
all.sf <- c(0.01, 0.1, 1, 10, 100)
col <- viridis::viridis(length(all.sf))

set.seed(1000)
pdf("pics/error_bound.pdf")
for (disp in c(0, 0.1, 1, 10)) {
    output <- matrix(0, length(all.sf), length(all.sf))

    for (i1 in seq_along(all.sf)) {
        for (i2 in seq_len(i1-1)) {
            pseudo <- abs(1/all.sf[i1] - 1/all.sf[i2])
            output[i1, i2] <- FIND_MAX_ERR(s1=all.sf[i1], s2=all.sf[i2], pseudo=pseudo, size=1/disp)
        }
    }

    output <- output/log(2)
    PLOT_MAX_ERR(all.sf, output, col=col, legend=FALSE, ylim=c(0, 0.3))
    abline(h=1/8/log(2), col="red")
    legend("topright", fill=col, legend=all.sf, title="Second size factor", bg="white") 
}
dev.off()
