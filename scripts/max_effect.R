# Determines the maximum error of the log-transformation under a variety of conditions.

source("functions.R")
all.sf <- c(0.01, 0.1, 1, 10, 100)
col <- viridis::viridis(length(all.sf))

set.seed(1000)
pdf("pics/max_effect.pdf")
for (disp in c(0, 0.1, 1, 10)) {
    output <- matrix(0, length(all.sf), length(all.sf))

    for (i1 in seq_along(all.sf)) {
        for (i2 in seq_len(i1)) {
            output[i1, i2] <- FIND_MAX_ERR(s1=all.sf[i1], s2=all.sf[i2], size=1/disp)
        }
    }

    output <- output/log(2)
    PLOT_MAX_ERR(all.sf, output, col=col)
}
dev.off()

