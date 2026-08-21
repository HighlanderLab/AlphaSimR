## Step 1 smoke test: read site positions via C API

library(AlphaSimRTmp)

tsPath <- "/Users/jliang2/Projects/test_TSK2ASR/data/simulations/normal/msprime_chr1.trees"

pos <- tsSitesPosition(tsPath)

cat("numSites =", length(pos), "\n")
cat("firstPositions =", paste(head(pos, 10), collapse = ", "), "\n")
cat("lastPositions  =", paste(tail(pos, 3), collapse = ", "), "\n")

library(reticulate)
use_virtualenv("~/r-reticulate-env", required = TRUE)
tskit <- import("tskit")
ts <- tskit$load(tsPath)
posPy <- ts$tables$sites$position
posPy <- as.numeric(py_to_r(ts$tables$sites$position))

cat("numSites =", length(posPy), "\n")
cat("firstPositions =", paste(head(posPy, 10), collapse = ", "), "\n")
cat("lastPositions  =", paste(tail(posPy, 3), collapse = ", "), "\n")
