library(jsonlite)

morgan2bpRate <- function(m, x0, breaks, rates, side=c("left","right")) {
  # turn breaks into Morgan
  segLen <- diff(breaks)
  mStart <- c(0, cumsum(rates * segLen))
  # position of the 1st SNP in Morgan
  i0 <- findInterval(x0, breaks, rightmost.closed = TRUE)
  i0 <- pmin(pmax(i0, 1), length(rates))
  mX0 <- mStart[i0] + rates[i0] * (x0 - breaks[i0])
  # recombination breakpoints in Morgan count from the 1st SNP
  M <- m + mX0
  i <- findInterval(M, mStart, rightmost.closed = TRUE)
  i <- pmin(pmax(i, 1), length(rates))
  # record zero-recombination-rate regions
  mEnd <- mStart[-1]
  plateau <- (rates[i] == 0) | (mEnd[i] == mStart[i])
  out <- numeric(length(M))

  # non-zero-recombination-rate regions
  ii <- which(!plateau)
  if (length(ii) > 0) {
    out[ii] <- breaks[i[ii]] + (M[ii] - mStart[i[ii]]) / rates[i[ii]]
  }

  # zero-recombination-rate regions
  jj <- which(plateau)
  if (length(jj) > 0) {
    out[jj] <- if (side == "left") breaks[i[jj]] else breaks[i[jj] + 1L]
  }

  out

}



recHistGenMatToSegDf <- function(histMat, x0, breaks, rates, seqLen) {

  origin <- as.integer(histMat[, 1])
  mStart <- as.numeric(histMat[, 2])
  mNext <- c(mStart[-1], NA_real_)

  left  <- morgan2bpRate(mStart, x0, breaks, rates, side="left")

  right <- numeric(length(mStart))
  if (length(mStart) > 1) {
    right[1:(length(mStart)-1)] <- morgan2bpRate(mNext[1:(length(mStart)-1)],
                                                 x0, breaks, rates, side="right")
  }
  right[length(mStart)] <- seqLen
  left[1] <- 0

  keep <- right > left

  data.frame(
    origin = origin[keep],
    left   = left[keep],
    right  = right[keep],
    stringsAsFactors = FALSE
  )
}

recHistGenToSegDfWithParents <- function(SP, offspringPop) {
  childIds <- offspringPop@id
  ped <- SP$pedigree[childIds, , drop = FALSE]

  out <- list()
  k <- 1

  for (childId in childIds) {
    x <- SP$recHistGen[[childId]]

    motherId <- ped[childId, "mother"]
    fatherId <- ped[childId, "father"]

    for (cc in seq_along(x)) {
      tc <- tc_load(chr_info[[cc]]$ts_path)
      seqLen <- as.numeric(tc$sequence_length())
      breaks <- chr_info[[cc]]$breaks
      rates <- chr_info[[cc]]$rates

      x0 <- chrKeptPosBpList[[cc]][1]

      haps <- as.vector(x[[cc]])
      nHap <- length(haps)

      for (h in seq_len(nHap)) {
        seg <- recHistGenMatToSegDf(haps[[h]], x0, breaks, rates, seqLen)

        parentId <- if (h <= nHap/2) motherId else fatherId

        seg$childId <- childId
        seg$chr <- cc
        seg$hap <- h
        seg$parentId <- parentId

        seg$parentHap <- seg$origin
        seg$parentGlobalHapId <- (parentId - 1) * nHap + seg$parentHap

        out[[k]] <- seg[, c("childId", "hap", "chr",
                            "left","right",
                            "parentId","parentHap","parentGlobalHapId")]
        k <- k + 1
      }
    }
  }

  do.call(rbind, out)
}

bridgeCollectSegGenFromSimOutput <- function(SP, simOutput) {
  bridgeSegDfListGen <<- list()

  for (k in 2:length(simOutput)) {
    segDf <- recHistGenToSegDfWithParents(SP, simOutput[[k]])
    bridgeSegDfListGen[[length(bridgeSegDfListGen) + 1]] <<- segDf
  }

  invisible(bridgeSegDfListGen)
}
