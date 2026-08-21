library(reticulate)
library(jsonlite)

use_virtualenv("~/r-reticulate-env", required = TRUE)
tskit <- import("tskit")

genLogRecord <- function(parentPop, offspringPop, simParam, genIndex) {
  list(
    genIndex = genIndex,
    offspringIds = offspringPop@id,
    ibd = pullIbdHaplo(offspringPop, simParam = simParam)
  )
}

ibdToSegDf <- function(ibdMat) {
  # chr_locus
  cn <- colnames(ibdMat)
  # ind_hap
  rn <- rownames(ibdMat)

  chr <- as.integer(sub("_.*$", "", cn))
  locus <- as.integer(sub("^.*_", "", cn))

  childId <- sub("_.*$", "", rn)
  hap <- as.integer(sub("^.*_", "", rn))

  out <- list()
  k <- 1

  num_hap <- length(hap)
  uniqueChr <- sort(unique(chr))

  for (r in seq_len(num_hap)) {
    # each row in ibdMat = every child hap
    v <- ibdMat[r, ]
    for (cc in uniqueChr) {
      # get index from chr to extract parent hap (vv) and position (ll)
      idx <- which(chr == cc)
      vv <- v[idx]
      ll <- locus[idx]

      # find breakpoints
      chg <- which(diff(vv) != 0)
      starts <- c(1, chg + 1)
      ends <- c(chg, length(vv))

      out[[k]] <- data.frame(
        childId = childId[r],
        hap = hap[r],
        chr = cc,
        # extract position of ibd change
        locusStart = ll[starts],
        locusEnd = ll[ends],
        origin = vv[starts],
        stringsAsFactors = FALSE
      )
      k <- k + 1
    }
  }
  do.call(rbind, out)
}

recHistMatToSegDfPy <- function(histMat, nLoci) {

  origin <- as.integer(histMat[, 1])
  starts <- as.integer(histMat[, 2])

  ends <- c(starts[-1] - 1L, nLoci)

  data.frame(
    origin = origin,
    locusStart = starts,
    locusEnd = ends,
    stringsAsFactors = FALSE
  )
}


recHistToSegDfWithParentsPy <- function(SP, offspringPop, nLociByChr) {
  childIds <- offspringPop@id
  ped <- SP$pedigree[childIds, , drop = FALSE]

  out <- list()
  k <- 1

  for (childId in childIds) {
    x <- SP$recHist[[childId]]

    motherId <- ped[childId, "mother"]
    fatherId <- ped[childId, "father"]

    for (cc in seq_along(x)) {
      nLoci <- nLociByChr[[cc]]

      haps <- as.vector(x[[cc]])
      nHap <- length(haps)

      for (h in seq_len(nHap)) {
        seg <- recHistMatToSegDfPy(haps[[h]], nLoci = nLoci)

        parentId <- if (h <= nHap/2) motherId else fatherId

        seg$childId <- childId
        seg$chr <- cc
        seg$hap <- h
        seg$parentId <- parentId

        seg$parentHap <- seg$origin
        seg$parentGlobalHapId <- (parentId - 1) * nHap + seg$parentHap

        out[[k]] <- seg[, c("childId", "hap", "chr",
                            "locusStart","locusEnd",
                            "parentId","parentHap","parentGlobalHapId")]
        k <- k + 1
      }
    }
  }

  do.call(rbind, out)
}

bridgeCollectSegFromSimOutputPy <- function(SP, simOutput) {
  bridgeSegDfList <<- list()

  nLociByChr <- lapply(chrKeptPosBpList, length)

  for (k in 2:length(simOutput)) {
    segDf <- recHistToSegDfWithParentsPy(SP, simOutput[[k]], nLociByChr)
    bridgeSegDfList[[length(bridgeSegDfList) + 1]] <<- segDf
  }

  invisible(bridgeSegDfList)
}


segDfToEdgeDfUsingBridge <- function(segDf, chr_info) {
  # segDF: childID, hap, chr, locusStart, locusEnd, origin
  out <- segDf
  out$left <- NA_real_
  out$right <- NA_real_

  for (cc in sort(unique(out$chr))) {
    #posBp <- bridgeEnv$chrKeptPosBpList[[cc]]
    #if (is.null(posBp)) stop("bridgeEnv$chrKeptPosBpList[[", cc, "]] is NULL.")
    posBp <- chrKeptPosBpList[[cc]]

    tsPath <- chr_info[[cc]]$ts_path
    ts <- tskit$load(tsPath)
    seqLen <- as.numeric(ts$sequence_length)

    idx <- which(out$chr == cc)

    for (i in idx) {
      s <- out$locusStart[i]
      e <- out$locusEnd[i]
      out$left[i] <- if (s == 1) 0 else posBp[s]
      out$right[i] <- if (e < length(posBp)) posBp[e + 1] else seqLen
    }
  }
  out
}

bridgeAllSegToEdgeDfPy <- function(chr_info) {
  allSeg <- do.call(rbind, bridgeSegDfList)

  out <- allSeg
  out$left <- NA
  out$right <- NA

  for (cc in sort(unique(out$chr))) {
    posBp <- chrKeptPosBpList[[cc]]

    tsPath <- chr_info[[cc]]$ts_path
    ts <- tskit$load(tsPath)
    seqLen <- ts$sequence_length

    idx <- which(out$chr == cc)
    for (i in idx) {
      s <- out$locusStart[i]
      e <- out$locusEnd[i]
      out$left[i]  <- if (s == 1) 0 else posBp[s]
      out$right[i] <- if (e < length(posBp)) posBp[e + 1] else seqLen
    }
  }

  out
}

bridgeComputeIndTimePy <- function(pedigree) {
  n <- nrow(pedigree)
  indTime <- rep(NA, n)

  for (i in 1:n) {
    m <- pedigree[i, "mother"]
    f <- pedigree[i, "father"]

    if (m == 0 && f == 0) {
      indTime[i] <- 0
    } else {
      indTime[i] <- min(indTime[m], indTime[f]) - 1
    }
  }

  indTime
}


bridgeWriteTreesPy <- function(chr_info, edgeDf, SP, out_dir = NULL,
                             out_basename = "AlphaSimR_extended") {
  nodeSchema <- tskit$MetadataSchema(list(
    codec = "json",
    type = "object",
    properties = list(
      alphaSimR = list(
        type = "object",
        properties = list(
          id = list(type = "string", description = "AlphaSimR node id (childId_hap)")
        ),
        required = list("id"),
        additionalProperties = FALSE
      )
    ),
    required = list("alphaSimR"),
    additionalProperties = FALSE
  ))

  indTime <- bridgeComputeIndTimePy(SP$pedigree)

  nodeIdMapByChr <<- vector("list", length(chr_info))
  indIdMapByChr  <<- vector("list", length(chr_info))

  for (cc in seq_along(chr_info)) {
    #nodeIdMap <<- list()
    nodeIdMapByChr[[cc]] <<- list()
    indIdMapByChr[[cc]]  <<- list()

    ts <- tskit$load(chr_info[[cc]]$ts_path)
    tables <- ts$dump_tables()
    reticulate::py_set_attr(tables$nodes, "metadata_schema", nodeSchema)

    # for metadata
    n <- tables$nodes$num_rows
    encoded <- vector("list", n)
    for (i in 0:(n - 1)) {
      md_i <- tables$nodes[i]$metadata
      encoded[[i + 1]] <- tables$nodes$metadata_schema$encode_row(md_i)
    }
    # ----


    df <- edgeDf[edgeDf$chr == cc, , drop = FALSE]
    if (nrow(df) == 0) next

    # get indIDs for sampled nodes
    sampNodeId <- ts$samples()
    sampIndRow <- integer(length(sampNodeId))
    for (i in seq_along(sampNodeId)) {
      sampIndRow[i] <- tables$nodes[sampNodeId[i]]$individual
    }
    if (any(sampIndRow < 0)) {
      bad <- which(sampIndRow < 0)[1]
      stop(
        "Sample node", sampNodeId[bad], "has individual = -1. ",
        "Cannot reuse founders' individuals. ",
      )
    }

    nFounder <- length(sampNodeId) / ploidy
    idx <- 1
    for (ind in 1:nFounder) {
      indRow <- sampIndRow[idx]
      indIdMapByChr[[cc]][[as.character(ind)]] <<- indRow

      for (h in 1:ploidy) {
        #nodeIdMap[[paste(ind, h, sep = "_")]] <<- as.integer(sampNodeId[[idx]])
        nodeId <- as.integer(unlist(sampNodeId[[idx]]))[1]
        key <- paste(ind, h, sep = "_")
        nodeIdMapByChr[[cc]][[key]] <<- nodeId
        #tables$nodes$metadata[[nodeId + 1]] <- tskit$pack_bytes(list(alphaSimR = list(id = key)))
        #tables$nodes[nodeId + 1] <- tables$nodes[nodeId + 1]$replace(metadata=list(alphaSimR = list(id = key)))
        encoded[[nodeId + 1]] <- tables$nodes$metadata_schema$encode_row(
          list(alphaSimR = list(id = key)))
        idx <- idx + 1
      }
    }
    tables$nodes$packset_metadata(encoded)

    # add indIDs for offSpring nodes
    nextInd <- as.integer(tables$individuals$num_rows)
    addNewIndividual <- function(alphaId) {
      key <- as.character(alphaId)
      if (!is.null(indIdMapByChr[[cc]][[key]])) return(indIdMapByChr[[cc]][[key]])

      m <- SP$pedigree[alphaId, "mother"]
      f <- SP$pedigree[alphaId, "father"]

      mRow <- addNewIndividual(m)
      fRow <- addNewIndividual(f)

      newId <- nextInd
      tables$individuals$add_row(
        parents = list(as.integer(mRow), as.integer(fRow)),
        metadata = list(file_id=as.integer(newId)))
      indIdMapByChr[[cc]][[key]] <<- as.integer(newId)

      nextInd <<- nextInd + 1L
      newId
    }

    childIdsNeeded <- sort(as.integer(unique(df$childId)))
    for (childId in childIdsNeeded) {
      addNewIndividual(childId)
    }

    # append child nodes
    childKeys <- unique(paste(df$childId, df$hap, sep = "_"))
    for (key in childKeys) {
      #if (is.null(nodeIdMap[[key]])) {
      if (is.null(nodeIdMapByChr[[cc]][[key]])) {
        childId <- as.integer(sub("_.*$", "", key))
        indRow  <- indIdMapByChr[[cc]][[as.character(childId)]]

        tables$nodes$add_row(
          flags = 0L,
          time  = indTime[[childId]],
          population = -1L,
          individual = indRow,
          metadata = list(alphaSimR = list(id = key))
        )
        #nodeIdMap[[key]] <<- as.integer(tables$nodes$num_rows - 1)
        nodeIdMapByChr[[cc]][[key]] <<- as.integer(tables$nodes$num_rows - 1)
      }
    }

    # append edges
    for (i in 1:nrow(df)) {
      parentKey <- paste(df$parentId[i], df$parentHap[i], sep = "_")
      childKey  <- paste(df$childId[i], df$hap[i], sep = "_")

      #if (is.null(nodeIdMap[[key]])) {
      if (is.null(nodeIdMapByChr[[cc]][[parentKey]])) {
        stop("Missing parent node for key=", parentKey,
             " on chr=", cc, ". Check founder mapping.")
      }

      tables$edges$add_row(
        left   = df$left[i],
        right  = df$right[i],
        #parent = as.integer(nodeIdMap[[parentKey]]),
        #child  = as.integer(nodeIdMap[[childKey]])
        parent = nodeIdMapByChr[[cc]][[parentKey]],
        child = nodeIdMapByChr[[cc]][[childKey]]
      )
    }

    tables$sort()
    newTs <- tables$tree_sequence()

    outDirCc <- if (is.null(out_dir)) dirname(chr_info[[cc]]$ts_path) else out_dir
    outPath <- file.path(outDirCc, paste0(out_basename, "_chr", cc - 1, ".trees"))

    newTs$dump(outPath)
    cat("Wrote:", outPath, "\n")
  }

  invisible(TRUE)
}
