library(jsonlite)

recHistMatToSegDf <- function(histMat, nLoci) {

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


recHistToSegDfWithParents <- function(SP, offspringPop, nLociByChr) {
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
        seg <- recHistMatToSegDf(haps[[h]], nLoci = nLoci)

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

bridgeCollectSegFromSimOutput <- function(SP, simOutput) {
  bridgeSegDfList <<- list()

  nLociByChr <- lapply(chrKeptPosBpList, length)

  for (k in 2:length(simOutput)) {
    segDf <- recHistToSegDfWithParents(SP, simOutput[[k]], nLociByChr)
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

bridgeAllSegToEdgeDf <- function(chr_info) {
  allSeg <- do.call(rbind, bridgeSegDfList)

  out <- allSeg
  out$left <- NA
  out$right <- NA

  for (cc in sort(unique(out$chr))) {
    posBp <- chrKeptPosBpList[[cc]]

    tsPath <- chr_info[[cc]]$ts_path
    tc <- tc_load(tsPath)
    seqLen <- tc$sequence_length()

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

bridgeComputeIndTime <- function(pedigree) {
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


bridgeWriteTrees <- function(chr_info, edgeDf, SP, out_dir = NULL,
                             out_basename = "AlphaSimR_extended") {

  indTime <- bridgeComputeIndTime(SP$pedigree)

  nodeIdMapByChr <<- vector("list", length(chr_info))
  indIdMapByChr  <<- vector("list", length(chr_info))

  for (cc in seq_along(chr_info)) {

    nodeIdMapByChr[[cc]] <<- list()
    indIdMapByChr[[cc]]  <<- list()

    ts <- ts_load(chr_info[[cc]]$ts_path)
    tc <- ts$dump_tables()

    df <- edgeDf[edgeDf$chr == cc, , drop = FALSE]
    if (nrow(df) == 0) next

    # get indIDs for sampled nodes
    sampNodeId <- ts$samples()
    sampIndRow <- integer(length(sampNodeId))
    for (i in seq_along(sampNodeId)) {
      sampIndRow[i] <- tc$node_table_get_row(sampNodeId[i])$individual
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
        nodeId <- as.integer(unlist(sampNodeId[[idx]]))[1]
        key <- paste(ind, h, sep = "_")
        nodeIdMapByChr[[cc]][[key]] <<- nodeId
        #  list(alphaSimR = list(id = key)))
        idx <- idx + 1
      }
    }

    # add indIDs for offSpring nodes
    nextInd <- as.integer(tc$num_individuals())
    addNewIndividual <- function(alphaId) {
      key <- as.character(alphaId)
      if (!is.null(indIdMapByChr[[cc]][[key]])) return(indIdMapByChr[[cc]][[key]])

      m <- SP$pedigree[alphaId, "mother"]
      f <- SP$pedigree[alphaId, "father"]

      mRow <- addNewIndividual(m)
      fRow <- addNewIndividual(f)

      newId <- nextInd
      tc$individual_table_add_row(
        #parents = list(as.integer(mRow), as.integer(fRow)),
        parents = c(as.integer(mRow), as.integer(fRow)),
        metadata = charToRaw(toJSON(
          list(file_id=as.integer(newId)),
        auto_unbox = TRUE)))

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
      if (is.null(nodeIdMapByChr[[cc]][[key]])) {
        childId <- as.integer(sub("_.*$", "", key))
        indRow  <- indIdMapByChr[[cc]][[as.character(childId)]]

        tc$node_table_add_row(
          flags = 0L,
          time  = indTime[[childId]],
          population = -1L,
          individual = indRow,
          metadata = as.character(toJSON(
            list(alphaSimR = list(id = key)),
            auto_unbox = TRUE, force = TRUE))
        )
        nodeIdMapByChr[[cc]][[key]] <<- as.integer(tc$num_nodes() - 1)
      }
    }

    # append edges
    for (i in 1:nrow(df)) {
      parentKey <- paste(df$parentId[i], df$parentHap[i], sep = "_")
      childKey  <- paste(df$childId[i], df$hap[i], sep = "_")

      if (is.null(nodeIdMapByChr[[cc]][[parentKey]])) {
        stop("Missing parent node for key=", parentKey,
             " on chr=", cc, ". Check founder mapping.")
      }

      tc$edge_table_add_row(
        left   = df$left[i],
        right  = df$right[i],
        parent = nodeIdMapByChr[[cc]][[parentKey]],
        child = nodeIdMapByChr[[cc]][[childKey]]
      )
    }

    tc$sort()
    newTs <- tc$tree_sequence()

    outDirCc <- if (is.null(out_dir)) dirname(chr_info[[cc]]$ts_path) else out_dir
    outPath <- file.path(outDirCc, paste0(out_basename, "_chr", cc - 1, ".trees"))

    newTs$dump(outPath)
    cat("Wrote:", outPath, "\n")
  }

  invisible(TRUE)
}
