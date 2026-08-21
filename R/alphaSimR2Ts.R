.bridgeMetadataString <- function(key) {
  key <- as.character(key)
  key <- gsub("\\\\", "\\\\\\\\", key)
  key <- gsub("\"", "\\\\\"", key)
  paste0("{\"alphaSimR\":{\"id\":\"", key, "\"}}")
}

.bridgeIndividualMetadataRaw <- function(file_id) {
  charToRaw(paste0("{\"file_id\":", as.integer(file_id), "}"))
}

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

bridgeCollectSegFromSimOutput <- function(SP, simOutput, pos_list) {
  out <- list()
  nLociByChr <- lapply(pos_list, length)

  if (length(simOutput) >= 2L) {
    for (k in seq.int(2L, length(simOutput))) {
      segDf <- recHistToSegDfWithParents(SP, simOutput[[k]], nLociByChr)
      out[[length(out) + 1L]] <- segDf
    }
  }

  invisible(out)
}

bridgeAllSegToEdgeDf <- function(chr_info, seg_list, pos_list) {
  allSeg <- do.call(rbind, seg_list)

  out <- allSeg
  out$left <- NA
  out$right <- NA

  for (cc in sort(unique(out$chr))) {
    posBp <- pos_list[[cc]]

    tsPath <- chr_info[[cc]]$ts_path
    tc <- RcppTskit::tc_load(tsPath)
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
                             out_basename = "AlphaSimR_extended",
                             ploidy = NULL) {

  indTime <- bridgeComputeIndTime(SP$pedigree)

  if (is.null(ploidy)) {
    ploidy <- suppressWarnings(max(edgeDf$hap, edgeDf$parentHap, na.rm = TRUE))
  }
  ploidy <- as.integer(ploidy)
  if (length(ploidy) != 1L || is.na(ploidy) || ploidy < 1L) {
    stop("ploidy must be supplied or inferable from edgeDf", call. = FALSE)
  }

  outPaths <- character(length(chr_info))

  for (cc in seq_along(chr_info)) {

    nodeIdMap <- new.env(parent = emptyenv())
    indIdMap <- new.env(parent = emptyenv())

    ts <- RcppTskit::ts_load(chr_info[[cc]]$ts_path)
    tc <- ts$dump_tables()

    df <- edgeDf[edgeDf$chr == cc, , drop = FALSE]
    if (nrow(df) == 0) {
      next
    }

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

    if (length(sampNodeId) %% ploidy != 0L) {
      stop("Sample node count is not divisible by ploidy on chr ", cc, call. = FALSE)
    }
    nFounder <- length(sampNodeId) / ploidy
    idx <- 1
    for (ind in 1:nFounder) {
      indRow <- sampIndRow[idx]
      assign(as.character(ind), indRow, envir = indIdMap)

      for (h in 1:ploidy) {
        nodeId <- as.integer(unlist(sampNodeId[[idx]]))[1]
        key <- paste(ind, h, sep = "_")
        assign(key, nodeId, envir = nodeIdMap)
        #  list(alphaSimR = list(id = key)))
        idx <- idx + 1
      }
    }

    # add indIDs for offSpring nodes
    nextInd <- as.integer(tc$num_individuals())
    addNewIndividual <- function(alphaId) {
      key <- as.character(alphaId)
      indRow <- get0(key, envir = indIdMap, inherits = FALSE)
      if (!is.null(indRow)) return(indRow)

      m <- SP$pedigree[alphaId, "mother"]
      f <- SP$pedigree[alphaId, "father"]

      mRow <- addNewIndividual(m)
      fRow <- addNewIndividual(f)

      newId <- nextInd
      tc$individual_table_add_row(
        #parents = list(as.integer(mRow), as.integer(fRow)),
        parents = c(as.integer(mRow), as.integer(fRow)),
        metadata = .bridgeIndividualMetadataRaw(newId))

      assign(key, as.integer(newId), envir = indIdMap)

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
      if (is.null(get0(key, envir = nodeIdMap, inherits = FALSE))) {
        childId <- as.integer(sub("_.*$", "", key))
        indRow  <- get(as.character(childId), envir = indIdMap, inherits = FALSE)

        tc$node_table_add_row(
          flags = 0L,
          time  = indTime[[childId]],
          population = -1L,
          individual = indRow,
          metadata = .bridgeMetadataString(key)
        )
        assign(key, as.integer(tc$num_nodes() - 1), envir = nodeIdMap)
      }
    }

    # append edges
    for (i in 1:nrow(df)) {
      parentKey <- paste(df$parentId[i], df$parentHap[i], sep = "_")
      childKey  <- paste(df$childId[i], df$hap[i], sep = "_")

      parentNode <- get0(parentKey, envir = nodeIdMap, inherits = FALSE)
      if (is.null(parentNode)) {
        stop("Missing parent node for key=", parentKey,
             " on chr=", cc, ". Check founder mapping.")
      }
      childNode <- get(childKey, envir = nodeIdMap, inherits = FALSE)

      tc$edge_table_add_row(
        left   = df$left[i],
        right  = df$right[i],
        parent = parentNode,
        child = childNode
      )
    }

    tc$sort()
    newTs <- tc$tree_sequence()

    outDirCc <- if (is.null(out_dir)) dirname(chr_info[[cc]]$ts_path) else out_dir
    outPath <- file.path(outDirCc, paste0(out_basename, "_chr", cc - 1, ".trees"))

    newTs$dump(outPath)
    outPaths[[cc]] <- outPath
    cat("Wrote:", outPath, "\n")
  }

  invisible(outPaths[nzchar(outPaths)])
}
