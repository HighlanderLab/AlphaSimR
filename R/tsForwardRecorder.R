.tsForwardNodeKey <- function(iid, hap) {
  paste0(as.integer(iid), "_", as.integer(hap))
}

.tsForwardAttachMetadata <- function() {
  isTRUE(getOption("AlphaSimR.tsForwardAttachMetadata", TRUE))
}

.tsForwardIndividualMetadataRaw <- function(file_id) {
  charToRaw(paste0("{\"file_id\":", as.integer(file_id), "}"))
}

.tsForwardUpdateIndTimeCache <- function(indTime, pedigree) {
  n <- nrow(pedigree)
  if (is.null(indTime)) {
    indTime <- numeric(0)
  }
  oldN <- length(indTime)
  if (oldN >= n) {
    return(indTime)
  }
  indTime <- c(indTime, rep(NA_real_, n - oldN))
  for (i in seq.int(oldN + 1L, n)) {
    m <- pedigree[i, "mother"]
    f <- pedigree[i, "father"]
    if (m == 0L && f == 0L) {
      indTime[i] <- 0
    } else {
      indTime[i] <- min(indTime[m], indTime[f]) - 1
    }
  }
  indTime
}

.tsForwardEnsureIndividual <- function(chrState, iid, pedigree) {
  key <- as.character(as.integer(iid))
  if (exists(key, envir = chrState$indMap, inherits = FALSE)) {
    return(get(key, envir = chrState$indMap, inherits = FALSE))
  }

  m <- pedigree[iid, "mother"]
  f <- pedigree[iid, "father"]
  nextInd <- as.integer(chrState$tc$num_individuals())
  md <- if (.tsForwardAttachMetadata()) .tsForwardIndividualMetadataRaw(nextInd) else NULL
  if (m > 0L) {
    mRow <- .tsForwardEnsureIndividual(chrState, m, pedigree)
    fRow <- .tsForwardEnsureIndividual(chrState, f, pedigree)
    indRow <- chrState$tc$individual_table_add_row(
      parents = c(as.integer(mRow), as.integer(fRow)),
      metadata = md
    )
  } else {
    indRow <- chrState$tc$individual_table_add_row(
      parents = NULL,
      metadata = md
    )
  }

  assign(key, as.integer(indRow), envir = chrState$indMap)
  as.integer(indRow)
}

.tsForwardMorgToTsPosVec <- function(m, gm, pos, seqLen, side = c("left", "right")) {
  side <- match.arg(side)
  out <- rep(NA_real_, length(m))
  ok <- is.finite(m)
  if (!any(ok)) {
    return(out)
  }
  if (length(gm) < 2L) {
    v <- min(max(pos[1], 0), seqLen)
    out[ok] <- v
    return(out)
  }

  out[ok & m <= 0] <- 0
  out[ok & m >= gm[length(gm)]] <- seqLen

  mid <- ok & m > 0 & m < gm[length(gm)]
  if (any(mid)) {
    mm <- m[mid]
    i <- findInterval(mm, gm, rightmost.closed = TRUE)
    i <- pmin.int(pmax.int(i, 1L), length(gm) - 1L)
    gmL <- gm[i]
    gmR <- gm[i + 1L]
    posL <- pos[i]
    posR <- pos[i + 1L]
    val <- posL + (mm - gmL) * (posR - posL) / (gmR - gmL)
    bad <- gmR <= gmL
    if (any(bad)) {
      val[bad] <- if (side == "left") posL[bad] else posR[bad]
    }
    out[mid] <- pmin(seqLen, pmax(0, val))
  }
  out
}

.tsForwardMorgToTsPosRateMapVec <- function(m, x0, breaks, rates, seqLen, maxM,
                                            side = c("left", "right")) {
  side <- match.arg(side)
  out <- rep(NA_real_, length(m))
  ok <- is.finite(m)
  if (!any(ok)) {
    return(out)
  }

  breaks <- as.numeric(breaks)
  rates <- as.numeric(rates)
  if (length(breaks) != length(rates) + 1L ||
      length(rates) < 1L ||
      any(!is.finite(breaks)) ||
      any(!is.finite(rates)) ||
      any(rates < 0) ||
      any(diff(breaks) < 0)) {
    stop("Invalid TS forward rate-map metadata", call. = FALSE)
  }

  if (!is.finite(maxM) || maxM <= 0) {
    out[ok] <- if (side == "left") 0 else seqLen
    return(out)
  }

  out[ok & m <= 0] <- 0
  out[ok & m >= maxM] <- seqLen

  mid <- ok & m > 0 & m < maxM
  if (any(mid)) {
    segLen <- diff(breaks)
    mStart <- c(0, cumsum(rates * segLen))
    mEnd <- mStart[-1L]

    i0 <- findInterval(x0, breaks, rightmost.closed = TRUE)
    i0 <- pmin.int(pmax.int(i0, 1L), length(rates))
    mX0 <- mStart[i0] + rates[i0] * (x0 - breaks[i0])

    M <- m[mid] + mX0
    i <- findInterval(M, mStart, rightmost.closed = TRUE)
    i <- pmin.int(pmax.int(i, 1L), length(rates))

    plateau <- (rates[i] == 0) | (mEnd[i] == mStart[i])
    val <- numeric(length(M))
    nonPlateau <- !plateau
    if (any(nonPlateau)) {
      val[nonPlateau] <- breaks[i[nonPlateau]] +
        (M[nonPlateau] - mStart[i[nonPlateau]]) / rates[i[nonPlateau]]
    }
    if (any(plateau)) {
      val[plateau] <- if (side == "left") breaks[i[plateau]] else breaks[i[plateau] + 1L]
    }
    out[mid] <- pmin(seqLen, pmax(0, val))
  }
  out
}

.tsForwardHasRateMapMeta <- function(posMeta, cc) {
  !is.null(posMeta$breaksList) &&
    !is.null(posMeta$ratesList) &&
    length(posMeta$breaksList) >= cc &&
    length(posMeta$ratesList) >= cc &&
    !is.null(posMeta$breaksList[[cc]]) &&
    !is.null(posMeta$ratesList[[cc]])
}

.tsForwardSegGenToTsIntervals <- function(chr, parentSide, leftGen, rightGen,
                                          femaleMap, maleMap, posMeta) {
  nSeg <- length(chr)
  left <- as.numeric(leftGen)
  right <- as.numeric(rightGen)
  if (nSeg == 0L || is.null(posMeta) || is.null(posMeta$posList)) {
    return(list(left = left, right = right))
  }

  outLeft <- numeric(nSeg)
  outRight <- numeric(nSeg)
  chrVals <- sort.int(unique(as.integer(chr)), method = "quick")
  for (cc in chrVals) {
    idxChr <- which(chr == cc)
    pos <- posMeta$posList[[cc]]
    seqLen <- as.numeric(posMeta$seqLenList[[cc]])
    useRateMap <- .tsForwardHasRateMapMeta(posMeta, cc)
    x0 <- if (useRateMap) as.numeric(pos[[1L]]) else NA_real_
    breaks <- if (useRateMap) posMeta$breaksList[[cc]] else NULL
    rates <- if (useRateMap) posMeta$ratesList[[cc]] else NULL

    idxM <- idxChr[parentSide[idxChr] == 1L]
    if (length(idxM) > 0L) {
      gm <- femaleMap[[cc]]
      if (useRateMap) {
        outLeft[idxM] <- .tsForwardMorgToTsPosRateMapVec(
          left[idxM], x0, breaks, rates, seqLen, max(gm), side = "left"
        )
        outRight[idxM] <- .tsForwardMorgToTsPosRateMapVec(
          right[idxM], x0, breaks, rates, seqLen, max(gm), side = "right"
        )
      } else {
        outLeft[idxM] <- .tsForwardMorgToTsPosVec(left[idxM], gm, pos, seqLen, side = "left")
        outRight[idxM] <- .tsForwardMorgToTsPosVec(right[idxM], gm, pos, seqLen, side = "right")
      }
    }

    idxF <- idxChr[parentSide[idxChr] == 2L]
    if (length(idxF) > 0L) {
      gm <- maleMap[[cc]]
      if (useRateMap) {
        outLeft[idxF] <- .tsForwardMorgToTsPosRateMapVec(
          left[idxF], x0, breaks, rates, seqLen, max(gm), side = "left"
        )
        outRight[idxF] <- .tsForwardMorgToTsPosRateMapVec(
          right[idxF], x0, breaks, rates, seqLen, max(gm), side = "right"
        )
      } else {
        outLeft[idxF] <- .tsForwardMorgToTsPosVec(left[idxF], gm, pos, seqLen, side = "left")
        outRight[idxF] <- .tsForwardMorgToTsPosVec(right[idxF], gm, pos, seqLen, side = "right")
      }
    }
  }

  list(left = outLeft, right = outRight)
}

.tsForwardResolvePosMeta <- function() {
  if (exists("chrKeptPosTsList", inherits = TRUE)) {
    posList <- get("chrKeptPosTsList", inherits = TRUE)
  } else if (exists("chrKeptPosBpList", inherits = TRUE)) {
    posList <- get("chrKeptPosBpList", inherits = TRUE)
  } else {
    return(list(posList = NULL, seqLenList = NULL))
  }
  if (exists("chrSeqLenTsList", inherits = TRUE)) {
    seqLenList <- get("chrSeqLenTsList", inherits = TRUE)
  } else if (exists("chrSeqLenBpList", inherits = TRUE)) {
    seqLenList <- get("chrSeqLenBpList", inherits = TRUE)
  } else {
    seqLenList <- lapply(posList, function(x) max(x, na.rm = TRUE))
  }
  list(posList = posList, seqLenList = seqLenList)
}

.tsForwardParseNodeKeys <- function(keys) {
  list(
    iids = as.integer(sub("_.*$", "", keys)),
    haps = as.integer(sub("^.*_", "", keys))
  )
}

.tsForwardAddMissingChildNodes <- function(chrState, missingKeys, nodeTimes) {
  if (length(missingKeys) == 0L) {
    return(invisible(NULL))
  }
  if (length(nodeTimes) != length(missingKeys)) {
    stop("Internal error: nodeTimes length does not match missingKeys length", call. = FALSE)
  }

  keyInfo <- .tsForwardParseNodeKeys(missingKeys)
  indRows <- vapply(keyInfo$iids, function(iid) {
    as.integer(get(as.character(iid), envir = chrState$indMap, inherits = FALSE))
  }, integer(1))

  if (!.tsForwardAttachMetadata()) {
    nodeIds <- tsForwardNodeTableAddRows(
      chrState$tc,
      flags = rep.int(0L, length(missingKeys)),
      time = as.numeric(nodeTimes),
      population = rep.int(-1L, length(missingKeys)),
      individual = as.integer(indRows)
    )
  } else {
    nodeIds <- tsForwardNodeTableAddRowsWithMetadata(
      chrState$tc,
      flags = rep.int(0L, length(missingKeys)),
      time = as.numeric(nodeTimes),
      population = rep.int(-1L, length(missingKeys)),
      individual = as.integer(indRows),
      nodeKey = as.character(missingKeys)
    )
  }

  for (j in seq_along(missingKeys)) {
    assign(missingKeys[[j]], as.integer(nodeIds[j]), envir = chrState$nodeMap)
  }
  invisible(NULL)
}

.tsForwardAppendEdgesByKeys <- function(chrState, left, right, parentKeys, childKeys, cc) {
  n <- length(left)
  if (length(right) != n || length(parentKeys) != n || length(childKeys) != n) {
    stop("Internal error: edge vector lengths do not match", call. = FALSE)
  }

  keep <- rep.int(FALSE, n)
  leftVec <- numeric(n)
  rightVec <- numeric(n)
  parentVec <- integer(n)
  childVec <- integer(n)
  for (i in seq_len(n)) {
    leftI <- max(0, as.numeric(left[[i]]))
    rightI <- min(chrState$seqLen, as.numeric(right[[i]]))
    if (!is.finite(leftI) || !is.finite(rightI) || rightI <= leftI) {
      next
    }

    pKey <- parentKeys[[i]]
    cKey <- childKeys[[i]]
    if (!exists(pKey, envir = chrState$nodeMap, inherits = FALSE)) {
      stop("Missing parent node for key ", pKey, " on chr ", cc, call. = FALSE)
    }
    if (!exists(cKey, envir = chrState$nodeMap, inherits = FALSE)) {
      stop("Missing child node for key ", cKey, " on chr ", cc, call. = FALSE)
    }

    keep[[i]] <- TRUE
    leftVec[[i]] <- leftI
    rightVec[[i]] <- rightI
    parentVec[[i]] <- as.integer(get(pKey, envir = chrState$nodeMap, inherits = FALSE))
    childVec[[i]] <- as.integer(get(cKey, envir = chrState$nodeMap, inherits = FALSE))
  }

  if (any(keep)) {
    idx <- which(keep)
    tsForwardEdgeTableAddRows(
      chrState$tc,
      left = leftVec[idx],
      right = rightVec[idx],
      parent = parentVec[idx],
      child = childVec[idx]
    )
  }
  invisible(NULL)
}

.tsForwardEnsureChildrenIndividuals <- function(chrState,
                                                children,
                                                simParam,
                                                cc,
                                                childIid = NULL,
                                                childMotherIid = NULL,
                                                childFatherIid = NULL) {
  useDirectParents <- !is.null(childMotherIid) && !is.null(childFatherIid)

  for (iid in children) {
    key <- as.character(as.integer(iid))
    if (exists(key, envir = chrState$indMap, inherits = FALSE)) {
      next
    }

    if (useDirectParents) {
      ii <- match(iid, childIid)
      if (is.na(ii)) {
        stop("Internal error: child iid not found in childIid map", call. = FALSE)
      }
      mIid <- childMotherIid[[ii]]
      fIid <- childFatherIid[[ii]]
      mKey <- as.character(as.integer(mIid))
      fKey <- as.character(as.integer(fIid))

      if (!exists(mKey, envir = chrState$indMap, inherits = FALSE)) {
        if (mIid <= nrow(simParam$pedigree)) {
          .tsForwardEnsureIndividual(chrState, mIid, simParam$pedigree)
        } else {
          stop("Missing mother individual row for iid ", mIid, " on chr ", cc, call. = FALSE)
        }
      }
      if (!exists(fKey, envir = chrState$indMap, inherits = FALSE)) {
        if (fIid <= nrow(simParam$pedigree)) {
          .tsForwardEnsureIndividual(chrState, fIid, simParam$pedigree)
        } else {
          stop("Missing father individual row for iid ", fIid, " on chr ", cc, call. = FALSE)
        }
      }

      mRow <- as.integer(get(mKey, envir = chrState$indMap, inherits = FALSE))
      fRow <- as.integer(get(fKey, envir = chrState$indMap, inherits = FALSE))
      nextInd <- as.integer(chrState$tc$num_individuals())
      md <- if (.tsForwardAttachMetadata()) .tsForwardIndividualMetadataRaw(nextInd) else NULL
      indRow <- chrState$tc$individual_table_add_row(
        parents = c(mRow, fRow),
        metadata = md
      )
      assign(key, as.integer(indRow), envir = chrState$indMap)
    } else {
      .tsForwardEnsureIndividual(chrState, iid, simParam$pedigree)
    }
  }

  invisible(NULL)
}

.tsSegGenToTsPos <- function(tsSegGen, femaleMap, maleMap) {
  if (is.null(tsSegGen)) {
    return(NULL)
  }
  if (is.null(dim(tsSegGen)) || nrow(tsSegGen) == 0L) {
    return(data.frame(
      childLocal = integer(),
      chr = integer(),
      hap = integer(),
      parentSide = integer(),
      parentIndex = integer(),
      parentHap = integer(),
      leftGen = numeric(),
      rightGen = numeric(),
      left = numeric(),
      right = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  df <- as.data.frame(tsSegGen, stringsAsFactors = FALSE)
  names(df) <- c("childLocal", "chr", "hap", "parentSide",
                 "parentIndex", "parentHap", "leftGen", "rightGen")
  intCols <- c("childLocal", "chr", "hap", "parentSide", "parentIndex", "parentHap")
  for (cc in intCols) {
    df[[cc]] <- as.integer(df[[cc]])
  }

  posMeta <- .tsForwardResolvePosMeta()
  intervals <- .tsForwardSegGenToTsIntervals(
    chr = df$chr,
    parentSide = df$parentSide,
    leftGen = df$leftGen,
    rightGen = df$rightGen,
    femaleMap = femaleMap,
    maleMap = maleMap,
    posMeta = posMeta
  )
  df$left <- intervals$left
  df$right <- intervals$right

  df[df$right > df$left, , drop = FALSE]
}

.attachTsSegGenIds <- function(tsSegGen, childPop, femalePop, malePop) {
  if (is.null(tsSegGen)) {
    return(NULL)
  }
  if (nrow(tsSegGen) == 0L) {
    tsSegGen$childIid <- integer()
    tsSegGen$parentIid <- integer()
    tsSegGen$childId <- character()
    tsSegGen$parentId <- character()
    return(tsSegGen)
  }

  nChild <- nInd(childPop)
  if (any(tsSegGen$childLocal < 1L | tsSegGen$childLocal > nChild)) {
    stop("Internal error: tsSegGen childLocal out of bounds", call. = FALSE)
  }

  motherRows <- which(tsSegGen$parentSide == 1L)
  fatherRows <- which(tsSegGen$parentSide == 2L)
  if (length(motherRows) > 0L) {
    if (any(tsSegGen$parentIndex[motherRows] < 1L |
            tsSegGen$parentIndex[motherRows] > nInd(femalePop))) {
      stop("Internal error: tsSegGen mother parentIndex out of bounds", call. = FALSE)
    }
  }
  if (length(fatherRows) > 0L) {
    if (any(tsSegGen$parentIndex[fatherRows] < 1L |
            tsSegGen$parentIndex[fatherRows] > nInd(malePop))) {
      stop("Internal error: tsSegGen father parentIndex out of bounds", call. = FALSE)
    }
  }

  tsSegGen$childIid <- childPop@iid[tsSegGen$childLocal]
  tsSegGen$childId <- childPop@id[tsSegGen$childLocal]

  tsSegGen$parentIid <- integer(nrow(tsSegGen))
  tsSegGen$parentId <- character(nrow(tsSegGen))
  if (length(motherRows) > 0L) {
    tsSegGen$parentIid[motherRows] <- femalePop@iid[tsSegGen$parentIndex[motherRows]]
    tsSegGen$parentId[motherRows] <- femalePop@id[tsSegGen$parentIndex[motherRows]]
  }
  if (length(fatherRows) > 0L) {
    tsSegGen$parentIid[fatherRows] <- malePop@iid[tsSegGen$parentIndex[fatherRows]]
    tsSegGen$parentId[fatherRows] <- malePop@id[tsSegGen$parentIndex[fatherRows]]
  }
  tsSegGen
}

.appendTsForwardRawIfActive <- function(simParam,
                                        tsSegGenRaw,
                                        childIid,
                                        motherIid,
                                        fatherIid,
                                        childMotherIid = NULL,
                                        childFatherIid = NULL,
                                        femaleMap,
                                        maleMap) {
  if (!simParam$isTrackRecGen) {
    return(FALSE)
  }
  rec <- .tsForwardGetRecorder(simParam)
  if (is.null(rec)) {
    return(FALSE)
  }
  rec <- tsForwardAppendSegGenRaw(
    recorder = rec,
    tsSegGen = tsSegGenRaw,
    childIid = childIid,
    motherIid = motherIid,
    fatherIid = fatherIid,
    childMotherIid = childMotherIid,
    childFatherIid = childFatherIid,
    femaleMap = femaleMap,
    maleMap = maleMap,
    simParam = simParam
  )
  .tsForwardSetRecorder(simParam, rec)
  TRUE
}

.useTsForwardDirectCross <- function(simParam) {
  if (!simParam$isTrackRecGen) {
    return(FALSE)
  }
  !is.null(.tsForwardGetRecorder(simParam))
}

.keepTsForwardSeg <- function() {
  isTRUE(getOption("AlphaSimR.tsForwardKeepSeg", FALSE))
}

.needRecHistGenFromCross <- function(simParam) {
  if (!isTRUE(simParam$isTrackRecGen)) {
    return(FALSE)
  }
  isTRUE(getOption("AlphaSimR.tsForwardKeepRecHistGen", FALSE))
}

.setupTsForwardCross <- function(simParam, crossPlan, motherIid, fatherIid, femaleMap, maleMap) {
  keepRecHistGen <- .needRecHistGenFromCross(simParam)
  useTsDirect <- .useTsForwardDirectCross(simParam)
  returnTsSegGen <- .keepTsForwardSeg()

  childIidPred <- NULL
  directAppendFn <- NULL
  if (useTsDirect) {
    childIidPred <- simParam$lastId + seq_len(nrow(crossPlan))
    childMotherIidPred <- motherIid[crossPlan[, 1]]
    childFatherIidPred <- fatherIid[crossPlan[, 2]]
    directAppendFn <- function(segChunk) {
      .appendTsForwardRawIfActive(
        simParam = simParam,
        tsSegGenRaw = segChunk,
        childIid = childIidPred,
        motherIid = motherIid,
        fatherIid = fatherIid,
        childMotherIid = childMotherIidPred,
        childFatherIid = childFatherIidPred,
        femaleMap = femaleMap,
        maleMap = maleMap
      )
    }
  }

  list(
    keepRecHistGen = keepRecHistGen,
    useTsDirect = useTsDirect,
    returnTsSegGen = returnTsSegGen,
    childIidPred = childIidPred,
    directAppendFn = directAppendFn
  )
}

.extractTsCrossRecGen <- function(tmp, simParam, keepRecHistGen) {
  if (simParam$isTrackRecGen) {
    return(list(
      histGen = if (keepRecHistGen) tmp$recHistGen else NULL,
      tsSegGenRaw = tmp$tsSegGen
    ))
  }
  list(histGen = NULL, tsSegGenRaw = NULL)
}

.finalizeTsForwardCross <- function(outPop,
                                    tsSegGenRaw,
                                    recorderAppended,
                                    childIidPred,
                                    femalePop,
                                    malePop,
                                    femaleMap,
                                    maleMap) {
  if (recorderAppended && !is.null(childIidPred) &&
      !identical(as.integer(outPop@iid), as.integer(childIidPred))) {
    stop("Internal error: predicted child iid does not match realized iid in direct TS mode", call. = FALSE)
  }

  if (is.null(tsSegGenRaw)) {
    return(outPop)
  }

  if (.keepTsForwardSeg()) {
    tsSegGen <- .tsSegGenToTsPos(tsSegGenRaw, femaleMap, maleMap)
    outPop@misc$tsSegGen <- .attachTsSegGenIds(tsSegGen, outPop, femalePop, malePop)
  }
  outPop
}

#' Append raw meiosis tsSegGen matrix to forward TS recorder
#'
#' @param recorder output of `tsForwardInit()`
#' @param tsSegGen numeric matrix with columns
#'   `childLocal, chr, hap, parentSide, parentIndex, parentHap, leftGen, rightGen`
#' @param childIid integer vector mapping child local index to iid
#' @param motherIid integer vector mapping mother parentIndex to iid
#' @param fatherIid integer vector mapping father parentIndex to iid
#' @param femaleMap female genetic maps
#' @param maleMap male genetic maps
#' @param simParam `SimParam` with pedigree
#' @return Updated recorder state
#' @keywords internal
#' @noRd
tsForwardAppendSegGenRaw <- function(recorder,
                                     tsSegGen,
                                     childIid,
                                     motherIid,
                                     fatherIid,
                                     childMotherIid = NULL,
                                     childFatherIid = NULL,
                                     femaleMap,
                                     maleMap,
                                     simParam) {
  stopifnot(inherits(recorder, "tsForwardRecorder"))
  if (is.null(tsSegGen) || length(tsSegGen) == 0L || nrow(tsSegGen) == 0L) {
    return(recorder)
  }
  seg <- as.matrix(tsSegGen)
  if (ncol(seg) != 8L) {
    stop("tsSegGen must have 8 columns", call. = FALSE)
  }

  childIid <- as.integer(childIid)
  motherIid <- as.integer(motherIid)
  fatherIid <- as.integer(fatherIid)
  if (!is.null(childMotherIid)) {
    childMotherIid <- as.integer(childMotherIid)
  }
  if (!is.null(childFatherIid)) {
    childFatherIid <- as.integer(childFatherIid)
  }

  childLocal <- as.integer(seg[, 1])
  chr <- as.integer(seg[, 2])
  hap <- as.integer(seg[, 3])
  parentSide <- as.integer(seg[, 4])
  parentIndex <- as.integer(seg[, 5])
  parentHap <- as.integer(seg[, 6])
  leftGen <- as.numeric(seg[, 7])
  rightGen <- as.numeric(seg[, 8])
  nSeg <- length(childLocal)

  if (any(chr < 1L | chr > recorder$nChr)) {
    stop("tsSegGen has chromosome index out of bounds", call. = FALSE)
  }
  if (any(hap < 1L | hap > recorder$ploidy)) {
    stop("tsSegGen has child haplotype out of bounds", call. = FALSE)
  }
  if (any(parentHap < 1L)) {
    stop("tsSegGen has invalid parent haplotype", call. = FALSE)
  }
  if (any(childLocal < 1L | childLocal > length(childIid))) {
    stop("tsSegGen childLocal out of bounds", call. = FALSE)
  }
  if (any(!(parentSide %in% c(1L, 2L)))) {
    stop("tsSegGen parentSide must be 1 or 2", call. = FALSE)
  }
  motherRows <- which(parentSide == 1L)
  fatherRows <- which(parentSide == 2L)
  if (any(parentIndex[motherRows] < 1L | parentIndex[motherRows] > length(motherIid))) {
    stop("tsSegGen mother parentIndex out of bounds", call. = FALSE)
  }
  if (any(parentIndex[fatherRows] < 1L | parentIndex[fatherRows] > length(fatherIid))) {
    stop("tsSegGen father parentIndex out of bounds", call. = FALSE)
  }

  childIidBySeg <- childIid[childLocal]
  if (!is.null(childMotherIid) || !is.null(childFatherIid)) {
    if (is.null(childMotherIid) || is.null(childFatherIid)) {
      stop("Both childMotherIid and childFatherIid must be supplied together", call. = FALSE)
    }
    if (length(childMotherIid) != length(childIid) ||
        length(childFatherIid) != length(childIid)) {
      stop("childMotherIid/childFatherIid must match length(childIid)", call. = FALSE)
    }
    if (anyDuplicated(childIid) > 0L) {
      stop("childIid must be unique when child parent iid vectors are provided", call. = FALSE)
    }
  }
  parentIid <- integer(nSeg)
  parentIid[motherRows] <- motherIid[parentIndex[motherRows]]
  parentIid[fatherRows] <- fatherIid[parentIndex[fatherRows]]

  if (is.null(recorder$posMeta)) {
    recorder$posMeta <- .tsForwardResolvePosMeta()
  }
  posMeta <- recorder$posMeta
  intervals <- .tsForwardSegGenToTsIntervals(
    chr = chr,
    parentSide = parentSide,
    leftGen = leftGen,
    rightGen = rightGen,
    femaleMap = femaleMap,
    maleMap = maleMap,
    posMeta = posMeta
  )
  left <- intervals$left
  right <- intervals$right

  recorder$indTime <- .tsForwardUpdateIndTimeCache(recorder$indTime, simParam$pedigree)
  indTime <- recorder$indTime
  for (cc in seq_len(recorder$nChr)) {
    chrState <- recorder$chr[[cc]]
    idxChr <- which(chr == cc)
    if (length(idxChr) == 0L) {
      next
    }

    # For direct append mode, child pedigree rows may not exist yet.
    # Derive a fallback child-node time from parent node times in this chunk.
    parentMinTimeByChildNodeKey <- new.env(parent = emptyenv(), hash = TRUE)
    pKeysChr <- .tsForwardNodeKey(parentIid[idxChr], parentHap[idxChr])
    cKeysChr <- .tsForwardNodeKey(childIidBySeg[idxChr], hap[idxChr])
    for (k in seq_along(idxChr)) {
      pKey <- pKeysChr[[k]]
      cKey <- cKeysChr[[k]]
      if (!exists(pKey, envir = chrState$nodeMap, inherits = FALSE)) {
        next
      }
      pNode <- as.integer(get(pKey, envir = chrState$nodeMap, inherits = FALSE))
      pTime <- as.numeric(chrState$tc$node_table_get_row(pNode)$time)
      if (!is.finite(pTime)) {
        next
      }
      if (!exists(cKey, envir = parentMinTimeByChildNodeKey, inherits = FALSE)) {
        assign(cKey, pTime, envir = parentMinTimeByChildNodeKey)
      } else {
        old <- as.numeric(get(cKey, envir = parentMinTimeByChildNodeKey, inherits = FALSE))
        if (pTime < old) {
          assign(cKey, pTime, envir = parentMinTimeByChildNodeKey)
        }
      }
    }

    children <- unique(childIidBySeg[idxChr])
    .tsForwardEnsureChildrenIndividuals(
      chrState = chrState,
      children = children,
      simParam = simParam,
      cc = cc,
      childIid = childIid,
      childMotherIid = childMotherIid,
      childFatherIid = childFatherIid
    )

    childNodeKeys <- unique(.tsForwardNodeKey(childIidBySeg[idxChr], hap[idxChr]))
    missingKeys <- childNodeKeys[!vapply(childNodeKeys, exists, logical(1),
                                         envir = chrState$nodeMap, inherits = FALSE)]
    if (length(missingKeys) > 0L) {
      keyInfo <- .tsForwardParseNodeKeys(missingKeys)
      nodeTimes <- numeric(length(missingKeys))
      for (j in seq_along(missingKeys)) {
        iid <- keyInfo$iids[[j]]
        t <- NA_real_
        if (iid <= length(indTime)) {
          t <- indTime[[iid]]
        }
        if (!is.finite(t) &&
            exists(missingKeys[[j]], envir = parentMinTimeByChildNodeKey, inherits = FALSE)) {
          t <- as.numeric(get(missingKeys[[j]],
                              envir = parentMinTimeByChildNodeKey,
                              inherits = FALSE)) - 1
        }
        if (!is.finite(t) && !is.null(childMotherIid) && !is.null(childFatherIid)) {
          ii <- match(iid, childIid)
          if (!is.na(ii)) {
            mIid <- childMotherIid[[ii]]
            fIid <- childFatherIid[[ii]]
            if (mIid <= length(indTime) && fIid <= length(indTime)) {
              t <- min(indTime[[mIid]], indTime[[fIid]]) - 1
            }
          }
        }
        if (!is.finite(t)) {
          stop("Could not resolve node time for child iid ", iid, " on chr ", cc, call. = FALSE)
        }
        nodeTimes[[j]] <- t
      }
      .tsForwardAddMissingChildNodes(chrState, missingKeys, nodeTimes)
    }

    .tsForwardAppendEdgesByKeys(
      chrState = chrState,
      left = left[idxChr],
      right = right[idxChr],
      parentKeys = pKeysChr,
      childKeys = cKeysChr,
      cc = cc
    )
  }

  recorder
}

.tsForwardHasSourceFields <- function(x) {
  is.list(x) && (
    !is.null(x[["ts_path"]]) ||
      !is.null(x[["ts"]]) ||
      !is.null(x[["ts_xptr"]]) ||
      !is.null(x[["tc_xptr"]]) ||
      !is.null(x[["table_xptr"]])
  )
}

.tsForwardSourceByChr <- function(founderPop) {
  nChr <- founderPop@nChr
  src <- attr(founderPop, "tsForwardSource", exact = TRUE)
  if (is.null(src)) {
    stop(
      "No TS source found on founderPop. Build founders with asMapPop() before setTrackTs().",
      call. = FALSE
    )
  }
  if (!is.list(src) || length(src) != nChr) {
    stop("Invalid TS source metadata on founderPop", call. = FALSE)
  }
  for (cc in seq_len(nChr)) {
    if (!.tsForwardHasSourceFields(src[[cc]])) {
      stop("Missing TS source for chromosome ", cc, call. = FALSE)
    }
  }
  src
}

.tsForwardPosMetaFromFounder <- function(founderPop) {
  meta <- attr(founderPop, "tsForwardPosMeta", exact = TRUE)
  if (is.null(meta)) {
    meta <- .tsForwardResolvePosMeta()
  }
  if (is.null(meta) || is.null(meta$posList) || is.null(meta$seqLenList)) {
    stop(
      "Missing TS position metadata on founderPop. Build founders with asMapPop() before setTrackTs().",
      call. = FALSE
    )
  }
  if (length(meta$posList) != founderPop@nChr || length(meta$seqLenList) != founderPop@nChr) {
    stop("TS position metadata length does not match number of chromosomes", call. = FALSE)
  }
  if ((!is.null(meta$breaksList) && length(meta$breaksList) != founderPop@nChr) ||
      (!is.null(meta$ratesList) && length(meta$ratesList) != founderPop@nChr)) {
    stop("TS rate-map metadata length does not match number of chromosomes", call. = FALSE)
  }
  meta
}

#' Initialize forward TS recorder state
#'
#' @param founderPop founder/sample population used to seed TS node maps
#' @return Recorder state object for `tsForwardAppendSeg()` / `tsForwardFinalize()`
#' @keywords internal
#' @noRd
tsForwardInit <- function(founderPop) {
  stopifnot(is(founderPop, "Pop") || is(founderPop, "MapPop"))

  nChr <- founderPop@nChr
  founderIid <- if (is(founderPop, "Pop")) {
    as.integer(founderPop@iid)
  } else {
    as.integer(seq_len(founderPop@nInd))
  }
  sourceByChr <- .tsForwardSourceByChr(founderPop)
  posMeta <- .tsForwardPosMetaFromFounder(founderPop)

  out <- list(
    chr = vector("list", nChr),
    ploidy = founderPop@ploidy,
    nChr = nChr,
    ts_path = lapply(sourceByChr, function(x) x[["ts_path"]]),
    posMeta = posMeta,
    indTime = numeric(0)
  )

  for (cc in seq_len(nChr)) {
    src <- sourceByChr[[cc]]
    tsObj <- asMapPop_load_ts(
      ts_path = src[["ts_path"]],
      ts = src[["ts"]],
      ts_xptr = src[["ts_xptr"]],
      tc_xptr = src[["tc_xptr"]],
      table_xptr = src[["table_xptr"]]
    )
    tc <- tsObj$dump_tables()
    samples <- tsObj$samples()
    expected <- founderPop@nInd * founderPop@ploidy
    if (length(samples) != expected) {
      stop("Sample node count does not match founderPop nInd*ploidy on chr ", cc, call. = FALSE)
    }

    nodeMap <- new.env(parent = emptyenv(), hash = TRUE)
    indMap <- new.env(parent = emptyenv(), hash = TRUE)

    idx <- 1L
    for (ind in seq_len(founderPop@nInd)) {
      iid <- founderIid[[ind]]
      indRow <- tc$node_table_get_row(as.integer(samples[[idx]]))$individual
      if (indRow < 0L) {
        stop("Sample node has individual = -1 on chr ", cc, call. = FALSE)
      }
      assign(as.character(iid), as.integer(indRow), envir = indMap)
      for (h in seq_len(founderPop@ploidy)) {
        key <- .tsForwardNodeKey(iid, h)
        assign(key, as.integer(samples[[idx]]), envir = nodeMap)
        idx <- idx + 1L
      }
    }

    out$chr[[cc]] <- list(
      tc = tc,
      nodeMap = nodeMap,
      indMap = indMap,
      seqLen = as.numeric(tc$sequence_length())
    )
  }

  class(out) <- "tsForwardRecorder"
  out
}

#' Append one segment table to forward TS recorder
#'
#' @param recorder output of `tsForwardInit()`
#' @param seg segment data frame with `tsSegGen`-style columns
#' @param simParam `SimParam` with pedigree
#' @return Updated recorder state
#' @keywords internal
#' @noRd
tsForwardAppendSeg <- function(recorder, seg, simParam) {
  stopifnot(inherits(recorder, "tsForwardRecorder"))

  if (is.null(seg) || nrow(seg) == 0L) {
    return(recorder)
  }

  requiredCols <- c("chr", "hap", "childIid", "parentIid", "parentHap", "left", "right")
  miss <- setdiff(requiredCols, names(seg))
  if (length(miss) > 0L) {
    stop("tsSegGen missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  recorder$indTime <- .tsForwardUpdateIndTimeCache(recorder$indTime, simParam$pedigree)
  indTime <- recorder$indTime

  for (cc in seq_len(recorder$nChr)) {
    chrState <- recorder$chr[[cc]]
    df <- seg[seg$chr == cc, , drop = FALSE]
    if (nrow(df) == 0L) {
      next
    }

    children <- unique(as.integer(df$childIid))
    .tsForwardEnsureChildrenIndividuals(
      chrState = chrState,
      children = children,
      simParam = simParam,
      cc = cc
    )

    childNodeKeys <- unique(.tsForwardNodeKey(df$childIid, df$hap))
    missingKeys <- childNodeKeys[!vapply(childNodeKeys, exists, logical(1),
                                         envir = chrState$nodeMap, inherits = FALSE)]
    if (length(missingKeys) > 0L) {
      keyInfo <- .tsForwardParseNodeKeys(missingKeys)
      nodeTimes <- as.numeric(indTime[keyInfo$iids])
      .tsForwardAddMissingChildNodes(chrState, missingKeys, nodeTimes)
    }

    pKeys <- .tsForwardNodeKey(df$parentIid, df$parentHap)
    cKeys <- .tsForwardNodeKey(df$childIid, df$hap)
    .tsForwardAppendEdgesByKeys(
      chrState = chrState,
      left = df$left,
      right = df$right,
      parentKeys = pKeys,
      childKeys = cKeys,
      cc = cc
    )
  }

  recorder
}

.tsForwardNormalizeSamplesByChr <- function(samples, nChr) {
  if (is.null(samples)) {
    return(NULL)
  }
  if (is.atomic(samples) && !is.list(samples)) {
    one <- as.integer(unique(samples))
    one <- one[is.finite(one) & !is.na(one)]
    return(rep(list(one), nChr))
  }
  if (!is.list(samples)) {
    stop("samples must be NULL, an integer vector, or a list of integer vectors", call. = FALSE)
  }
  if (length(samples) == 1L && nChr > 1L) {
    one <- as.integer(unique(samples[[1L]]))
    one <- one[is.finite(one) & !is.na(one)]
    return(rep(list(one), nChr))
  }
  if (length(samples) != nChr) {
    stop("samples list length must equal number of chromosomes (", nChr, ")", call. = FALSE)
  }
  out <- vector("list", nChr)
  for (cc in seq_len(nChr)) {
    x <- as.integer(unique(samples[[cc]]))
    x <- x[is.finite(x) & !is.na(x)]
    out[[cc]] <- x
  }
  out
}

.tsForwardDefaultSamplesByChr <- function(recorder, indTime) {
  has_ind_time <- !(is.null(indTime) || length(indTime) == 0L || !any(is.finite(indTime)))
  lastIid <- integer(0)
  if (has_ind_time) {
    lastTime <- min(indTime[is.finite(indTime)])
    lastIid <- which(indTime == lastTime)
  }

  out <- vector("list", recorder$nChr)
  for (cc in seq_len(recorder$nChr)) {
    nodeMap <- recorder$chr[[cc]]$nodeMap
    nNode <- as.integer(recorder$chr[[cc]]$tc$num_nodes())

    if (length(lastIid) > 0L) {
      keys <- as.vector(outer(lastIid, seq_len(recorder$ploidy), paste, sep = "_"))
      keep <- vapply(keys, exists, logical(1), envir = nodeMap, inherits = FALSE)
      keys <- keys[keep]
      if (length(keys) > 0L) {
        ids <- as.integer(vapply(keys, function(k) {
          get(k, envir = nodeMap, inherits = FALSE)
        }, integer(1)))
        ids <- ids[ids >= 0L & ids < nNode]
        if (length(ids) > 0L) {
          out[[cc]] <- sort(unique(ids))
          next
        }
      }
    }

    # Fallback: derive "last generation" from node times in recorder maps.
    allKeys <- ls(nodeMap, all.names = TRUE)
    if (length(allKeys) == 0L) {
      stop("Cannot infer default samples: empty node map on chromosome ", cc, call. = FALSE)
    }
    nodeIds <- as.integer(vapply(allKeys, function(k) {
      get(k, envir = nodeMap, inherits = FALSE)
    }, integer(1)))
    nodeIds <- nodeIds[nodeIds >= 0L & nodeIds < nNode]
    nodeIds <- unique(nodeIds)
    if (length(nodeIds) == 0L) {
      stop("Cannot infer default samples: no valid node ids in node map on chromosome ", cc, call. = FALSE)
    }
    times <- as.numeric(vapply(nodeIds, function(nid) {
      recorder$chr[[cc]]$tc$node_table_get_row(as.integer(nid))$time
    }, numeric(1)))
    ok <- is.finite(times)
    if (!any(ok)) {
      stop("Cannot infer default samples: no finite node times on chromosome ", cc, call. = FALSE)
    }
    tmin <- min(times[ok])
    out[[cc]] <- sort(unique(nodeIds[ok & times == tmin]))
    if (length(out[[cc]]) == 0L) {
      stop("Cannot infer default samples from node times on chromosome ", cc, call. = FALSE)
    }
  }
  out
}

.tsForwardResolveSamplesByChr <- function(recorder, samples, indTime) {
  out <- .tsForwardNormalizeSamplesByChr(samples, recorder$nChr)
  if (is.null(out)) {
    out <- .tsForwardDefaultSamplesByChr(recorder, indTime)
  }
  for (cc in seq_len(recorder$nChr)) {
    x <- as.integer(out[[cc]])
    if (length(x) == 0L) {
      stop("samples for chromosome ", cc, " cannot be empty", call. = FALSE)
    }
    nNode <- as.integer(recorder$chr[[cc]]$tc$num_nodes())
    if (any(x < 0L | x >= nNode)) {
      stop("samples out of bounds on chromosome ", cc, call. = FALSE)
    }
    out[[cc]] <- sort(unique(x))
  }
  out
}

.tsForwardCurrentSampleNodes <- function(tc) {
  nNode <- as.integer(tc$num_nodes())
  if (nNode <= 0L) {
    return(integer(0))
  }
  out <- integer(0)
  for (i in seq_len(nNode) - 1L) {
    fl <- as.integer(tc$node_table_get_row(i)$flags)
    if (isTRUE(bitwAnd(fl, 1L) != 0L)) {
      out <- c(out, i)
    }
  }
  out
}

.tsForwardPastIndividualNodes <- function(tc) {
  nNode <- as.integer(tc$num_nodes())
  if (nNode <= 0L) {
    return(integer(0))
  }
  out <- integer(0)
  for (i in seq_len(nNode) - 1L) {
    node <- tc$node_table_get_row(i)
    t <- as.numeric(node$time)
    ind <- as.integer(node$individual)
    if (is.finite(t) && t > 0 && !is.na(ind) && ind != -1L) {
      out <- c(out, i)
    }
  }
  out
}

.tsForwardOutputDirForChr <- function(recorder, cc, out_dir = NULL) {
  if (!is.null(out_dir)) {
    return(out_dir)
  }
  ts_path <- recorder$ts_path[[cc]]
  if (!is.null(ts_path) && is.character(ts_path) && nzchar(ts_path)) {
    return(dirname(ts_path))
  }
  getwd()
}

#' Finalize forward recorder into per-chromosome TreeSequence objects
#'
#' @param recorder output of `tsForwardInit()`
#' @param simplify logical; if TRUE, call table-collection simplify before dump
#' @param keep_unary logical; passed to `tc$simplify(keep_unary=...)` when `simplify=TRUE`.
#'   Default is TRUE.
#' @param update_sample_flags logical; passed to `tc$simplify(update_sample_flags=...)`
#'   when `simplify=TRUE`. Default is FALSE.
#' @param samples NULL (default: last generation nodes), integer vector (all chromosomes),
#'   or list of integer vectors by chromosome
#' @param keep_existing_samples logical; if TRUE (default), keep founder/original sample
#'   flags and add requested samples on top. If FALSE, replace sample flags with requested samples.
#' @param indTime optional individual-time vector; used only for default sample inference
#' @param update_samples logical; if TRUE apply sample selection/flags before sort/simplify.
#'   Default is FALSE for raw finalize behavior.
#' @return list of `RcppTskit::TreeSequence` objects (one per chromosome)
#' @keywords internal
#' @noRd
tsForwardFinalize <- function(recorder,
                              simplify = FALSE,
                              keep_unary = TRUE,
                              update_sample_flags = FALSE,
                              samples = NULL,
                              keep_existing_samples = TRUE,
                              indTime = NULL,
                              update_samples = FALSE) {
  stopifnot(inherits(recorder, "tsForwardRecorder"))
  stopifnot(is.logical(simplify), length(simplify) == 1L, !is.na(simplify))
  stopifnot(is.logical(keep_unary), length(keep_unary) == 1L, !is.na(keep_unary))
  stopifnot(is.logical(update_sample_flags), length(update_sample_flags) == 1L, !is.na(update_sample_flags))
  stopifnot(is.logical(keep_existing_samples), length(keep_existing_samples) == 1L, !is.na(keep_existing_samples))
  stopifnot(is.logical(update_samples), length(update_samples) == 1L, !is.na(update_samples))

  useDefaultSamples <- is.null(samples)
  samplesByChr <- NULL
  if (isTRUE(update_samples)) {
    samplesByChr <- .tsForwardResolveSamplesByChr(recorder, samples, indTime)
  }

  ts_list <- vector("list", recorder$nChr)
  for (cc in seq_len(recorder$nChr)) {
    chrState <- recorder$chr[[cc]]

    sampleNodes <- NULL
    if (isTRUE(update_samples)) {
      sampleNodes <- as.integer(samplesByChr[[cc]])
      if (isTRUE(simplify) && isTRUE(useDefaultSamples)) {
        sampleNodes <- c(
          .tsForwardCurrentSampleNodes(chrState$tc),
          sampleNodes
        )
      }
      if (isTRUE(keep_existing_samples)) {
        extraSampleNodes <- integer(0)
        if (isTRUE(simplify) && isTRUE(useDefaultSamples)) {
          extraSampleNodes <- c(
            extraSampleNodes,
            .tsForwardPastIndividualNodes(chrState$tc)
          )
        }
        sampleNodes <- sort(unique(c(extraSampleNodes, sampleNodes)))
      } else {
        sampleNodes <- sort(unique(sampleNodes))
      }
    }

    tcOut <- chrState$tc$clone()
    if (isTRUE(update_samples) && !isTRUE(simplify)) {
      tsForwardSetSampleFlags(tcOut, sampleNodes, clearExisting = !isTRUE(keep_existing_samples))
    }
    tcOut$sort()
    if (isTRUE(simplify)) {
      if (is.null(sampleNodes)) {
        sampleNodes <- .tsForwardCurrentSampleNodes(tcOut)
      }
      tcOut$simplify(
        samples = sampleNodes,
        keep_unary = keep_unary,
        update_sample_flags = update_sample_flags
      )
      tcOut$sort()
    }

    ts_list[[cc]] <- tcOut$tree_sequence()
  }
  ts_list
}

#' Convert table/TS objects to tree sequences and write all chromosomes
#'
#' @param ts_list list of `RcppTskit::TreeSequence`, table collections, or table pointers
#' @param out_dir output directory, default uses `recorder` source dirs
#' @param out_basename output basename prefix
#' @param recorder optional recorder for default output-directory resolution
#' @return character vector of written files
#' @keywords internal
#' @noRd
tsForwardWriteTreeSequences <- function(ts_list,
                                        out_dir = NULL,
                                        out_basename = "AlphaSimR_forward",
                                        recorder = NULL) {
  stopifnot(is.list(ts_list), length(ts_list) > 0L)
  nChr <- length(ts_list)
  if (is.null(out_dir) && is.null(recorder)) {
    stop("Provide out_dir or recorder for output-path resolution", call. = FALSE)
  }
  if (!is.null(recorder)) {
    stopifnot(inherits(recorder, "tsForwardRecorder"))
    if (recorder$nChr != nChr) {
      stop("recorder chromosome count does not match ts_list length", call. = FALSE)
    }
  }

  .as_ts <- function(x) {
    if (inherits(x, "externalptr")) {
      tc <- RcppTskit::TableCollection$new(xptr = x)
      return(tc$tree_sequence())
    }
    if (!is.null(x$dump) && is.function(x$dump)) {
      return(x)
    }
    if (!is.null(x$tree_sequence) && is.function(x$tree_sequence)) {
      return(x$tree_sequence())
    }
    stop("Each entry in ts_list must be a TreeSequence, table collection, or table pointer", call. = FALSE)
  }

  out_paths <- character(nChr)
  for (cc in seq_len(nChr)) {
    outDirCc <- if (!is.null(recorder)) {
      .tsForwardOutputDirForChr(recorder, cc, out_dir = out_dir)
    } else {
      out_dir
    }
    outPath <- file.path(outDirCc, paste0(out_basename, "_chr", cc - 1, ".trees"))
    ts_obj <- .as_ts(ts_list[[cc]])
    ts_obj$dump(outPath)
    out_paths[[cc]] <- outPath
  }
  out_paths
}

#' Get recorder from SimParam attribute
#' @keywords internal
#' @noRd
.tsForwardGetRecorder <- function(simParam) {
  rec <- attr(simParam, "tsForwardRecorder", exact = TRUE)
  if (is.null(rec)) {
    return(NULL)
  }
  if (!inherits(rec, "tsForwardRecorder")) {
    stop("simParam attr 'tsForwardRecorder' is not a tsForwardRecorder object", call. = FALSE)
  }
  rec
}

#' Set recorder on SimParam attribute
#' @keywords internal
#' @noRd
.tsForwardSetRecorder <- function(simParam, recorder) {
  if (!is.null(recorder) && !inherits(recorder, "tsForwardRecorder")) {
    stop("recorder must be NULL or a tsForwardRecorder object", call. = FALSE)
  }
  attr(simParam, "tsForwardRecorder") <- recorder
  invisible(simParam)
}

#' Check whether a forward recorder is attached to SimParam
#' @keywords internal
#' @noRd
tsForwardHasRecorder <- function(simParam) {
  !is.null(.tsForwardGetRecorder(simParam))
}

#' Initialize and attach a forward recorder to SimParam
#' @param simParam SimParam object
#' @param founderPop founder/sample population used to seed the recorder node map
#' @return invisibly returns simParam
#' @keywords internal
#' @noRd
tsForwardInitOnSimParam <- function(simParam, founderPop) {
  stopifnot(inherits(simParam, "SimParam"))
  stopifnot(is(founderPop, "Pop") || is(founderPop, "MapPop"))
  rec <- tsForwardInit(founderPop = founderPop)
  .tsForwardSetRecorder(simParam, rec)
}

.tsForwardFinalizeWriteFromSimParam <- function(simParam,
                                                out_dir = NULL,
                                                out_basename = "AlphaSimR_forward",
                                                simplify = FALSE,
                                                keep_unary = TRUE,
                                                update_sample_flags = FALSE,
                                                samples = NULL,
                                                keep_existing_samples = TRUE,
                                                clear = FALSE,
                                                update_samples = FALSE) {
  rec <- .tsForwardGetRecorder(simParam)
  if (is.null(rec)) {
    stop("No tsForwardRecorder attached to simParam", call. = FALSE)
  }
  indTime <- .tsForwardUpdateIndTimeCache(rec$indTime, simParam$pedigree)
  rec$indTime <- indTime
  ts_list <- tsForwardFinalize(
    recorder = rec,
    simplify = simplify,
    keep_unary = keep_unary,
    update_sample_flags = update_sample_flags,
    samples = samples,
    keep_existing_samples = keep_existing_samples,
    indTime = indTime,
    update_samples = update_samples
  )
  out_paths <- tsForwardWriteTreeSequences(
    ts_list = ts_list,
    out_dir = out_dir,
    out_basename = out_basename,
    recorder = rec
  )
  if (isTRUE(clear)) {
    .tsForwardSetRecorder(simParam, NULL)
  } else {
    .tsForwardSetRecorder(simParam, rec)
  }
  out_paths
}

#' Finalize recorder attached to SimParam and write trees
#' @param simParam SimParam object
#' @param out_dir output directory
#' @param out_basename output basename
#' @param simplify logical; if TRUE, call table-collection simplify before dump
#' @param keep_unary logical; passed to `tc$simplify(keep_unary=...)` when `simplify=TRUE`.
#' @param update_sample_flags logical; passed to `tc$simplify(update_sample_flags=...)`
#'   when `simplify=TRUE`. Default is FALSE.
#' @param samples optional sample set (used only when simplifying and/or when
#'   `keep_existing_samples = FALSE`)
#' @param keep_existing_samples logical; if TRUE (default), keep founder/original sample
#'   flags when sample updates are requested
#' @param clear clear recorder on simParam after writing
#' @return character vector of output tree file paths
#' @keywords internal
#' @noRd
tsForwardFinalizeFromSimParam <- function(simParam,
                                          out_dir = NULL,
                                          out_basename = "AlphaSimR_forward",
                                          simplify = FALSE,
                                          keep_unary = TRUE,
                                          update_sample_flags = FALSE,
                                          samples = NULL,
                                          keep_existing_samples = TRUE,
                                          clear = FALSE) {
  use_sample_updates <- !is.null(samples) || !isTRUE(keep_existing_samples) || isTRUE(simplify)
  .tsForwardFinalizeWriteFromSimParam(
    simParam = simParam,
    out_dir = out_dir,
    out_basename = out_basename,
    simplify = simplify,
    keep_unary = keep_unary,
    update_sample_flags = update_sample_flags,
    samples = samples,
    keep_existing_samples = keep_existing_samples,
    clear = clear,
    update_samples = use_sample_updates
  )
}

#' Write recorder attached to SimParam with custom samples and optional simplify
#' @param simParam SimParam object
#' @param out_dir output directory
#' @param out_basename output basename
#' @param simplify logical; if TRUE, call table-collection simplify before dump
#' @param keep_unary logical; passed to `tc$simplify(keep_unary=...)` when `simplify=TRUE`.
#'   Default is TRUE.
#' @param update_sample_flags logical; passed to `tc$simplify(update_sample_flags=...)`
#'   when `simplify=TRUE`. Default is FALSE.
#' @param samples NULL (default: last generation nodes), integer vector (all chromosomes),
#'   or list of integer vectors by chromosome
#' @param keep_existing_samples logical; if TRUE (default), keep founder/original sample
#'   flags and add requested samples on top.
#' @param clear clear recorder on simParam after writing
#' @return character vector of output tree file paths
#' @keywords internal
#' @noRd
tsForwardWriteTreesFromSimParam <- function(simParam,
                                            out_dir = NULL,
                                            out_basename = "AlphaSimR_forward",
                                            simplify = FALSE,
                                            keep_unary = TRUE,
                                            update_sample_flags = FALSE,
                                            samples = NULL,
                                            keep_existing_samples = TRUE,
                                            clear = FALSE) {
  .tsForwardFinalizeWriteFromSimParam(
    simParam = simParam,
    out_dir = out_dir,
    out_basename = out_basename,
    simplify = simplify,
    keep_unary = keep_unary,
    update_sample_flags = update_sample_flags,
    samples = samples,
    keep_existing_samples = keep_existing_samples,
    clear = clear,
    update_samples = TRUE
  )
}
