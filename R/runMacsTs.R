#' Parse Scaled Mutation Rate (`dTheta`) from a MaCS Argument String
#'
#' @param args Character scalar MaCS command string with sample size and
#'   sequence length as the first two tokens.
#'
#' @return Numeric scalar `dTheta` used by MaCS-style mutation placement.
#' @keywords internal
#' @noRd
.simAnc_parse_dTheta <- function(args) {
  tokens <- strsplit(as.character(args), "[,[:space:]]+", perl = TRUE)[[1L]]
  tokens <- tokens[nzchar(tokens)]
  if (length(tokens) < 2L) {
    stop("args must contain at least sample size and sequence length")
  }
  seqLen <- suppressWarnings(as.numeric(tokens[2L]))
  if (!is.finite(seqLen) || seqLen <= 0) {
    stop("Failed to parse sequence length from args")
  }
  idx <- match("-t", tokens)
  if (is.na(idx) || idx >= length(tokens)) {
    return(0)
  }
  thetaScaled <- suppressWarnings(as.numeric(tokens[idx + 1L]))
  if (!is.finite(thetaScaled) || thetaScaled < 0) {
    stop("Failed to parse -t value from args")
  }
  seqLen * thetaScaled
}

#' Simulate MaCS Ancestry as Tree-Sequence Tables
#'
#' @description
#' Runs the MaCS ancestry simulation path and returns one tskit table collection
#' per chromosome, without adding post-ancestry mutations. This is the staged
#' ancestry step used by \code{\link{runMacsTS}}.
#'
#' @param args Character MaCS command prefix with trailing \code{-s}.
#' @param nChr Integer number of chromosomes.
#' @param inbred Logical.
#' @param ploidy Integer ploidy.
#' @param nThreads Integer thread count.
#' @param seed Integer vector of chromosome seeds.
#' @param usePhysicalPositions Logical; use bp positions in TS if \code{TRUE}.
#' @param Nref Optional numeric reference \code{Ne} for time scaling.
#'
#' @details
#' This function is useful when ancestry simulation and mutation placement need
#' to be controlled separately. For the standard founder-population workflow,
#' use \code{\link{runMacsTS}}.
#'
#' @return List with ancestry table collections and metadata.
#'
#' @examples
#' \dontrun{
#' cmd = runMacs2(nInd=10, nChr=1, segSites=NULL, returnCommand=TRUE)
#' anc = simAnc(args=paste0(20, " ", cmd, " -s "),
#'              nChr=1, inbred=FALSE, ploidy=2L,
#'              nThreads=1L, seed=123L)
#' }
#'
#' @export
simAnc <- function(args, nChr, inbred, ploidy, nThreads, seed,
                   usePhysicalPositions = FALSE, Nref = NA_real_) {
  nChr <- as.integer(nChr)
  if (length(nChr) != 1L || is.na(nChr) || nChr <= 0L) {
    stop("nChr must be a positive integer scalar")
  }
  anc <- MaCSTS(
    args = args,
    nChr = nChr,
    inbred = inbred,
    ploidy = ploidy,
    nThreads = nThreads,
    seed = seed,
    usePhysicalPositions = usePhysicalPositions,
    useMacsMut = FALSE,
    Nref = Nref,
    expandInbredSamples = FALSE
  )
  anc$dTheta <- .simAnc_parse_dTheta(args)
  anc$seed <- as.integer(seed)
  anc$ploidy <- as.integer(ploidy)
  anc$inbred <- isTRUE(inbred)
  anc$stage <- "simAnc"
  anc
}

#' Add Mutations to Ancestry Tree-Sequence Tables
#'
#' @description
#' Adds MaCS-style mutations to table collections returned by
#' \code{\link{simAnc}}. This is the staged mutation step used by
#' \code{\link{runMacsTS}} when \code{mutationMode = "postTs"}.
#'
#' @param x List returned by \code{\link{simAnc}}, or a list of
#' table-collection pointers.
#' @param dTheta Optional scalar/vector mutation-rate parameter in MaCS units.
#' @param seed Optional scalar/vector integer seeds for mutation sampling.
#'
#' @details
#' If \code{dTheta} is omitted and \code{x} is the full list returned by
#' \code{\link{simAnc}}, the mutation-rate metadata stored by \code{simAnc()} is
#' used. If \code{seed} is omitted, mutation seeds are derived from the ancestry
#' seeds when available.
#'
#' @return List with mutated table collections and metadata.
#'
#' @examples
#' \dontrun{
#' cmd = runMacs2(nInd=10, nChr=1, segSites=NULL, returnCommand=TRUE)
#' anc = simAnc(args=paste0(20, " ", cmd, " -s "),
#'              nChr=1, inbred=FALSE, ploidy=2L,
#'              nThreads=1L, seed=123L)
#' mut = simMut(anc, seed=456L)
#' }
#'
#' @export
simMut <- function(x, dTheta = NULL, seed = NULL) {
  tables <- if (is.list(x) && !is.null(x$tables)) x$tables else x
  if (!is.list(tables) || length(tables) == 0L) {
    stop("simMut requires a non-empty list of table collections")
  }
  nChr <- length(tables)
  
  if (is.null(dTheta)) {
    if (is.list(x) && !is.null(x$dTheta)) {
      dTheta <- x$dTheta
    } else {
      stop("dTheta is required when x has no dTheta metadata")
    }
  }
  if (length(dTheta) == 1L) {
    dTheta <- rep(as.numeric(dTheta), nChr)
  }
  if (length(dTheta) != nChr) {
    stop("dTheta length must be 1 or number of chromosomes")
  }
  
  if (is.null(seed)) {
    if (is.list(x) && !is.null(x$seed)) {
      seed <- as.integer(x$seed) + 104729L
    } else {
      seed <- sample.int(1e8, nChr)
    }
  }
  if (length(seed) == 1L) {
    seed <- rep(as.integer(seed), nChr)
  }
  if (length(seed) != nChr) {
    stop("seed length must be 1 or number of chromosomes")
  }
  
  for (chr in seq_len(nChr)) {
    tsMutateTableCollection(tables[[chr]], as.numeric(dTheta[[chr]]), as.numeric(seed[[chr]]))
  }
  
  out <- if (is.list(x) && !is.null(x$tables)) x else list()
  out$tables <- tables
  out$dTheta <- as.numeric(dTheta)
  out$mutationSeed <- as.integer(seed)
  out$mutationMode <- "postTs"
  out$stage <- "simMut"
  out
}

#' Finalize Inbred TS by Duplicating Sample Leaves per Individual
#'
#' @param x List returned by `simAnc`/`simMut`, or a list of table pointers.
#' @param inbred Logical.
#' @param ploidy Integer ploidy.
#'
#' @return List with finalized table collections and metadata.
#' @keywords internal
#' @noRd
finalizeInbredTs <- function(x, inbred = FALSE, ploidy = 2L) {
  tables <- if (is.list(x) && !is.null(x$tables)) x$tables else x
  if (!is.list(tables) || length(tables) == 0L) {
    stop("finalizeInbredTs requires a non-empty list of table collections")
  }
  ploidy <- as.integer(ploidy)
  if (ploidy <= 0L) {
    stop("ploidy must be a positive integer")
  }
  
  if (isTRUE(inbred) && ploidy > 1L) {
    for (chr in seq_along(tables)) {
      tsFinalizeInbredTableCollection(tables[[chr]], ploidy)
    }
  }
  
  out <- if (is.list(x) && !is.null(x$tables)) x else list()
  out$tables <- tables
  out$inbred <- isTRUE(inbred)
  out$ploidy <- ploidy
  out$stage <- "finalizeInbredTs"
  out
}

#' Build runMacs-style MaCS Command for TS Workflow
#'
#' @param nInd Integer number of individuals.
#' @param inbred Logical.
#' @param species Character species preset name.
#' @param split Optional split time in generations.
#' @param ploidy Integer ploidy.
#' @param manualCommand Optional user-provided MaCS command tail.
#' @param manualGenLen Optional user-provided chromosome genetic length(s) in Morgan.
#' @param nChr Integer number of chromosomes.
#'
#' @return List with `command`, `genLen`, and `seqLen`.
#' @keywords internal
#' @noRd
.runMacsTS_build_command <- function(nInd, inbred, species, split, ploidy,
                                    manualCommand, manualGenLen, nChr) {
  popSize <- ifelse(inbred, nInd, ploidy * nInd)
  if (!is.null(manualCommand)) {
    if (is.null(manualGenLen)) {
      stop("You must define manualGenLen when using manualCommand")
    }
    command <- paste0(popSize, " ", manualCommand, " -s ")
    genLen <- manualGenLen
  } else {
    species <- toupper(species)
    if (species == "GENERIC") {
      genLen <- 1.0
      Ne <- 100
      speciesParams <- "1E8 -t 1E-5 -r 4E-6"
      speciesHist <- "-eN 0.25 5.0 -eN 2.50 15.0 -eN 25.00 60.0 -eN 250.00 120.0 -eN 2500.00 1000.0"
    } else if (species == "CATTLE") {
      cattleChrSum <- 2.8e9
      cattleChrBp <- cattleChrSum / 30
      recRate <- 9.26e-9
      genLen <- recRate * cattleChrBp
      mutRate <- 9.4e-9
      Ne <- 90
      histNe <- c(120, 250, 350, 1000, 1500, 2000, 2500, 3500, 7000, 10000, 17000, 62000)
      histGen <- c(3, 6, 12, 18, 24, 154, 454, 654, 1754, 2354, 3354, 33154)
      speciesParams <- paste(c(round(cattleChrBp), "-t", mutRate * 4 * Ne, "-r", recRate * 4 * Ne),
                             collapse = " ")
      histNe <- histNe / Ne
      histGen <- histGen / (4 * Ne)
      speciesHist <- NULL
      for (i in seq_len(length(histNe))) {
        speciesHist <- paste(speciesHist, "-eN", histGen[i], histNe[i])
      }
    } else if (species == "WHEAT") {
      genLen <- 1.43
      Ne <- 50
      speciesParams <- "8E8 -t 4E-7 -r 3.6E-7"
      speciesHist <- "-eN 0.03 1 -eN 0.05 2 -eN 0.10 4 -eN 0.15 6 -eN 0.20 8 -eN 0.25 10 -eN 0.30 12 -eN 0.35 14 -eN 0.40 16 -eN 0.45 18 -eN 0.50 20 -eN 1.00 40 -eN 2.00 60 -eN 3.00 80 -eN 4.00 100 -eN 5.00 120 -eN 10.00 140 -eN 20.00 160 -eN 30.00 180 -eN 40.00 200 -eN 50.00 240 -eN 100.00 320 -eN 200.00 400 -eN 300.00 480 -eN 400.00 560 -eN 500.00 640"
    } else if (species == "MAIZE") {
      genLen <- 2.0
      Ne <- 100
      speciesParams <- "2E8 -t 5E-6 -r 4E-6"
      speciesHist <- "-eN 0.03 1 -eN 0.05 2 -eN 0.10 4 -eN 0.15 6 -eN 0.20 8 -eN 0.25 10 -eN 0.30 12 -eN 0.35 14 -eN 0.40 16 -eN 0.45 18 -eN 0.50 20 -eN 2.00 40 -eN 3.00 60 -eN 4.00 80 -eN 5.00 100"
    } else {
      stop("No rules for species ", species)
    }
    if (is.null(split)) {
      splitI <- ""
      splitJ <- ""
    } else {
      stopifnot(popSize %% 2 == 0)
      splitI <- paste(" -I 2", popSize %/% 2, popSize %/% 2)
      splitJ <- paste(" -ej", split / (4 * Ne) + 0.000001, "2 1")
    }
    command <- paste0(popSize, " ", speciesParams, splitI, " ", speciesHist, splitJ, " -s ")
  }
  if (!is.null(manualGenLen)) {
    genLen <- manualGenLen
  }
  if (length(genLen) == 1L) {
    genLen <- rep(genLen, nChr)
  }
  if (length(genLen) != nChr) {
    stop("genLen must have length 1 or nChr")
  }
  tokens <- strsplit(command, "[,[:space:]]+", perl = TRUE)[[1L]]
  tokens <- tokens[nzchar(tokens)]
  if (length(tokens) < 2L) {
    stop("Failed to parse sequence length from command")
  }
  seqLen <- suppressWarnings(as.numeric(tokens[2L]))
  if (!is.finite(seqLen) || seqLen <= 0) {
    stop("Invalid sequence length parsed from command")
  }
  list(command = command, genLen = as.numeric(genLen), seqLen = seqLen)
}

.runMacsTS_hotspot_path <- function(args) {
  tokens <- strsplit(as.character(args), "[,[:space:]]+", perl = TRUE)[[1L]]
  tokens <- tokens[nzchar(tokens)]
  idx <- match("-R", tokens)
  if (is.na(idx) || idx >= length(tokens)) {
    return(NULL)
  }
  tokens[[idx + 1L]]
}

.runMacsTS_read_hotspot_map <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("Invalid MaCS -R hotspot file path", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("MaCS -R hotspot file does not exist: ", path, call. = FALSE)
  }
  hot <- utils::read.table(path, header = FALSE, col.names = c("start", "end", "ratio"))
  if (ncol(hot) != 3L || nrow(hot) == 0L) {
    stop("MaCS -R hotspot file must contain rows of: start end ratio", call. = FALSE)
  }
  hot$start <- as.numeric(hot$start)
  hot$end <- as.numeric(hot$end)
  hot$ratio <- as.numeric(hot$ratio)
  bad <- !is.finite(hot$start) | !is.finite(hot$end) | !is.finite(hot$ratio) |
    hot$start < 0 | hot$end > 1 | hot$start >= hot$end | hot$ratio < 0
  if (any(bad)) {
    stop("Invalid MaCS -R hotspot row. Expected 0 <= start < end <= 1 and ratio >= 0.", call. = FALSE)
  }
  hot <- hot[order(hot$start, hot$end), , drop = FALSE]
  if (nrow(hot) > 1L && any(hot$start[-1L] < hot$end[-nrow(hot)])) {
    stop("Overlapping MaCS -R hotspot intervals are not supported", call. = FALSE)
  }
  hot
}

.runMacsTS_map_from_hotspots <- function(path, nChr, seqLen, genLen, usePhysicalPositions) {
  hot <- .runMacsTS_read_hotspot_map(path)
  coordLen <- if (isTRUE(usePhysicalPositions)) as.numeric(seqLen) else 1
  starts <- coordLen * hot$start
  ends <- coordLen * hot$end
  recBreaks <- sort(unique(c(0, coordLen, starts, ends)))
  ratio <- rep(1, length(recBreaks) - 1L)
  for (i in seq_len(nrow(hot))) {
    idx <- recBreaks[-length(recBreaks)] >= starts[[i]] & recBreaks[-1L] <= ends[[i]]
    ratio[idx] <- hot$ratio[[i]]
  }
  recRates <- lapply(as.numeric(genLen), function(g) {
    (g / coordLen) * ratio
  })
  list(
    breaks = rep(list(recBreaks), as.integer(nChr)),
    rates = recRates
  )
}

.runMacsTS_resolve_rec_map <- function(args, nChr, seqLen, genLen, usePhysicalPositions) {
  hotspotPath <- .runMacsTS_hotspot_path(args)
  if (!is.null(hotspotPath)) {
    return(.runMacsTS_map_from_hotspots(
      path = hotspotPath,
      nChr = nChr,
      seqLen = seqLen,
      genLen = genLen,
      usePhysicalPositions = usePhysicalPositions
    ))
  }

  recBreaks <- if (usePhysicalPositions) {
    rep(list(c(0, seqLen)), nChr)
  } else {
    rep(list(c(0, 1)), nChr)
  }
  recRates <- if (usePhysicalPositions) {
    lapply(genLen, function(g) c(g / seqLen))
  } else {
    lapply(genLen, function(g) c(g))
  }
  list(breaks = recBreaks, rates = recRates)
}

#' Create Founder Haplotypes and Tree Sequences using MaCS
#'
#' @description
#' Tree-sequence analogue of \code{\link{runMacs}}. This function simulates
#' MaCS ancestry, optionally adds mutations, converts the resulting tree
#' sequence data to a \code{\link{MapPop-class}}, and can return the founder
#' table collections for forward tree-sequence recording.
#'
#' @param nInd Integer number of individuals to simulate.
#' @param nChr Integer number of chromosomes.
#' @param segSites Optional site-count cap per chromosome (scalar or vector).
#' @param inbred Logical.
#' @param species Species preset used by \code{\link{runMacs}}.
#' @param split Optional population split time in generations.
#' @param ploidy Integer ploidy.
#' @param manualCommand Optional MaCS command tail (advanced users).
#' @param manualGenLen Optional genetic length(s) in Morgan.
#' @param nThreads Optional thread count.
#' @param mutationMode One of \code{"postTs"}, \code{"macs"}, or
#' \code{"none"}.
#' @param usePhysicalPositions Logical; TS coordinates in bp if \code{TRUE}.
#' @param Nref Optional reference \code{Ne} for time scaling.
#' @param seed Optional integer vector (length 1 or \code{nChr}) for ancestry.
#' @param mutSeed Optional integer vector (length 1 or \code{nChr}) for
#' post-TS mutation.
#' @param mutSeedOffset Integer offset used when deriving post-TS mutation seeds.
#' @param siteSamplingSeed Integer seed for \code{\link{asMapPop}} site
#' sampling.
#' @param expandInbredTs Logical; whether to expand inbred TS sample leaves before conversion.
#' @param returnTs Logical; return TS tables and metadata alongside
#' \code{\link{MapPop-class}}.
#'
#' @details
#' \code{mutationMode = "postTs"} first creates ancestry tables and then calls
#' \code{\link{simMut}} to place mutations. \code{mutationMode = "macs"} keeps
#' the mutation placement from the native MaCS path. \code{mutationMode = "none"}
#' returns ancestry-only tables and requires \code{returnTs = TRUE}.
#'
#' When \code{returnTs = TRUE}, the returned population keeps
#' \code{tsForwardSource} and \code{tsForwardPosMeta} attributes that can be used
#' by \code{SimParam$setTrackTs(TRUE, founderPop=...)}.
#'
#' @return \code{\link{MapPop-class}} by default; otherwise a list with
#' \code{pop}, \code{tables}, and metadata.
#'
#' @examples
#' \dontrun{
#' founderPop = runMacsTS(nInd=10, nChr=1, segSites=100)
#'
#' out = runMacsTS(nInd=10, nChr=1, segSites=100, returnTs=TRUE)
#' SP = SimParam$new(out$pop)
#' SP$setTrackTs(TRUE, founderPop=out$pop)
#' paths = tsForwardFinalizeFromSimParam(SP, out_dir=tempdir())
#' }
#'
#' @export
runMacsTS <- function(nInd, nChr = 1, segSites = NULL, inbred = FALSE,
                     species = "GENERIC", split = NULL, ploidy = 2L,
                     manualCommand = NULL, manualGenLen = NULL, nThreads = NULL,
                     mutationMode = c("postTs", "macs", "none"),
                     usePhysicalPositions = FALSE, Nref = NA_real_,
                     seed = NULL, mutSeed = NULL, mutSeedOffset = 104729L,
                     siteSamplingSeed = 42L, expandInbredTs = FALSE,
                     returnTs = FALSE) {
  mutationMode <- match.arg(mutationMode)
  nInd <- as.integer(nInd)
  nChr <- as.integer(nChr)
  ploidy <- as.integer(ploidy)
  if (is.null(nThreads)) {
    nThreads <- getNumThreads()
  }
  nThreads <- as.integer(nThreads)
  if (nChr < nThreads) {
    nThreads <- nChr
  }
  if (nInd <= 0L || nChr <= 0L || ploidy <= 0L) {
    stop("nInd, nChr, and ploidy must be positive integers")
  }
  if (!is.null(segSites)) {
    segSites <- as.integer(segSites)
    if (length(segSites) == 1L) {
      segSites <- rep(segSites, nChr)
    }
    if (length(segSites) != nChr) {
      stop("segSites must have length 1 or nChr")
    }
  }
  
  setup <- .runMacsTS_build_command(
    nInd = nInd,
    inbred = inbred,
    species = species,
    split = split,
    ploidy = ploidy,
    manualCommand = manualCommand,
    manualGenLen = manualGenLen,
    nChr = nChr
  )
  args <- setup$command
  genLen <- setup$genLen
  seqLen <- setup$seqLen
  
  if (is.null(seed)) {
    seed <- sample.int(n = 1e8, size = nChr)
  }
  seed <- as.integer(seed)
  if (length(seed) == 1L) {
    seed <- rep(seed, nChr)
  }
  if (length(seed) != nChr) {
    stop("seed must have length 1 or nChr")
  }
  
  runOut <- NULL
  if (mutationMode == "macs") {
    runOut <- MaCSTS(
      args = args,
      nChr = nChr,
      inbred = inbred,
      ploidy = ploidy,
      nThreads = nThreads,
      seed = seed,
      usePhysicalPositions = usePhysicalPositions,
      useMacsMut = TRUE,
      Nref = Nref,
      expandInbredSamples = FALSE
    )
  } else {
    runOut <- simAnc(
      args = args,
      nChr = nChr,
      inbred = inbred,
      ploidy = ploidy,
      nThreads = nThreads,
      seed = seed,
      usePhysicalPositions = usePhysicalPositions,
      Nref = Nref
    )
    if (mutationMode == "postTs") {
      if (is.null(mutSeed)) {
        mutSeed <- as.integer(seed + as.integer(mutSeedOffset))
      }
      mutSeed <- as.integer(mutSeed)
      if (length(mutSeed) == 1L) {
        mutSeed <- rep(mutSeed, nChr)
      }
      if (length(mutSeed) != nChr) {
        stop("mutSeed must have length 1 or nChr")
      }
      timeScale <- if (!is.null(runOut$timeScale)) as.numeric(runOut$timeScale) else 1
      dThetaPost <- as.numeric(runOut$dTheta) / timeScale
      runOut <- simMut(runOut, dTheta = dThetaPost, seed = mutSeed)
    }
  }
  
  if (isTRUE(expandInbredTs) && isTRUE(inbred) && ploidy > 1L) {
    runOut <- finalizeInbredTs(runOut, inbred = inbred, ploidy = ploidy)
  }
  
  if (mutationMode == "none") {
    if (!isTRUE(returnTs)) {
      stop("mutationMode='none' produces ancestry-only TS with zero sites; set returnTs=TRUE or use mutationMode='postTs'/'macs'.")
    }
    return(list(
      pop = NULL,
      tables = runOut$tables,
      args = args,
      seed = seed,
      mutationMode = mutationMode,
      mutSeed = NA_integer_,
      usePhysicalPositions = usePhysicalPositions,
      timeScale = if (!is.null(runOut$timeScale)) runOut$timeScale else 1,
      Nref = if (!is.null(runOut$Nref)) runOut$Nref else NA_real_
    ))
  }
  
  siteCounts <- vapply(runOut$tables, function(tc_xptr) {
    as.integer(rtsk_table_collection_summary2(tc_xptr)$num_sites)
  }, integer(1))
  if (any(siteCounts <= 0L)) {
    badChr <- which(siteCounts <= 0L)
    stop("No segregating sites on chromosome(s): ",
         paste(badChr, collapse = ", "),
         ". Increase mutation rate or inspect TS via returnTs=TRUE.")
  }
  
  recMap <- .runMacsTS_resolve_rec_map(
    args = args,
    nChr = nChr,
    seqLen = seqLen,
    genLen = genLen,
    usePhysicalPositions = usePhysicalPositions
  )
  
  popOut <- asMapPop(
    chr_info = list(
      tables = runOut$tables,
      breaks = recMap$breaks,
      rates = recMap$rates
    ),
    ploidy = ploidy,
    inbred = inbred,
    segSites = segSites,
    site_sampling_seed = as.integer(siteSamplingSeed),
    nThreads = nThreads,
    returnMeta = FALSE
  )
  
  if (!isTRUE(returnTs)) {
    return(popOut)
  }
  list(
    pop = popOut,
    tables = runOut$tables,
    args = args,
    seed = seed,
    mutationMode = mutationMode,
    mutSeed = if (!is.null(runOut$mutationSeed)) runOut$mutationSeed else NA_integer_,
    usePhysicalPositions = usePhysicalPositions,
    timeScale = if (!is.null(runOut$timeScale)) runOut$timeScale else 1,
    Nref = if (!is.null(runOut$Nref)) runOut$Nref else NA_real_
  )
}

#' Build file-backed bridge chromosome info from a runMacsTS result
#'
#' `runMacsTS(..., returnTs = TRUE)` keeps founder tables in memory, while the
#' old bridge validation path expects `chr_info[[cc]]$ts_path`. This helper
#' writes those tables to `.trees` files and reuses the map metadata stored by
#' `asMapPop()`.
#'
#' @param x List returned by `runMacsTS(..., returnTs = TRUE)`.
#' @param out_dir Directory where founder `.trees` files should be written.
#' @param out_basename Prefix for written founder tree files.
#' @param segSites Optional scalar or per-chromosome site counts for returned
#'   `chr_info`. Defaults to `x$pop@nLoci`.
#'
#' @return List of chromosome info entries with `ts_path`, `breaks`, `rates`,
#'   and `segSites`.
#' @keywords internal
#' @noRd
runMacsTSBridgeChrInfo <- function(x, out_dir, out_basename = "runMacsTS_founder",
                                  segSites = NULL) {
  if (!is.list(x) || is.null(x$pop) || is.null(x$tables)) {
    stop("x must be the list returned by runMacsTS(..., returnTs = TRUE)", call. = FALSE)
  }
  pop <- x$pop
  tables <- x$tables
  if (!is(pop, "MapPop") && !is(pop, "Pop")) {
    stop("x$pop must be a Pop or MapPop object", call. = FALSE)
  }
  if (!is.list(tables) || length(tables) == 0L) {
    stop("x$tables must be a non-empty list of table collection pointers", call. = FALSE)
  }
  nChr <- length(tables)
  if (pop@nChr != nChr) {
    stop("length(x$tables) must match x$pop@nChr", call. = FALSE)
  }
  if (!is.character(out_dir) || length(out_dir) != 1L || is.na(out_dir) || !nzchar(out_dir)) {
    stop("out_dir must be a non-empty character scalar", call. = FALSE)
  }
  if (!is.character(out_basename) || length(out_basename) != 1L ||
      is.na(out_basename) || !nzchar(out_basename)) {
    stop("out_basename must be a non-empty character scalar", call. = FALSE)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(out_dir)) {
    stop("Failed to create out_dir: ", out_dir, call. = FALSE)
  }

  posMeta <- attr(pop, "tsForwardPosMeta", exact = TRUE)
  if (is.null(posMeta) || is.null(posMeta$posList) || length(posMeta$posList) != nChr ||
      is.null(posMeta$breaksList) || length(posMeta$breaksList) != nChr ||
      is.null(posMeta$ratesList) || length(posMeta$ratesList) != nChr) {
    stop("x$pop is missing TS position metadata from asMapPop()", call. = FALSE)
  }
  if (length(pop@genMap) != nChr) {
    stop("x$pop@genMap length must match x$pop@nChr", call. = FALSE)
  }

  if (is.null(segSites)) {
    segSites <- as.integer(pop@nLoci)
  } else {
    segSites <- as.integer(segSites)
  }
  if (length(segSites) == 1L) {
    segSites <- rep(segSites, nChr)
  }
  if (length(segSites) != nChr || any(is.na(segSites)) || any(segSites < 0L)) {
    stop("segSites must have length 1 or number of chromosomes", call. = FALSE)
  }

  lapply(seq_len(nChr), function(cc) {
    tc <- RcppTskit::TableCollection$new(xptr = tables[[cc]])
    seqLen <- as.numeric(tc$sequence_length())
    if (!is.finite(seqLen) || seqLen <= 0) {
      stop("Invalid sequence length for chromosome ", cc, call. = FALSE)
    }

    tsPath <- file.path(out_dir, paste0(out_basename, "_chr", cc - 1L, ".trees"))
    tc$tree_sequence()$dump(tsPath)

    list(
      ts_path = tsPath,
      breaks = as.numeric(posMeta$breaksList[[cc]]),
      rates = as.numeric(posMeta$ratesList[[cc]]),
      segSites = segSites[[cc]]
    )
  })
}

# Backward-compatible development aliases for the earlier runMacTS spelling.
runMacTS <- function(...) {
  runMacsTS(...)
}

runMacTSBridgeChrInfo <- function(x, out_dir, out_basename = "runMacTS_founder",
                                  segSites = NULL) {
  runMacsTSBridgeChrInfo(x=x,
                         out_dir=out_dir,
                         out_basename=out_basename,
                         segSites=segSites)
}

.runMacTS_build_command <- .runMacsTS_build_command
.runMacTS_hotspot_path <- .runMacsTS_hotspot_path
.runMacTS_read_hotspot_map <- .runMacsTS_read_hotspot_map
.runMacTS_map_from_hotspots <- .runMacsTS_map_from_hotspots
.runMacTS_resolve_rec_map <- .runMacsTS_resolve_rec_map
