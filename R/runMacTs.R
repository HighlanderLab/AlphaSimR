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

#' Simulate Ancestry as TS Tables without Post-Ancestry Mutations
#'
#' @param args Character MaCS command prefix with trailing `-s`.
#' @param nChr Integer number of chromosomes.
#' @param inbred Logical.
#' @param ploidy Integer ploidy.
#' @param nThreads Integer thread count.
#' @param seed Integer vector of chromosome seeds.
#' @param usePhysicalPositions Logical; use bp positions in TS if `TRUE`.
#' @param Nref Optional numeric reference `Ne` for time scaling.
#'
#' @return List with ancestry table collections and metadata.
#' @keywords internal
#' @noRd
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

#' Add Mutations to Ancestry TS Tables
#'
#' @param x List returned by `simAnc`, or a list of table-collection pointers.
#' @param dTheta Optional scalar/vector mutation-rate parameter in MaCS units.
#' @param seed Optional scalar/vector integer seeds for mutation sampling.
#'
#' @return List with mutated table collections and metadata.
#' @keywords internal
#' @noRd
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
