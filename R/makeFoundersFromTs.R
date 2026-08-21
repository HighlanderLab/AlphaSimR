#' Sample Biallelic Variants from a Tree Sequence
#'
#' @param ts A `RcppTskit::TreeSequence` object.
#' @param segSites Integer number of biallelic variants to sample.
#' @param seed Integer random seed for reservoir sampling.
#'
#' @return A list with `H` (haplotypes matrix; rows are samples) and
#'   `P` (numeric vector of site positions).
#' @keywords internal
#' @noRd
sample_segregating_variants <- function(ts, segSites, seed) {

  # Sample segregating variants from the tree sequence.
  #
  # Parameters
  # ==========
  # ts: tskit.TreeSequence
  #     The tree sequence to sample from.
  # segSites: int
  #     The number of segregating sites to sample.
  # seed: int
  #     The random seed to use for sampling.
  #
  # Returns
  # =======
  # list of int
  #     The positions of the sampled segregating sites.
  # Set the random seed for reproducibility.
  set.seed(seed)
  num_samples <- as.integer(ts$num_samples())

  # 2. Pre-allocate H matrix and P vector based on required sample size (segSites)
  # We only need space for 'segSites' number of variants
  H <- matrix(NA_integer_, nrow = num_samples, ncol = segSites)
  P <- numeric(segSites)

  it <- ts$variants()

  # k tracks how many biallelic variants we have encountered so far
  k <- 0
  # current_size tracks how many variants are currently in our reservoir
  current_size <- 0
  # 3. Iterate through variants
  repeat {
    v <- it$next_variant()
    if (is.null(v)) break

    g <- v$genotypes

    # Filter for biallelic sites
    if (length(unique(g)) == 2) {
      k <- k + 1

      if (current_size < segSites) {
        # Case A: Reservoir is not full yet
        current_size <- current_size + 1
        H[, current_size] <- g
        P[current_size] <- v$position
      } else {
        # Case B: Reservoir is full, use Prob. entry: j/k
        # sample.int(k, 1) returns a value from 1 to k
        j <- sample.int(k, 1)

        if (j <= segSites) {
          # Replace the existing variant at index j
          H[, j] <- g
          P[j] <- v$position
        }
      }
    }
  }

  # 4. Final check: if we found fewer biallelic sites than segSites, trim the output
  if (k < segSites) {
    if (k > 0) {
      H <- H[, 1:k, drop = FALSE]
      P <- P[1:k]
    } else {
      H <- matrix(nrow = num_samples, ncol = 0)
      P <- numeric(0)
    }
  }

  return(list(H = H, P = P))
}


#' Extract All Biallelic Variants from a Tree Sequence
#'
#' @param ts A `RcppTskit::TreeSequence` object.
#' @param debug Logical; if `TRUE`, print diagnostics while scanning variants.
#'
#' @return A list with `H` (haplotypes matrix; rows are samples) and
#'   `P` (numeric vector of site positions).
#' @keywords internal
#' @noRd
segregating_variants <- function(ts, debug = FALSE) {
  # 1. Get dimensions for pre-allocation
  max_sites <- as.integer(ts$num_sites())
  num_samples <- as.integer(ts$num_samples())
  if (debug) {
    message("Expected max sites: ", max_sites)
    message("Expected num samples (from ts): ", num_samples)
  }

  # 2. Pre-allocate H matrix (Rows: samples, Cols: sites)
  # Using integer matrix to save memory (similar to np.int8)
  H_full <- matrix(NA_integer_, nrow = num_samples, ncol = max_sites)
  # Pre-allocate P vector for positions
  P_full <- numeric(max_sites)

  it <- ts$variants()
  count <- 0

  # 3. Iterate through variants
  repeat {
    v <- it$next_variant()
    if (is.null(v)) break

    g <- v$genotypes
    if (debug && count == 0L) {
      message("Actual length of genotype vector (g): ", length(g))
      message("Matrix H_full has ", nrow(H_full), " rows")
      if (length(g) != nrow(H_full)) {
        stop("DIMENSION MISMATCH: The genotype vector length does not match matrix rows!")
      }
    }

    # Filter for biallelic sites (exactly 2 unique alleles)
    if (length(unique(g)) == 2) {
      count <- count + 1
      if (debug && count > max_sites) {
        stop("INDEX OVERFLOW: count (", count, ") exceeded max_sites (", max_sites, ")")
      }
      # Fill the matrix column directly
      H_full[, count] <- g
      P_full[count] <- v$position
    }
  }

  # 4. Trim the results to the actual number of kept variants
  if (count > 0) {
    H <- H_full[, 1:count, drop = FALSE]
    P <- P_full[1:count]
  } else {
    H <- matrix(nrow = num_samples, ncol = 0)
    P <- numeric(0)
  }
  if (debug) {
    message("Success! Final count of biallelic variants: ", count)
  }

  return(list(H = H, P = P))
}

#' Convert Physical Positions to Cumulative Morgan Positions
#'
#' @param x Numeric vector of physical positions.
#' @param breaks Numeric vector of recombination map breakpoints.
#' @param rates Numeric vector of per-bp recombination rates for each interval.
#'
#' @return Numeric vector of cumulative Morgan positions.
#' @keywords internal
#' @noRd
rateMap2cumMorgan <- function(x, breaks, rates) {
  stopifnot(length(breaks) == length(rates) + 1)

  o <- order(breaks)
  breaks <- breaks[o]

  # M_i = m(breaks[i])
  seg_len <- diff(breaks)
  M_start <- c(0, cumsum(rates * seg_len))  # length = length(breaks)

  i <- findInterval(x, breaks, rightmost.closed = FALSE)
  i <- pmin(pmax(i, 1), length(rates))

  m <- M_start[i] + rates[i] * (x - breaks[i])
  return(m)
}


#' Load Tree Sequence from Supported Sources
#'
#' @param ts_path Character path to a `.trees` file.
#' @param ts Optional `TreeSequence` object or external pointer.
#' @param ts_xptr Optional external pointer to `tsk_treeseq_t`.
#' @param tc_xptr Optional external pointer to `tsk_table_collection_t`.
#' @param table_xptr Alias for `tc_xptr`.
#'
#' @return A `RcppTskit::TreeSequence` object.
#' @keywords internal
#' @noRd
asMapPop_load_ts <- function(ts_path = NULL, ts = NULL, ts_xptr = NULL,
                             tc_xptr = NULL, table_xptr = NULL) {
  if (!is.null(ts)) {
    if (inherits(ts, "externalptr")) {
      try_ts <- try(RcppTskit::TreeSequence$new(xptr = ts), silent = TRUE)
      if (!inherits(try_ts, "try-error")) {
        return(try_ts)
      }
      tc <- RcppTskit::TableCollection$new(xptr = ts)
      return(tc$tree_sequence())
    }
    return(ts)
  }
  if (!is.null(ts_xptr)) {
    return(RcppTskit::TreeSequence$new(xptr = ts_xptr))
  }
  tc_ptr <- if (!is.null(tc_xptr)) tc_xptr else table_xptr
  if (!is.null(tc_ptr)) {
    tc <- RcppTskit::TableCollection$new(xptr = tc_ptr)
    return(tc$tree_sequence())
  }
  if (!is.null(ts_path)) {
    return(RcppTskit::ts_load(ts_path))
  }
  stop("No tree-sequence source provided. Provide one of ts_path, ts, ts_xptr, tc_xptr, or table_xptr.")
}

#' Convert One Chromosome Tree Sequence to Map/Haplotypes
#'
#' @param ts_path Character path to `.trees` file (optional).
#' @param breaks Numeric vector of recombination map breakpoints.
#' @param rates Numeric vector of per-bp recombination rates.
#' @param segSites Optional integer number of biallelic variants to keep.
#' @param site_sampling_seed Integer seed for site sampling.
#' @param ts Optional in-memory tree sequence object.
#' @param ts_xptr Optional tree-sequence external pointer.
#' @param tc_xptr Optional table-collection external pointer.
#' @param table_xptr Alias for `tc_xptr`.
#'
#' @return A list with `genMap`, `haplotypes`, and `keptPosBp`.
#' @keywords internal
#' @noRd
ts2chrData <- function(ts_path = NULL, breaks, rates, segSites, site_sampling_seed,
                       ts = NULL, ts_xptr = NULL, tc_xptr = NULL, table_xptr = NULL) {
  ts <- asMapPop_load_ts(
    ts_path = ts_path,
    ts = ts,
    ts_xptr = ts_xptr,
    tc_xptr = tc_xptr,
    table_xptr = table_xptr
  )
  seqLen <- as.numeric(ts$sequence_length())
  num_pos <- ts$num_sites()

  if (!is.null(segSites)) {

    if (num_pos < segSites) {
      stop("Insufficient sites (only ", num_pos, " sites in the tree sequence).")
    }
    message(segSites, " variants sampled ", "(Random seed: ", site_sampling_seed, ")")
    out <- sample_segregating_variants(ts, segSites, site_sampling_seed)

    if (length(out[[2]]) < segSites) {
      stop("Insufficient sites (only ", length(out[[2]]), " sites after filtering non-biallelic sites).")
    }
  }
  else {
    out <- segregating_variants(ts)
  }

  H <- out[[1]]
  pos <- out[[2]]

  ord <- order(pos)
  pos <- pos[ord]
  H   <- H[, ord, drop = FALSE]
  mpos <- rateMap2cumMorgan(pos, breaks, rates)

  # relative position, so the 1st element is 0
  mpos <- mpos - min(mpos)

  list(
    genMap = list(mpos),
    haplotypes = list(H),
    keptPosBp = pos,
    seqLen = seqLen
  )
}

#' Expand Recombination Component to Per-Chromosome List
#'
#' @param x Scalar/list recombination component (`breaks` or `rates`).
#' @param nChr Integer number of chromosomes.
#' @param name Character label for error messages.
#'
#' @return A list of length `nChr`.
#' @keywords internal
#' @noRd
.asMapPop_expand_rec_component <- function(x, nChr, name) {
  if (is.null(x)) {
    stop("Missing `", name, "` for asMapPop input.")
  }
  values <- if (is.list(x)) x else list(x)
  if (length(values) == 1L) {
    return(rep(values, nChr))
  }
  if (length(values) != nChr) {
    stop("`", name, "` must have length 1 or nChr.")
  }
  values
}

#' Extract First Non-NULL Alias from a Named List
#'
#' @param x Named list.
#' @param keys Character vector of alias keys to try in order.
#'
#' @return First non-`NULL` component found, or `NULL`.
#' @keywords internal
#' @noRd
.asMapPop_extract_component <- function(x, keys) {
  for (k in keys) {
    if (!is.null(x[[k]])) {
      return(x[[k]])
    }
  }
  NULL
}

#' Safely Get Exact Named Component
#'
#' @param x Named list.
#' @param key Character scalar key.
#'
#' @return Value for `key` or `NULL` if absent.
#' @keywords internal
#' @noRd
.asMapPop_get <- function(x, key) {
  if (!is.list(x) || is.null(names(x))) {
    return(NULL)
  }
  if (!(key %in% names(x))) {
    return(NULL)
  }
  x[[key]]
}

#' Expand segSites to Per-Chromosome Specification
#'
#' @param segSites Optional scalar/vector/list of site counts.
#' @param nChr Integer number of chromosomes.
#' @param defaults Optional fallback segSites specification.
#'
#' @return List of length `nChr` with integer values or `NULL`.
#' @keywords internal
#' @noRd
.asMapPop_expand_seg_sites <- function(segSites, nChr, defaults = NULL) {
  if (is.null(segSites)) {
    values <- defaults
  } else {
    values <- segSites
  }
  if (is.null(values)) {
    return(rep(list(NULL), nChr))
  }
  if (!is.list(values)) {
    values <- as.list(as.integer(values))
  }
  if (length(values) == 1L) {
    values <- rep(values, nChr)
  } else if (length(values) != nChr) {
    stop("`segSites` must have length 1 or nChr.")
  }
  lapply(values, function(x) {
    if (is.null(x) || length(x) == 0) {
      return(NULL)
    }
    x <- as.integer(x[1])
    if (is.na(x) || x <= 0L) {
      return(NULL)
    }
    x
  })
}

#' Check Whether an Entry Looks Like Chromosome TS Specification
#'
#' @param x List candidate chromosome specification.
#'
#' @return Logical scalar.
#' @keywords internal
#' @noRd
.asMapPop_is_chr_info <- function(x) {
  is.list(x) && (
    !is.null(x$ts_path) ||
      !is.null(x$ts) ||
      !is.null(x$ts_xptr) ||
      !is.null(x$tc_xptr) ||
      !is.null(x$table_xptr)
  )
}

#' Normalize asMapPop Inputs to Per-Chromosome Specs
#'
#' @param chr_info Either explicit per-chromosome list or bundle style input
#'   containing `tables`/`ts` and map metadata.
#' @param segSites Optional override for per-chromosome site counts.
#'
#' @return A normalized list of per-chromosome specs.
#' @keywords internal
#' @noRd
.asMapPop_prepare_specs <- function(chr_info, segSites = NULL) {
  if (!is.list(chr_info)) {
    stop("`chr_info` must be a list.")
  }

  root_info <- chr_info
  if (!is.null(chr_info$chr_info)) {
    chr_info <- chr_info$chr_info
  }

  is_explicit_chr_info <- length(chr_info) > 0 &&
    all(vapply(chr_info, .asMapPop_is_chr_info, logical(1)))

  if (is_explicit_chr_info) {
    nChr <- length(chr_info)
    default_seg <- lapply(chr_info, function(x) x$segSites)
    seg_by_chr <- .asMapPop_expand_seg_sites(segSites, nChr, defaults = default_seg)
    root_breaks <- .asMapPop_extract_component(root_info, c("breaks", "rec_breaks", "recBreaks"))
    root_rates <- .asMapPop_extract_component(root_info, c("rates", "rec_rates", "recRates"))
    if (!is.null(root_breaks)) {
      root_breaks <- .asMapPop_expand_rec_component(root_breaks, nChr, "breaks")
    }
    if (!is.null(root_rates)) {
      root_rates <- .asMapPop_expand_rec_component(root_rates, nChr, "rates")
    }
    out <- vector("list", nChr)
    for (i in seq_len(nChr)) {
      info <- chr_info[[i]]
      if (is.null(info$breaks) && !is.null(root_breaks)) {
        info$breaks <- root_breaks[[i]]
      }
      if (is.null(info$rates) && !is.null(root_rates)) {
        info$rates <- root_rates[[i]]
      }
      if (is.null(info$breaks) || is.null(info$rates)) {
        stop("Each chromosome entry must include `breaks` and `rates`.")
      }
      info$segSites <- seg_by_chr[[i]]
      out[[i]] <- info
    }
    return(out)
  }

  tables <- .asMapPop_extract_component(chr_info, c("tables", "table_collections"))
  ts_list <- .asMapPop_extract_component(chr_info, c("ts", "tree_sequences"))
  if (is.null(tables) && is.null(ts_list)) {
    stop("Unsupported `chr_info` format. Provide a list of chromosome specs or a bundle with `tables`/`ts`.")
  }
  if (!is.null(tables) && !is.list(tables)) {
    stop("`tables` must be a list.")
  }
  if (!is.null(ts_list) && !is.list(ts_list)) {
    stop("`ts` must be a list.")
  }

  if (!is.null(tables)) {
    nChr <- length(tables)
  } else {
    nChr <- length(ts_list)
  }
  if (!is.null(tables) && !is.null(ts_list) && length(ts_list) != nChr) {
    stop("`tables` and `ts` must have the same length when both are supplied.")
  }

  breaks <- .asMapPop_expand_rec_component(
    .asMapPop_extract_component(chr_info, c("breaks", "rec_breaks", "recBreaks")),
    nChr, "breaks"
  )
  rates <- .asMapPop_expand_rec_component(
    .asMapPop_extract_component(chr_info, c("rates", "rec_rates", "recRates")),
    nChr, "rates"
  )
  default_seg <- .asMapPop_extract_component(chr_info, c("segSites", "seg_sites"))
  seg_by_chr <- .asMapPop_expand_seg_sites(segSites, nChr, defaults = default_seg)

  out <- vector("list", nChr)
  for (i in seq_len(nChr)) {
    out[[i]] <- list(
      ts = if (!is.null(ts_list)) ts_list[[i]] else NULL,
      tc_xptr = if (!is.null(tables)) tables[[i]] else NULL,
      breaks = breaks[[i]],
      rates = rates[[i]],
      segSites = seg_by_chr[[i]]
    )
  }
  out
}

#' Resolve and Validate Thread Count for asMapPop
#'
#' @param nThreads Optional requested thread count.
#'
#' @return Integer thread count >= 1.
#' @keywords internal
#' @noRd
.asMapPop_get_num_threads <- function(nThreads) {
  if (is.null(nThreads)) {
    if (exists("getNumThreads", mode = "function")) {
      nThreads <- getNumThreads()
    } else {
      nThreads <- 1L
    }
  }
  nThreads <- as.integer(nThreads)
  if (length(nThreads) != 1L || is.na(nThreads) || nThreads < 1L) {
    stop("`nThreads` must be a single positive integer.")
  }
  nThreads
}

#' Apply Chromosome Worker with Optional Parallelism
#'
#' @param chr_specs Normalized per-chromosome specs.
#' @param worker Function applied to each chromosome spec.
#' @param nThreads Integer requested threads.
#'
#' @return List of worker outputs.
#' @keywords internal
#' @noRd
.asMapPop_apply <- function(chr_specs, worker, nThreads) {
  if (length(chr_specs) <= 1L || nThreads <= 1L) {
    return(lapply(chr_specs, worker))
  }
  has_non_file_source <- any(vapply(chr_specs, function(x) {
    is.null(x$ts_path)
  }, logical(1)))
  if (has_non_file_source) {
    warning(
      "asMapPop: parallel conversion currently uses file-backed TS only. ",
      "Falling back to serial for in-memory TS/tables.",
      call. = FALSE
    )
    return(lapply(chr_specs, worker))
  }
  if (.Platform$OS.type == "unix") {
    return(parallel::mclapply(chr_specs, worker, mc.cores = nThreads))
  }
  warning(
    "asMapPop: parallel conversion is only enabled on unix via mclapply. ",
    "Falling back to serial on this platform.",
    call. = FALSE
  )
  lapply(chr_specs, worker)
}

.asMapPop_ts_source_from_spec <- function(info) {
  out <- list(
    ts_path = .asMapPop_get(info, "ts_path"),
    ts = .asMapPop_get(info, "ts"),
    ts_xptr = .asMapPop_get(info, "ts_xptr"),
    tc_xptr = .asMapPop_get(info, "tc_xptr"),
    table_xptr = .asMapPop_get(info, "table_xptr")
  )
  has_source <- vapply(out, function(x) !is.null(x), logical(1))
  out[has_source]
}

#' Build a MapPop from Tree Sequence Data
#'
#' @description
#' Converts one or more tree sequences to an AlphaSimR
#' \code{\link{MapPop-class}} by extracting biallelic segregating variants and
#' mapping tree-sequence coordinates through a recombination map. The resulting
#' population keeps the metadata needed for forward tree-sequence recording.
#'
#' @param chr_info Input tree-sequence data. Supports either:
#'   1) explicit per-chromosome list entries with \code{ts_path},
#'   \code{ts}, or \code{tc_xptr} plus \code{breaks} and \code{rates}; or
#'   2) bundle style list containing \code{tables} or \code{ts} plus map
#'   metadata.
#' @param ploidy Integer ploidy used to construct the resulting
#' \code{\link{MapPop-class}}.
#' @param inbred Logical; whether resulting individuals are inbred.
#' @param segSites Optional site-count override (scalar or per chromosome).
#' @param site_sampling_seed Integer seed used when downsampling segregating sites.
#' @param nThreads Optional chromosome-level worker count.
#' @param returnMeta Logical; if \code{TRUE}, return list with \code{pop},
#' \code{keptPosBp}, and \code{chrData}; otherwise return
#' \code{\link{MapPop-class}} only.
#'
#' @details
#' Each chromosome must provide tree-sequence input and a recombination map.
#' The map is supplied as \code{breaks} and \code{rates}, where
#' \code{breaks} are tree-sequence coordinate breakpoints and \code{rates} are
#' recombination rates for the corresponding intervals. The physical or
#' tree-sequence coordinates of retained variants are stored on the returned
#' population as \code{tsForwardPosMeta}.
#'
#' If \code{segSites} is supplied, biallelic segregating variants are sampled
#' with reservoir sampling using \code{site_sampling_seed}. Non-biallelic and
#' non-segregating variants are ignored.
#'
#' @return A \code{\link{MapPop-class}} object, or metadata list if
#' \code{returnMeta = TRUE}.
#'
#' @examples
#' \dontrun{
#' chr_info = list(list(
#'   ts_path="dev/testData/msprime_chr0.trees",
#'   breaks=c(0, 1),
#'   rates=c(1)
#' ))
#' founderPop = asMapPop(chr_info=chr_info, ploidy=2L)
#' }
#'
#' @export
asMapPop <- function(chr_info, ploidy = 2L, inbred = FALSE, segSites = NULL,
                     site_sampling_seed = 42L, nThreads = NULL,
                     returnMeta = FALSE) {
  ploidy <- as.integer(ploidy)
  nThreads <- .asMapPop_get_num_threads(nThreads)
  chr_specs <- .asMapPop_prepare_specs(chr_info, segSites = segSites)

  worker <- function(info) {
    ts2chrData(
      ts_path = .asMapPop_get(info, "ts_path"),
      ts = .asMapPop_get(info, "ts"),
      ts_xptr = .asMapPop_get(info, "ts_xptr"),
      tc_xptr = .asMapPop_get(info, "tc_xptr"),
      table_xptr = .asMapPop_get(info, "table_xptr"),
      breaks = .asMapPop_get(info, "breaks"),
      rates = .asMapPop_get(info, "rates"),
      segSites = .asMapPop_get(info, "segSites"),
      site_sampling_seed = site_sampling_seed
    )
  }

  chr_data <- .asMapPop_apply(chr_specs, worker, nThreads = nThreads)

  # save pos in bp for tskit tables
  chrKeptPosBp <- lapply(chr_data, `[[`, "keptPosBp")
  chrSeqLenBp <- lapply(chr_data, `[[`, "seqLen")
  chrKeptPosBpList <<- chrKeptPosBp
  chrSeqLenBpList <<- chrSeqLenBp
  # generic aliases: positions are in tree-sequence coordinate space
  # (bp when physical positions are used, otherwise normalized coordinates)
  chrKeptPosTsList <<- chrKeptPosBp
  chrSeqLenTsList <<- chrSeqLenBp

  genMap <- do.call(c, lapply(chr_data, `[[`, "genMap"))
  haplotypes <- do.call(c, lapply(chr_data, `[[`, "haplotypes"))

  pop <- newMapPop(genMap = genMap, haplotypes = haplotypes, inbred = inbred, ploidy = ploidy)
  attr(pop, "tsForwardSource") <- lapply(chr_specs, .asMapPop_ts_source_from_spec)
  attr(pop, "tsForwardPosMeta") <- list(
    posList = chrKeptPosBp,
    seqLenList = chrSeqLenBp,
    breaksList = lapply(chr_specs, function(x) as.numeric(.asMapPop_get(x, "breaks"))),
    ratesList = lapply(chr_specs, function(x) as.numeric(.asMapPop_get(x, "rates")))
  )
  if (isTRUE(returnMeta)) {
    return(list(
      pop = pop,
      keptPosBp = chrKeptPosBp,
      chrData = chr_data
    ))
  }
  pop
}
