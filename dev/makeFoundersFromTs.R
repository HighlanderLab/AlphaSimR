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


segregating_variants <- function(ts) {
  # 1. Get dimensions for pre-allocation
  max_sites <- as.integer(ts$num_sites())
  num_samples <- as.integer(ts$num_samples())

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

    # Filter for biallelic sites (exactly 2 unique alleles)
    if (length(unique(g)) == 2) {
      count <- count + 1
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

  return(list(H = H, P = P))
}

segregating_variants_debug <- function(ts) {
  # 1. Get dimensions for pre-allocation
  max_sites <- ts$num_sites()
  num_samples <- ts$num_samples()

  # DEBUG: Print initial metadata
  message(paste("Expected max sites:", max_sites))
  message(paste("Expected num samples (from ts):", num_samples))

  # 2. Pre-allocate H matrix
  H_full <- matrix(NA_integer_, nrow = num_samples, ncol = max_sites)
  P_full <- numeric(max_sites)

  it <- ts$variants()
  count <- 0

  # 3. Iterate through variants
  repeat {
    v <- it$next_variant()
    if (is.null(v)) break

    g <- v$genotypes

    # DEBUG: Check dimensions on the first iteration
    if (count == 0) {
      message(paste("Actual length of genotype vector (g):", length(g)))
      message(paste("Matrix H_full has", nrow(H_full), "rows"))

      if (length(g) != nrow(H_full)) {
        stop("DIMENSION MISMATCH: The genotype vector length does not match matrix rows!")
      }
    }

    # Filter for biallelic sites
    if (length(unique(g)) == 2) {
      count <- count + 1

      # DEBUG: Check for column overflow
      if (count > max_sites) {
        stop(paste("INDEX OVERFLOW: count (", count, ") exceeded max_sites (", max_sites, ")"))
      }

      # Fill the matrix column directly
      H_full[, count] <- g
      P_full[count] <- v$position
    }
  }

  # 4. Trim the results
  if (count > 0) {
    H <- H_full[, 1:count, drop = FALSE]
    P <- P_full[1:count]
  } else {
    H <- matrix(nrow = num_samples, ncol = 0)
    P <- numeric(0)
  }

  message(paste("Success! Final count of biallelic variants:", count))
  return(list(H = H, P = P))
}

# rec map used in msprime:
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


ts2chrData <- function(ts_path, breaks, rates, segSites, site_sampling_seed) {
  ts = ts_load(ts_path)
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
    message(segSites, " variants sampled ", "(Random seed: ", site_sampling_seed, ")")
  }
  else {
    out <- segregating_variants(ts)
  }

  H <- out[[1]]
  pos <- out[[2]]

  ordPos <- order(pos)

  pos <- pos[ordPos]

  mpos <- rateMap2cumMorgan(pos, breaks, rates)

  # relative position, so the 1st element is 0
  mpos <- mpos - min(mpos)

  ordMap <- order(mpos)
  mpos <- mpos[ordMap]
  pos <- pos[ordMap]
  H   <- H[, ordMap, drop = FALSE]


  list(
    genMap = list(mpos),
    haplotypes = list(H),
    keptPosBp = pos
  )
}

asMapPop <- function(chr_info, ploidy = 2, inbred = FALSE, segSites = NULL, site_sampling_seed = 42) {
  ploidy <<- ploidy
  chr_data <- lapply(chr_info, function(info) {
    ts2chrData(
      ts_path = info$ts_path,
      breaks  = info$breaks,
      rates   = info$rates,
      segSites = info$segSites,
      site_sampling_seed = site_sampling_seed
    )
  })

  # save pos in bp for tskit tables
  chrKeptPosBpList <<- lapply(chr_data, `[[`, "keptPosBp")

  genMap <- do.call(c, lapply(chr_data, `[[`, "genMap"))
  haplotypes <- do.call(c, lapply(chr_data, `[[`, "haplotypes"))

  newMapPop(genMap = genMap, haplotypes = haplotypes, inbred = inbred, ploidy = ploidy)
}
