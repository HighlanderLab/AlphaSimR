#library(reticulate)
library(jsonlite)

#use_virtualenv("~/r-reticulate-env", required = TRUE)
#msprime <- import("msprime")
#tskit <- import("tskit")

reticulate::py_run_string("
import numpy as np

def sample_segregating_variants(ts, segSites, seed):
    rng = np.random.default_rng(int(seed))

    kept_pos = []
    kept_g = []
    k = 0

    for v in ts.variants():
        g = v.genotypes

        if len(np.unique(g)) != 2:
            continue

        # reservoir sampling
        k += 1
        if len(kept_g) < segSites:
            kept_g.append(g.copy())
            kept_pos.append(v.site.position)
        else:
            # Prob. entry: j/k
            j = rng.integers(0, k)
            if j < segSites:
                kept_g[j] = g.copy()
                kept_pos[j] = v.site.position

    H = np.stack(kept_g, axis=1).astype(np.int8)
    P = np.array(kept_pos, dtype=float)
    return H, P
")

reticulate::py_run_string("
import numpy as np

def segregating_variants(ts):

    kept_pos = []
    kept_g = []

    for v in ts.variants():
        g = v.genotypes

        if len(np.unique(g)) != 2:
            continue

        kept_g.append(g.copy())
        kept_pos.append(v.site.position)

    H = np.stack(kept_g, axis=1).astype(np.int8)
    P = np.array(kept_pos, dtype=float)
    return H, P
")


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


ts2chrDataPy <- function(ts_path, breaks, rates, segSites, site_sampling_seed) {
  ts = tskit$load(ts_path)

  pos <- ts$tables$sites$position

  if (!is.null(segSites)) {
    # stopifnot(length(pos) >= segSites)
    if (length(pos) < segSites) {
      stop("Insufficient sites (only ", length(pos), " sites in the tree sequence).")
    }
    message(segSites, " variants sampled ", "(Random seed: ", site_sampling_seed, ")")
    out <- py$sample_segregating_variants(ts, segSites, seed=site_sampling_seed)
    #if (length(out[[2]]) < segSites) {
    #  warning("Not enough sites kept after filtering non-biallelic sites.")
    #}
    if (length(out[[2]]) < segSites) {
      stop("Insufficient sites (only ", length(out[[2]]), " sites after filtering non-biallelic sites).")
    }
    message(segSites, " variants sampled ", "(Random seed: ", site_sampling_seed, ")")
  }
  else {
    out <- py$segregating_variants(ts)
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
    # haplotypes <- list(H)
    haplotypes = list(H),
    keptPosBp = pos
  )
}

asMapPopPy <- function(chr_info, ploidy = 2, inbred = FALSE, segSites = NULL, site_sampling_seed = 42) {
  ploidy <<- ploidy
  chr_data <- lapply(chr_info, function(info) {
    ts2chrDataPy(
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
