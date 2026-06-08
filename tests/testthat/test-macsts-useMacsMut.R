context("MaCSTS useMacsMut compatibility")

to_int01_matrix <- function(x) {
  matrix(as.integer(x), nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
}

ts_get <- function(ts, name) {
  value <- ts[[name]]
  if (is.function(value)) {
    value()
  } else {
    value
  }
}

ts_variants_iterator <- function(ts) {
  variants <- ts$variants
  if (is.function(variants)) {
    variants()
  } else {
    variants
  }
}

extract_macs_chr <- function(macs_out, chr = 1L, nThreads = 1L) {
  pos <- as.numeric(macs_out$genMap[[chr]])
  n_sites <- length(pos)
  if (n_sites == 0L) {
    stop("zero sites produced; choose higher mutation settings for this test")
  }
  hap_raw <- AlphaSimR:::getHaplo(
    geno = macs_out$geno[chr],
    lociPerChr = as.integer(n_sites),
    lociLoc = as.integer(seq_len(n_sites)),
    nThreads = as.integer(nThreads)
  )
  list(pos = pos, hap = to_int01_matrix(hap_raw))
}

extract_ts_chr <- function(tc_xptr) {
  tc <- RcppTskit::TableCollection$new(xptr = tc_xptr)
  ts <- tc$tree_sequence()
  n_samples <- as.integer(ts_get(ts, "num_samples"))
  it <- ts_variants_iterator(ts)
  pos <- numeric(0)
  cols <- list()
  repeat {
    v <- it$next_variant()
    if (is.null(v)) {
      break
    }
    pos <- c(pos, as.numeric(v$position))
    cols[[length(cols) + 1L]] <- as.integer(v$genotypes)
  }
  hap <- if (length(cols) == 0L) {
    matrix(integer(0), nrow = n_samples, ncol = 0L)
  } else {
    do.call(cbind, cols)
  }
  list(
    pos = pos,
    hap = hap,
    num_sites = as.integer(ts$num_sites()),
    num_mutations = as.integer(ts$num_mutations())
  )
}

site_hap_keys <- function(pos, hap) {
  if (length(pos) == 0L) {
    return(character(0))
  }
  p <- format(signif(pos, 15), scientific = FALSE, trim = TRUE)
  vapply(
    seq_len(ncol(hap)),
    function(j) paste0(p[j], "|", paste(hap[, j], collapse = "")),
    character(1)
  )
}

compare_chr <- function(macs_chr, ts_chr) {
  ord_m <- order(macs_chr$pos)
  ord_t <- order(ts_chr$pos)
  m_pos <- macs_chr$pos[ord_m]
  t_pos <- ts_chr$pos[ord_t]
  m_hap <- macs_chr$hap[, ord_m, drop = FALSE]
  t_hap <- ts_chr$hap[, ord_t, drop = FALSE]

  same_nsites <- ncol(m_hap) == ncol(t_hap)
  same_positions_strict <- isTRUE(all.equal(m_pos, t_pos, tolerance = 0))
  same_hap_strict <- identical(m_hap, t_hap)

  keys_m <- sort(site_hap_keys(macs_chr$pos, macs_chr$hap))
  keys_t <- sort(site_hap_keys(ts_chr$pos, ts_chr$hap))
  same_site_hap_multiset <- identical(keys_m, keys_t)

  list(
    same_nsites = same_nsites,
    same_positions_strict = same_positions_strict,
    same_hap_strict = same_hap_strict,
    same_site_hap_multiset = same_site_hap_multiset,
    ts_num_sites = ts_chr$num_sites,
    ts_num_mutations = ts_chr$num_mutations
  )
}

run_case <- function(args, nChr, inbred, ploidy, seed, nThreads = 1L) {
  seed_vec <- rep(as.integer(seed), as.integer(nChr))
  macs <- AlphaSimR:::MaCS(
    args = args,
    maxSites = rep(0L, as.integer(nChr)),
    inbred = inbred,
    ploidy = as.integer(ploidy),
    nThreads = as.integer(nThreads),
    seed = seed_vec
  )
  ts_out <- AlphaSimR:::MaCSTS(
    args = args,
    nChr = as.integer(nChr),
    inbred = inbred,
    ploidy = as.integer(ploidy),
    nThreads = as.integer(nThreads),
    seed = seed_vec,
    usePhysicalPositions = FALSE,
    useMacsMut = TRUE
  )

  out <- lapply(seq_len(as.integer(nChr)), function(chr) {
    m <- extract_macs_chr(macs, chr = chr, nThreads = nThreads)
    t <- extract_ts_chr(ts_out$tables[[chr]])
    compare_chr(m, t)
  })
  out
}

test_that("MaCSTS(useMacsMut=TRUE) matches MaCS across representative scenarios", {
  scenarios <- list(
    list(
      name = "base_outbred_ploidy2",
      args = "8 50000 -t 1e-3 -r 1e-4 -s ",
      inbred = FALSE,
      ploidy = 2L
    ),
    list(
      name = "base_outbred_ploidy1",
      args = "8 50000 -t 1e-3 -r 1e-4 -s ",
      inbred = FALSE,
      ploidy = 1L
    ),
    list(
      name = "base_inbred_ploidy2",
      args = "8 50000 -t 1e-3 -r 1e-4 -s ",
      inbred = TRUE,
      ploidy = 2L
    ),
    list(
      name = "demography_eN",
      args = "8 50000 -t 1e-3 -r 1e-4 -eN 0.2 2.0 -eN 0.9 0.5 -s ",
      inbred = FALSE,
      ploidy = 2L
    ),
    list(
      name = "multipop_with_migration_change",
      args = "8 50000 -t 1e-3 -r 1e-4 -I 2 4 4 1e-2 -eM 0.5 5e-3 -s ",
      inbred = FALSE,
      ploidy = 2L
    ),
    list(
      name = "multipop_en_plus_join",
      args = "8 50000 -t 1e-3 -r 1e-4 -I 2 4 4 1e-2 -en 0.2 2 0.5 -ej 1.0 2 1 -s ",
      inbred = FALSE,
      ploidy = 2L
    )
  )

  for (sc in scenarios) {
    res <- run_case(
      args = sc$args,
      nChr = 1L,
      inbred = sc$inbred,
      ploidy = sc$ploidy,
      seed = 12345L,
      nThreads = 1L
    )[[1]]

    expect_true(res$same_nsites, info = sc$name)
    expect_true(res$same_site_hap_multiset, info = sc$name)
    expect_true(res$same_positions_strict, info = sc$name)
    expect_true(res$same_hap_strict, info = sc$name)
    expect_equal(res$ts_num_sites, res$ts_num_mutations, info = sc$name)
  }
})

test_that("MaCSTS(useMacsMut=TRUE) is reproducible across chromosomes for fixed seeds", {
  args <- "8 50000 -t 1e-3 -r 1e-4 -s "
  nChr <- 3L
  seed_vec <- as.integer(c(101, 202, 303))

  out_a <- AlphaSimR:::MaCSTS(
    args = args,
    nChr = nChr,
    inbred = FALSE,
    ploidy = 2L,
    nThreads = 1L,
    seed = seed_vec,
    usePhysicalPositions = FALSE,
    useMacsMut = TRUE
  )
  out_b <- AlphaSimR:::MaCSTS(
    args = args,
    nChr = nChr,
    inbred = FALSE,
    ploidy = 2L,
    nThreads = 1L,
    seed = seed_vec,
    usePhysicalPositions = FALSE,
    useMacsMut = TRUE
  )
  out_c <- AlphaSimR:::MaCSTS(
    args = args,
    nChr = nChr,
    inbred = FALSE,
    ploidy = 2L,
    nThreads = 1L,
    seed = seed_vec + 1L,
    usePhysicalPositions = FALSE,
    useMacsMut = TRUE
  )

  keys_from <- function(out) {
    lapply(seq_len(nChr), function(chr) {
      x <- extract_ts_chr(out$tables[[chr]])
      sort(site_hap_keys(x$pos, x$hap))
    })
  }

  keys_a <- keys_from(out_a)
  keys_b <- keys_from(out_b)
  keys_c <- keys_from(out_c)

  expect_identical(keys_a, keys_b)
  expect_true(any(!vapply(seq_len(nChr), function(i) identical(keys_a[[i]], keys_c[[i]]), logical(1))))
})

test_that("usePhysicalPositions changes coordinate scale only", {
  args <- "8 50000 -t 1e-3 -r 1e-4 -s "
  nChr <- 1L
  seed <- as.integer(777)
  seq_len_bp <- as.numeric(strsplit(args, "[,[:space:]]+", perl = TRUE)[[1]][2])

  out_unit <- AlphaSimR:::MaCSTS(
    args = args,
    nChr = nChr,
    inbred = FALSE,
    ploidy = 2L,
    nThreads = 1L,
    seed = seed,
    usePhysicalPositions = FALSE,
    useMacsMut = TRUE
  )
  out_bp <- AlphaSimR:::MaCSTS(
    args = args,
    nChr = nChr,
    inbred = FALSE,
    ploidy = 2L,
    nThreads = 1L,
    seed = seed,
    usePhysicalPositions = TRUE,
    useMacsMut = TRUE
  )

  chr_unit <- extract_ts_chr(out_unit$tables[[1]])
  chr_bp <- extract_ts_chr(out_bp$tables[[1]])

  expect_true(all(chr_unit$pos >= 0 & chr_unit$pos <= 1))
  expect_true(all(chr_bp$pos >= 0 & chr_bp$pos <= seq_len_bp))

  # Haplotypes should match exactly; positions should match after rescaling.
  expect_identical(sort(site_hap_keys(chr_unit$pos, chr_unit$hap)),
                   sort(site_hap_keys(chr_bp$pos / seq_len_bp, chr_bp$hap)))
  expect_true(isTRUE(all.equal(sort(chr_unit$pos), sort(chr_bp$pos / seq_len_bp), tolerance = 1e-12)))
})

test_that("MaCSTS validates key inputs", {
  args <- "8 50000 -t 1e-3 -r 1e-4 -s "

  expect_error(
    AlphaSimR:::MaCSTS(
      args = args,
      nChr = 0L,
      inbred = FALSE,
      ploidy = 2L,
      nThreads = 1L,
      seed = as.integer(1),
      usePhysicalPositions = FALSE,
      useMacsMut = TRUE
    ),
    "nChr must be a positive integer"
  )

  expect_error(
    AlphaSimR:::MaCSTS(
      args = args,
      nChr = 1L,
      inbred = FALSE,
      ploidy = 0L,
      nThreads = 1L,
      seed = as.integer(1),
      usePhysicalPositions = FALSE,
      useMacsMut = TRUE
    ),
    "ploidy must be a positive integer"
  )

  expect_error(
    AlphaSimR:::MaCSTS(
      args = args,
      nChr = 2L,
      inbred = FALSE,
      ploidy = 2L,
      nThreads = 1L,
      seed = as.integer(1),
      usePhysicalPositions = FALSE,
      useMacsMut = TRUE
    ),
    "seed length must match number of chromosomes"
  )

  expect_error(
    AlphaSimR:::MaCSTS(
      args = args,
      nChr = 1L,
      inbred = FALSE,
      ploidy = 2L,
      nThreads = 1L,
      seed = as.integer(1),
      usePhysicalPositions = FALSE,
      useMacsMut = TRUE,
      Nref = 0
    ),
    "Nref must be positive when provided"
  )
})
