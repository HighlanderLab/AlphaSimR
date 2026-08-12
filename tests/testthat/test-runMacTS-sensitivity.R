context("runMacTS wrapper and parameter sensitivity")

skip_if_not_installed("RcppTskit")

tc_summary <- function(tc_xptr) {
  AlphaSimR::rtsk_table_collection_summary2(tc_xptr)
}

node_times <- function(tc_xptr) {
  tc <- RcppTskit::TableCollection$new(xptr = tc_xptr)
  n_nodes <- as.integer(tc$num_nodes())
  vapply(seq_len(n_nodes), function(i) {
    as.numeric(tc$node_table_get_row(i - 1L)$time)
  }, numeric(1))
}

sample_node_count <- function(tc_xptr) {
  tc <- RcppTskit::TableCollection$new(xptr = tc_xptr)
  n_nodes <- as.integer(tc$num_nodes())
  flags <- vapply(seq_len(n_nodes), function(i) {
    as.integer(tc$node_table_get_row(i - 1L)$flags)
  }, integer(1))
  sum(bitwAnd(flags, 1L) != 0L)
}

site_position_range <- function(tc_xptr) {
  tc <- RcppTskit::TableCollection$new(xptr = tc_xptr)
  n_sites <- as.integer(tc$num_sites())
  if (n_sites == 0L) {
    return(c(min = NA_real_, max = NA_real_))
  }
  pos <- vapply(seq_len(n_sites), function(i) {
    as.numeric(tc$site_table_get_row(i - 1L)$position)
  }, numeric(1))
  c(min = min(pos), max = max(pos))
}

parse_seq_len <- function(args) {
  tokens <- strsplit(as.character(args), "[,[:space:]]+", perl = TRUE)[[1L]]
  tokens <- tokens[nzchar(tokens)]
  as.numeric(tokens[2L])
}

run_staged_from_wrapper <- function(out, inbred, ploidy, segSites,
                                    nThreads = 1L,
                                    usePhysicalPositions = FALSE,
                                    expandInbredTs = FALSE,
                                    siteSamplingSeed = 42L) {
  nChr <- length(out$tables)
  if (isTRUE(usePhysicalPositions)) {
    stop("run_staged_from_wrapper helper currently supports usePhysicalPositions = FALSE only")
  }
  
  anc <- AlphaSimR:::simAnc(
    args = out$args,
    nChr = nChr,
    inbred = inbred,
    ploidy = ploidy,
    nThreads = nThreads,
    seed = out$seed,
    usePhysicalPositions = usePhysicalPositions,
    Nref = NA_real_
  )
  
  dThetaPost <- as.numeric(anc$dTheta) / as.numeric(anc$timeScale)
  mutSeed <- if (!all(is.na(out$mutSeed))) {
    as.integer(out$mutSeed)
  } else {
    as.integer(out$seed) + 104729L
  }
  runOut <- AlphaSimR:::simMut(anc, dTheta = dThetaPost, seed = mutSeed)
  
  if (isTRUE(expandInbredTs) && isTRUE(inbred) && ploidy > 1L) {
    runOut <- AlphaSimR:::finalizeInbredTs(runOut, inbred = inbred, ploidy = ploidy)
  }
  
  breaks <- rep(list(c(0, 1)), nChr)
  rates <- rep(list(c(1)), nChr)
  pop <- AlphaSimR:::asMapPop(
    chr_info = list(tables = runOut$tables, breaks = breaks, rates = rates),
    ploidy = ploidy,
    inbred = inbred,
    segSites = segSites,
    site_sampling_seed = as.integer(siteSamplingSeed),
    nThreads = as.integer(nThreads),
    returnMeta = FALSE
  )
  
  list(tables = runOut$tables, pop = pop)
}

expect_wrapper_staged_equal <- function(out, staged) {
  expect_identical(out$pop@nLoci, staged$pop@nLoci)
  for (chr in seq_along(out$tables)) {
    expect_identical(out$pop@geno[[chr]], staged$pop@geno[[chr]])
    expect_true(isTRUE(all.equal(out$pop@genMap[[chr]], staged$pop@genMap[[chr]], tolerance = 0)))
    
    sw <- tc_summary(out$tables[[chr]])
    ss <- tc_summary(staged$tables[[chr]])
    expect_equal(sw$num_nodes, ss$num_nodes)
    expect_equal(sw$num_edges, ss$num_edges)
    expect_equal(sw$num_sites, ss$num_sites)
    expect_equal(sw$num_mutations, ss$num_mutations)
  }
}

test_that("runMacTS(postTs) wrapper matches staged workflow (outbred)", {
  out <- AlphaSimR:::runMacTS(
    nInd = 4,
    nChr = 2,
    segSites = 60,
    inbred = FALSE,
    ploidy = 2L,
    species = "GENERIC",
    mutationMode = "postTs",
    usePhysicalPositions = FALSE,
    nThreads = 1L,
    seed = as.integer(c(11, 22)),
    mutSeed = as.integer(c(111, 222)),
    siteSamplingSeed = 42L,
    returnTs = TRUE
  )
  
  staged <- run_staged_from_wrapper(
    out = out,
    inbred = FALSE,
    ploidy = 2L,
    segSites = 60,
    nThreads = 1L,
    usePhysicalPositions = FALSE,
    expandInbredTs = FALSE,
    siteSamplingSeed = 42L
  )
  
  expect_wrapper_staged_equal(out, staged)
})

test_that("runMacTS(postTs) wrapper matches staged workflow (inbred, ploidy > 1)", {
  out <- AlphaSimR:::runMacTS(
    nInd = 4,
    nChr = 2,
    segSites = 60,
    inbred = TRUE,
    ploidy = 2L,
    species = "GENERIC",
    mutationMode = "postTs",
    usePhysicalPositions = FALSE,
    expandInbredTs = TRUE,
    nThreads = 1L,
    seed = as.integer(c(11, 22)),
    mutSeed = as.integer(c(111, 222)),
    siteSamplingSeed = 42L,
    returnTs = TRUE
  )
  
  staged <- run_staged_from_wrapper(
    out = out,
    inbred = TRUE,
    ploidy = 2L,
    segSites = 60,
    nThreads = 1L,
    usePhysicalPositions = FALSE,
    expandInbredTs = TRUE,
    siteSamplingSeed = 42L
  )
  
  expect_wrapper_staged_equal(out, staged)
})

test_that("mutationMode='none' returns ancestry-only TS", {
  out <- AlphaSimR:::runMacTS(
    nInd = 4,
    nChr = 2,
    inbred = FALSE,
    ploidy = 2L,
    species = "GENERIC",
    mutationMode = "none",
    usePhysicalPositions = FALSE,
    nThreads = 1L,
    seed = as.integer(c(101, 202)),
    returnTs = TRUE
  )
  
  expect_null(out$pop)
  expect_identical(out$mutationMode, "none")
  expect_true(all(is.na(out$mutSeed)))
  
  for (chr in seq_along(out$tables)) {
    s <- tc_summary(out$tables[[chr]])
    expect_equal(s$num_sites, 0L)
    expect_equal(s$num_mutations, 0L)
    expect_gt(s$num_nodes, 0L)
    expect_gt(s$num_edges, 0L)
  }
})

test_that("usePhysicalPositions changes coordinate scale while keeping sampled output", {
  common <- list(
    nInd = 4,
    nChr = 1,
    segSites = 60,
    inbred = FALSE,
    ploidy = 2L,
    species = "GENERIC",
    mutationMode = "postTs",
    nThreads = 1L,
    seed = as.integer(123),
    mutSeed = as.integer(456),
    siteSamplingSeed = 42L,
    returnTs = TRUE
  )
  out_unit <- do.call(AlphaSimR:::runMacTS, c(common, list(usePhysicalPositions = FALSE)))
  out_bp <- do.call(AlphaSimR:::runMacTS, c(common, list(usePhysicalPositions = TRUE)))
  
  s_unit <- tc_summary(out_unit$tables[[1]])
  s_bp <- tc_summary(out_bp$tables[[1]])
  seq_len_bp <- parse_seq_len(out_bp$args)
  
  expect_equal(as.numeric(s_unit$sequence_length), 1.0, tolerance = 0)
  expect_equal(as.numeric(s_bp$sequence_length), seq_len_bp, tolerance = 0)
  expect_equal(s_unit$num_sites, s_bp$num_sites)
  expect_equal(s_unit$num_mutations, s_bp$num_mutations)
  
  rng_unit <- site_position_range(out_unit$tables[[1]])
  rng_bp <- site_position_range(out_bp$tables[[1]])
  expect_true(rng_unit["min"] >= 0 && rng_unit["max"] < 1)
  expect_true(rng_bp["min"] >= 0 && rng_bp["max"] < seq_len_bp)
  
  expect_identical(out_unit$pop@nLoci, out_bp$pop@nLoci)
  expect_identical(out_unit$pop@geno[[1]], out_bp$pop@geno[[1]])
  expect_true(isTRUE(all.equal(out_unit$pop@genMap[[1]], out_bp$pop@genMap[[1]], tolerance = 1e-12)))
})

test_that("runMacTSBridgeChrInfo reuses runMacTS map metadata", {
  out <- AlphaSimR:::runMacTS(
    nInd = 4,
    nChr = 2,
    segSites = c(40L, 40L),
    inbred = FALSE,
    ploidy = 2L,
    species = "GENERIC",
    mutationMode = "postTs",
    usePhysicalPositions = TRUE,
    nThreads = 1L,
    seed = as.integer(c(123, 456)),
    mutSeed = as.integer(c(789, 987)),
    siteSamplingSeed = 42L,
    returnTs = TRUE
  )
  out_dir <- tempfile("mac_bridge_chr_info_")
  chr_info <- AlphaSimR:::runMacTSBridgeChrInfo(out, out_dir = out_dir)

  pos_meta <- attr(out$pop, "tsForwardPosMeta", exact = TRUE)

  expect_length(chr_info, 2L)
  expect_true(all(file.exists(vapply(chr_info, `[[`, character(1), "ts_path"))))
  expect_equal(lapply(chr_info, `[[`, "breaks"), pos_meta$breaksList)
  expect_equal(lapply(chr_info, `[[`, "rates"), pos_meta$ratesList)
  expect_equal(vapply(chr_info, `[[`, integer(1), "segSites"), out$pop@nLoci)
})

test_that("runMacTS converts MaCS -R hotspot file into breaks/rates metadata", {
  hotspot_path <- tempfile("macs_hotspot_", fileext = ".txt")
  writeLines("0.25 0.5 3", hotspot_path)
  on.exit(unlink(hotspot_path, force = TRUE), add = TRUE)

  out <- AlphaSimR:::runMacTS(
    nInd = 4,
    nChr = 1,
    segSites = 30L,
    inbred = FALSE,
    ploidy = 2L,
    manualCommand = paste("1e5 -t 1e-3 -r 4e-6 -R", hotspot_path),
    manualGenLen = 1,
    mutationMode = "postTs",
    usePhysicalPositions = FALSE,
    nThreads = 1L,
    seed = as.integer(321),
    mutSeed = as.integer(654),
    siteSamplingSeed = 42L,
    returnTs = TRUE
  )

  pos_meta <- attr(out$pop, "tsForwardPosMeta", exact = TRUE)
  expect_equal(pos_meta$breaksList[[1L]], c(0, 0.25, 0.5, 1))
  expect_equal(pos_meta$ratesList[[1L]], c(1, 3, 1))
})

test_that("runMacs genMap conversion uses MaCS -R hotspot map", {
  hotspot_path <- tempfile("macs_hotspot_", fileext = ".txt")
  writeLines("0.25 0.5 3", hotspot_path)
  on.exit(unlink(hotspot_path, force = TRUE), add = TRUE)

  macs_pos <- list(c(0.10, 0.30, 0.60, 0.90))
  manual_command <- paste("1e5 -t 1e-3 -r 4e-6 -R", hotspot_path)
  out <- AlphaSimR:::.runMacsGenMapFromMacs(
    macsGenMap = macs_pos,
    genLen = 1,
    manualCommand = manual_command
  )

  expected <- AlphaSimR:::rateMap2cumMorgan(
    macs_pos[[1]],
    breaks = c(0, 0.25, 0.5, 1),
    rates = c(1, 3, 1)
  )
  expected <- expected - expected[[1L]]

  expect_equal(unname(out[[1]]), expected, tolerance = 1e-12)
  expect_identical(names(out[[1]]), paste(1, seq_along(macs_pos[[1]]), sep = "_"))
})

test_that("asMapPop requires breaks/rates for external tree input", {
  ts_path <- testthat::test_path("..", "..", "dev", "testData", "msprime_chr0.trees")
  skip_if(!file.exists(ts_path), "Missing test fixture dev/testData/msprime_chr0.trees")

  expect_error(
    AlphaSimR:::asMapPop(chr_info = list(list(ts_path = ts_path, segSites = 10L))),
    regexp = "breaks.*rates|rates.*breaks"
  )
})

test_that("Nref rescales TS times in runMacTS ancestry-only mode", {
  seed <- as.integer(42)
  nref <- 10000
  scale <- 4 * nref
  
  out_unit <- AlphaSimR:::runMacTS(
    nInd = 4,
    nChr = 1,
    inbred = FALSE,
    ploidy = 2L,
    species = "GENERIC",
    mutationMode = "none",
    usePhysicalPositions = FALSE,
    Nref = NA_real_,
    nThreads = 1L,
    seed = seed,
    returnTs = TRUE
  )
  out_gen <- AlphaSimR:::runMacTS(
    nInd = 4,
    nChr = 1,
    inbred = FALSE,
    ploidy = 2L,
    species = "GENERIC",
    mutationMode = "none",
    usePhysicalPositions = FALSE,
    Nref = nref,
    nThreads = 1L,
    seed = seed,
    returnTs = TRUE
  )
  
  expect_equal(out_unit$timeScale, 1)
  expect_equal(out_gen$timeScale, scale)
  
  s_unit <- tc_summary(out_unit$tables[[1]])
  s_gen <- tc_summary(out_gen$tables[[1]])
  expect_identical(s_unit$time_units, "unknown")
  expect_identical(s_gen$time_units, "generations")
  expect_equal(s_unit$num_nodes, s_gen$num_nodes)
  expect_equal(s_unit$num_edges, s_gen$num_edges)
  expect_equal(s_unit$num_trees, s_gen$num_trees)
  
  t_unit <- node_times(out_unit$tables[[1]])
  t_gen <- node_times(out_gen$tables[[1]])
  expect_equal(length(t_unit), length(t_gen))
  expect_true(isTRUE(all.equal(t_gen, t_unit * scale, tolerance = 1e-8)))
})

test_that("expandInbredTs toggles inbred leaf expansion in TS", {
  out_no_expand <- AlphaSimR:::runMacTS(
    nInd = 4,
    nChr = 1,
    inbred = TRUE,
    ploidy = 2L,
    species = "GENERIC",
    mutationMode = "none",
    usePhysicalPositions = FALSE,
    expandInbredTs = FALSE,
    nThreads = 1L,
    seed = as.integer(7),
    returnTs = TRUE
  )
  out_expand <- AlphaSimR:::runMacTS(
    nInd = 4,
    nChr = 1,
    inbred = TRUE,
    ploidy = 2L,
    species = "GENERIC",
    mutationMode = "none",
    usePhysicalPositions = FALSE,
    expandInbredTs = TRUE,
    nThreads = 1L,
    seed = as.integer(7),
    returnTs = TRUE
  )
  
  s0 <- tc_summary(out_no_expand$tables[[1]])
  s1 <- tc_summary(out_expand$tables[[1]])
  expect_equal(s0$num_sites, 0L)
  expect_equal(s1$num_sites, 0L)
  expect_equal(s0$num_mutations, 0L)
  expect_equal(s1$num_mutations, 0L)
  expect_gt(s1$num_nodes, s0$num_nodes)
  expect_gt(s1$num_edges, s0$num_edges)
  
  n_sample_0 <- sample_node_count(out_no_expand$tables[[1]])
  n_sample_1 <- sample_node_count(out_expand$tables[[1]])
  expect_equal(n_sample_0, 4L)
  expect_equal(n_sample_1, 8L)
})
