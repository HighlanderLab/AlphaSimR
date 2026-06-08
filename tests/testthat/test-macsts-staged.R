context("MaCSTS staged simAnc/simMut checks")

skip_if_not_installed("RcppTskit")

ts_get <- function(ts, name) {
  value <- ts[[name]]
  if (is.function(value)) {
    value()
  } else {
    value
  }
}

ts_variants_iterator <- function(ts) {
  variants <- ts[["variants"]]
  if (is.function(variants)) {
    variants()
  } else {
    variants
  }
}

ts_next_variant <- function(it) {
  if (is.null(it)) {
    return(NULL)
  }
  if (is.function(it)) {
    return(it())
  }
  
  nxt <- it[["next_variant"]]
  if (!is.null(nxt)) {
    if (is.function(nxt)) {
      return(nxt())
    }
    return(nxt)
  }
  
  nxt <- it[["next"]]
  if (!is.null(nxt)) {
    if (is.function(nxt)) {
      return(nxt())
    }
    return(nxt)
  }
  
  return(NULL)
}

ts_variant_keys <- function(tc_xptr) {
  tc <- RcppTskit::TableCollection$new(xptr = tc_xptr)
  ts <- tc$tree_sequence()
  it <- ts_variants_iterator(ts)
  keys <- character(0)
  repeat {
    v <- ts_next_variant(it)
    if (is.null(v)) {
      break
    }
    pos <- format(signif(as.numeric(v$position), 15),
                  scientific = FALSE, trim = TRUE)
    keys <- c(keys, paste0(pos, "|", paste(as.integer(v$genotypes), collapse = "")))
  }
  sort(keys)
}

table_counts <- function(tc_xptr) {
  tc <- RcppTskit::TableCollection$new(xptr = tc_xptr)
  ts <- tc$tree_sequence()
  list(
    num_sites = as.integer(tc$num_sites()),
    num_mutations = as.integer(tc$num_mutations()),
    num_nodes = as.integer(tc$num_nodes()),
    num_edges = as.integer(tc$num_edges()),
    num_trees = as.integer(ts_get(ts, "num_trees"))
  )
}

node_times <- function(tc) {
  n_nodes <- as.integer(tc$num_nodes())
  vapply(seq_len(n_nodes), function(i) {
    as.numeric(tc$node_table_get_row(i - 1L)$time)
  }, numeric(1))
}

test_that("simAnc + simMut staged workflow is reproducible for fixed seeds", {
  args <- "8 5000 -t 1e-3 -r 1e-4 -s "
  nChr <- 2L
  seed <- as.integer(c(101, 202))
  mut_seed <- as.integer(c(555, 666))
  dTheta <- c(80, 80)
  
  anc_a <- AlphaSimR:::simAnc(
    args = args,
    nChr = nChr,
    inbred = FALSE,
    ploidy = 2L,
    nThreads = 1L,
    seed = seed,
    usePhysicalPositions = FALSE,
    Nref = NA_real_
  )
  anc_b <- AlphaSimR:::simAnc(
    args = args,
    nChr = nChr,
    inbred = FALSE,
    ploidy = 2L,
    nThreads = 1L,
    seed = seed,
    usePhysicalPositions = FALSE,
    Nref = NA_real_
  )
  
  expect_identical(anc_a$stage, "simAnc")
  expect_identical(anc_a$mutationMode, "none")
  expect_identical(anc_a$useMacsMut, FALSE)
  expect_equal(length(anc_a$tables), nChr)
  
  # simAnc is ancestry-only: no sites/mutations should be present yet.
  for (chr in seq_len(nChr)) {
    c0 <- table_counts(anc_a$tables[[chr]])
    expect_equal(c0$num_sites, 0L)
    expect_equal(c0$num_mutations, 0L)
    expect_gt(c0$num_trees, 0L)
    expect_gt(c0$num_nodes, 0L)
    expect_gt(c0$num_edges, 0L)
  }
  
  mut_a <- AlphaSimR:::simMut(anc_a, dTheta = dTheta, seed = mut_seed)
  mut_b <- AlphaSimR:::simMut(anc_b, dTheta = dTheta, seed = mut_seed)
  
  expect_identical(mut_a$stage, "simMut")
  expect_identical(mut_a$mutationMode, "postTs")
  expect_equal(as.integer(mut_a$mutationSeed), mut_seed)
  
  keys_a <- lapply(mut_a$tables, ts_variant_keys)
  keys_b <- lapply(mut_b$tables, ts_variant_keys)
  expect_identical(keys_a, keys_b)
  
  for (chr in seq_len(nChr)) {
    c1 <- table_counts(mut_a$tables[[chr]])
    expect_gt(c1$num_mutations, 0L)
    expect_equal(c1$num_sites, c1$num_mutations)
  }
})

test_that("simMut with zero dTheta leaves ancestry tables unchanged", {
  args <- "8 5000 -t 1e-3 -r 1e-4 -s "
  anc <- AlphaSimR:::simAnc(
    args = args,
    nChr = 1L,
    inbred = FALSE,
    ploidy = 2L,
    nThreads = 1L,
    seed = as.integer(42),
    usePhysicalPositions = FALSE,
    Nref = NA_real_
  )
  
  before <- table_counts(anc$tables[[1]])
  out <- AlphaSimR:::simMut(anc, dTheta = 0, seed = as.integer(99))
  after <- table_counts(out$tables[[1]])
  
  expect_identical(before, after)
  expect_equal(after$num_sites, 0L)
  expect_equal(after$num_mutations, 0L)
})

test_that("post-TS mutations are placed within edge span and branch-time bounds", {
  args <- "8 5000 -t 1e-3 -r 1e-4 -s "
  anc <- AlphaSimR:::simAnc(
    args = args,
    nChr = 1L,
    inbred = FALSE,
    ploidy = 2L,
    nThreads = 1L,
    seed = as.integer(88),
    usePhysicalPositions = TRUE,
    Nref = NA_real_
  )
  
  out <- AlphaSimR:::simMut(anc, dTheta = 150, seed = as.integer(77))
  tc <- RcppTskit::TableCollection$new(xptr = out$tables[[1]])
  
  n_mut <- as.integer(tc$num_mutations())
  n_edge <- as.integer(tc$num_edges())
  n_node <- as.integer(tc$num_nodes())
  n_site <- as.integer(tc$num_sites())
  seq_len <- as.numeric(tc$sequence_length())
  
  expect_gt(n_mut, 0L)
  expect_equal(n_site, n_mut)
  
  edge_child <- integer(n_edge)
  edge_parent <- integer(n_edge)
  edge_left <- numeric(n_edge)
  edge_right <- numeric(n_edge)
  for (i in seq_len(n_edge)) {
    e <- tc$edge_table_get_row(i - 1L)
    edge_child[i] <- as.integer(e$child)
    edge_parent[i] <- as.integer(e$parent)
    edge_left[i] <- as.numeric(e$left)
    edge_right[i] <- as.numeric(e$right)
  }
  
  node_time <- numeric(n_node)
  for (i in seq_len(n_node)) {
    node_time[i] <- as.numeric(tc$node_table_get_row(i - 1L)$time)
  }
  
  site_pos <- numeric(n_site)
  for (i in seq_len(n_site)) {
    site_pos[i] <- as.numeric(tc$site_table_get_row(i - 1L)$position)
  }
  
  for (i in seq_len(n_mut)) {
    m <- tc$mutation_table_get_row(i - 1L)
    child <- as.integer(m$node)
    site <- as.integer(m$site)
    mut_time <- as.numeric(m$time)
    pos <- site_pos[site + 1L]
    
    expect_true(pos >= 0)
    expect_true(pos < seq_len)
    
    idx <- which(edge_child == child & edge_left <= pos & pos < edge_right)
    expect_true(length(idx) > 0L)
    
    child_time <- node_time[child + 1L]
    parent_time <- node_time[edge_parent[idx] + 1L]
    expect_true(any(mut_time > child_time & mut_time < parent_time))
  }
})

test_that("Nref rescales node times and sets TS time_units to generations", {
  args <- "8 5000 -t 1e-3 -r 1e-4 -s "
  seed <- as.integer(12345)
  nref <- 10000
  scale <- 4 * nref
  
  anc_unit <- AlphaSimR:::simAnc(
    args = args,
    nChr = 1L,
    inbred = FALSE,
    ploidy = 2L,
    nThreads = 1L,
    seed = seed,
    usePhysicalPositions = FALSE,
    Nref = NA_real_
  )
  anc_gen <- AlphaSimR:::simAnc(
    args = args,
    nChr = 1L,
    inbred = FALSE,
    ploidy = 2L,
    nThreads = 1L,
    seed = seed,
    usePhysicalPositions = FALSE,
    Nref = nref
  )
  
  tc_unit <- RcppTskit::TableCollection$new(xptr = anc_unit$tables[[1]])
  tc_gen <- RcppTskit::TableCollection$new(xptr = anc_gen$tables[[1]])
  
  expect_equal(anc_unit$timeScale, 1)
  expect_equal(anc_gen$timeScale, scale)
  expect_identical(tc_unit$time_units(), "unknown")
  expect_identical(tc_gen$time_units(), "generations")
  
  t_unit <- node_times(tc_unit)
  t_gen <- node_times(tc_gen)
  expect_equal(length(t_unit), length(t_gen))
  expect_true(isTRUE(all.equal(t_gen, t_unit * scale, tolerance = 1e-8)))
  
  c_unit <- table_counts(anc_unit$tables[[1]])
  c_gen <- table_counts(anc_gen$tables[[1]])
  expect_equal(c_unit$num_nodes, c_gen$num_nodes)
  expect_equal(c_unit$num_edges, c_gen$num_edges)
  expect_equal(c_unit$num_trees, c_gen$num_trees)
})
