context("TS forward recorder")

skip_if_not_installed("RcppTskit")

.ts_forward_test_chr_info <- function() {
  ts_path <- testthat::test_path("..", "..", "dev", "testData", "msprime_chr0.trees")
  skip_if(!file.exists(ts_path), "Missing test fixture dev/testData/msprime_chr0.trees")
  list(
    list(
      ts_path = ts_path,
      breaks = c(0, 1e6),
      rates = c(1e-8),
      segSites = 60L
    )
  )
}

.ts_forward_test_setup <- function() {
  chr_info <- .ts_forward_test_chr_info()
  founder <- AlphaSimR:::asMapPop(chr_info = chr_info, inbred = FALSE, ploidy = 2L)
  SP <- SimParam$new(founder)
  SP$nThreads <- 1L
  SP$setTrackRecGen(TRUE)
  SP$quadProb <- 0
  pop0 <- newPop(founder, simParam = SP)
  list(chr_info = chr_info, SP = SP, pop0 = pop0)
}

.normalize_edge_df <- function(df) {
  out <- data.frame(
    chr = as.integer(df$chr),
    child = as.integer(df$child),
    parent = as.integer(df$parent),
    left = as.numeric(df$left),
    right = as.numeric(df$right),
    stringsAsFactors = FALSE
  )
  out <- out[out$right > out$left, , drop = FALSE]
  if (nrow(out) == 0L) {
    return(out)
  }
  ord <- order(out$chr, out$child, out$parent, out$left, out$right)
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL
  out
}

.collect_forward_edges_for_children <- function(recorder, chr, child_iids) {
  chrState <- recorder$chr[[chr]]
  tc <- chrState$tc
  n_edge <- as.integer(tc$num_edges())
  out <- vector("list", n_edge)
  k <- 0L

  iid_keys <- ls(chrState$indMap, all.names = TRUE)
  if (length(iid_keys) == 0L) {
    return(data.frame(chr = integer(), child = integer(), parent = integer(),
                      left = numeric(), right = numeric(), stringsAsFactors = FALSE))
  }
  ind_rows <- vapply(iid_keys, function(key) {
    as.integer(get(key, envir = chrState$indMap, inherits = FALSE))
  }, integer(1))
  ind_row_to_iid <- setNames(as.integer(iid_keys), as.character(ind_rows))

  for (i in seq_len(n_edge)) {
    e <- tc$edge_table_get_row(i - 1L)
    child_node <- as.integer(e$child)
    parent_node <- as.integer(e$parent)
    child_ind_row <- as.integer(tc$node_table_get_row(child_node)$individual)
    parent_ind_row <- as.integer(tc$node_table_get_row(parent_node)$individual)

    child_iid <- ind_row_to_iid[as.character(child_ind_row)]
    parent_iid <- ind_row_to_iid[as.character(parent_ind_row)]
    if (length(child_iid) == 0L || length(parent_iid) == 0L ||
        is.na(child_iid) || is.na(parent_iid)) {
      next
    }
    child_iid <- as.integer(child_iid)
    if (!(child_iid %in% child_iids)) {
      next
    }

    k <- k + 1L
    out[[k]] <- data.frame(
      chr = as.integer(chr),
      child = child_iid,
      parent = as.integer(parent_iid),
      left = as.numeric(e$left),
      right = as.numeric(e$right),
      stringsAsFactors = FALSE
    )
  }

  if (k == 0L) {
    return(data.frame(chr = integer(), child = integer(), parent = integer(),
                      left = numeric(), right = numeric(), stringsAsFactors = FALSE))
  }
  do.call(rbind, out[seq_len(k)])
}

.collect_bridge_edges_for_children <- function(SP, sim_output, chr_info, child_iids) {
  pos_list <- attr(sim_output[[1L]], "tsForwardPosMeta", exact = TRUE)$posList
  bridge_list <- AlphaSimR:::bridgeCollectSegGenFromSimOutput(
    SP,
    sim_output,
    chr_info = chr_info,
    pos_list = pos_list
  )
  if (length(bridge_list) == 0L) {
    return(data.frame(chr = integer(), child = integer(), parent = integer(),
                      left = numeric(), right = numeric(), stringsAsFactors = FALSE))
  }
  bridge_df <- do.call(rbind, bridge_list)
  bridge_edges <- data.frame(
    chr = as.integer(bridge_df$chr),
    child = as.integer(bridge_df$childId),
    parent = as.integer(bridge_df$parentId),
    left = as.numeric(bridge_df$left),
    right = as.numeric(bridge_df$right),
    stringsAsFactors = FALSE
  )
  bridge_edges <- bridge_edges[bridge_edges$child %in% child_iids, , drop = FALSE]
  .normalize_edge_df(bridge_edges)
}

.ts_forward_metadata_text <- function(x) {
  if (is.raw(x)) {
    return(rawToChar(x))
  }
  paste0(as.character(x), collapse = "")
}

test_that("tsForward recorder lifecycle init->append->finalize works", {
  old_opts <- options(AlphaSimR.tsForwardKeepSeg = FALSE)
  on.exit(options(old_opts), add = TRUE)

  chr_info <- .ts_forward_test_chr_info()
  chr_info[[1L]]$breaks <- c(0, 1e6 / 3, 2e6 / 3, 1e6)
  chr_info[[1L]]$rates <- c(5e-7, 5e-6, 5e-7)
  founder <- AlphaSimR:::asMapPop(chr_info = chr_info, inbred = FALSE, ploidy = 2L)
  SP <- SimParam$new(founder)
  SP$nThreads <- 1L
  SP$setTrackRecGen(TRUE)
  SP$quadProb <- 0
  pop0 <- newPop(founder, simParam = SP)

  expect_false(AlphaSimR:::tsForwardHasRecorder(SP))
  expect_error(
    AlphaSimR:::tsForwardFinalizeFromSimParam(SP),
    regexp = "No tsForwardRecorder attached to simParam"
  )

  expect_silent(AlphaSimR:::tsForwardInitOnSimParam(SP, founderPop = pop0))
  expect_true(AlphaSimR:::tsForwardHasRecorder(SP))

  crossPlan <- matrix(c(1, 2), ncol = 2, byrow = TRUE)
  pop1 <- makeCross(pop0, crossPlan, nProgeny = 2, simParam = SP)
  expect_true(is.null(pop1@misc$tsSegGen))

  out <- AlphaSimR:::tsForwardFinalizeFromSimParam(
    SP,
    out_dir = tempdir(),
    out_basename = paste0("ts_forward_lifecycle_", Sys.getpid()),
    clear = TRUE
  )
  expect_length(out, 1L)
  expect_true(file.exists(out[[1]]))
  expect_false(AlphaSimR:::tsForwardHasRecorder(SP))
  expect_error(
    AlphaSimR:::tsForwardFinalizeFromSimParam(SP),
    regexp = "No tsForwardRecorder attached to simParam"
  )
  unlink(out, force = TRUE)
})

test_that("metadata-on child nodes use AlphaSimR node metadata", {
  old_opts <- options(
    AlphaSimR.tsForwardKeepSeg = FALSE,
    AlphaSimR.tsForwardKeepRecHistGen = FALSE,
    AlphaSimR.tsForwardAttachMetadata = TRUE
  )
  on.exit(options(old_opts), add = TRUE)

  expect_true(exists("tsForwardNodeTableAddRowsWithMetadata", where = asNamespace("AlphaSimR")))

  st <- .ts_forward_test_setup()
  SP <- st$SP
  pop0 <- st$pop0

  AlphaSimR:::tsForwardInitOnSimParam(SP, founderPop = pop0)
  pop1 <- makeCross(pop0, matrix(c(1, 2), ncol = 2, byrow = TRUE), nProgeny = 1, simParam = SP)

  rec <- attr(SP, "tsForwardRecorder", exact = TRUE)
  child_key <- AlphaSimR:::.tsForwardNodeKey(pop1@iid[[1]], 1L)
  child_node <- get(child_key, envir = rec$chr[[1]]$nodeMap, inherits = FALSE)
  node_row <- rec$chr[[1]]$tc$node_table_get_row(as.integer(child_node))

  expect_equal(
    .ts_forward_metadata_text(node_row$metadata),
    paste0("{\"alphaSimR\":{\"id\":\"", child_key, "\"}}")
  )
})

test_that("setTrackTs convenience method initializes and clears recorder", {
  old_opts <- options(AlphaSimR.tsForwardKeepSeg = FALSE)
  on.exit(options(old_opts), add = TRUE)

  st <- .ts_forward_test_setup()
  SP <- st$SP
  pop0 <- st$pop0
  chr_info <- st$chr_info

  expect_error(
    SP$setTrackTs(TRUE),
    regexp = "requires founderPop"
  )

  expect_silent(SP$setTrackTs(TRUE, founderPop = pop0))
  expect_true(SP$isTrackRecGen)
  expect_true(AlphaSimR:::tsForwardHasRecorder(SP))

  pop1 <- makeCross(pop0, matrix(c(1, 2), ncol = 2, byrow = TRUE), nProgeny = 2, simParam = SP)
  expect_true(is.null(pop1@misc$tsSegGen))

  out <- AlphaSimR:::tsForwardFinalizeFromSimParam(
    SP,
    out_dir = tempdir(),
    out_basename = paste0("ts_forward_setTrackTs_", Sys.getpid()),
    clear = TRUE
  )
  expect_length(out, 1L)
  expect_true(file.exists(out[[1]]))
  expect_false(AlphaSimR:::tsForwardHasRecorder(SP))

  expect_silent(SP$setTrackTs(FALSE))
  expect_false(AlphaSimR:::tsForwardHasRecorder(SP))

  unlink(out, force = TRUE)
})

test_that("direct recorder path works with nThreads > 1", {
  old_opts <- options(
    AlphaSimR.tsForwardKeepSeg = FALSE,
    AlphaSimR.tsForwardKeepRecHistGen = FALSE
  )
  on.exit(options(old_opts), add = TRUE)

  st <- .ts_forward_test_setup()
  SP <- st$SP
  pop0 <- st$pop0
  chr_info <- st$chr_info

  SP$nThreads <- 2L
  expect_silent(AlphaSimR:::tsForwardInitOnSimParam(SP, founderPop = pop0))

  crossPlan <- matrix(c(1, 2), ncol = 2, byrow = TRUE)
  pop1 <- makeCross(pop0, crossPlan, nProgeny = 2, simParam = SP, nThreads = 2L)
  expect_true(is.null(pop1@misc$tsSegGen))
  expect_true(all(vapply(SP$recHistGen[pop1@iid], is.null, logical(1))))

  recorder <- attr(SP, "tsForwardRecorder", exact = TRUE)
  expect_true(inherits(recorder, "tsForwardRecorder"))
  forward_edges <- .collect_forward_edges_for_children(recorder, chr = 1L, child_iids = pop1@iid)
  forward_edges <- .normalize_edge_df(forward_edges)
  expect_gt(nrow(forward_edges), 0L)

  out <- AlphaSimR:::tsForwardFinalizeFromSimParam(
    SP,
    out_dir = tempdir(),
    out_basename = paste0("ts_forward_threads_", Sys.getpid()),
    clear = TRUE
  )
  expect_length(out, 1L)
  expect_true(file.exists(out[[1]]))
  expect_false(AlphaSimR:::tsForwardHasRecorder(SP))
  unlink(out, force = TRUE)
})

test_that("forward recorder edges match recHistGen bridge edges for one generation", {
  old_opts <- options(
    AlphaSimR.tsForwardKeepSeg = FALSE,
    AlphaSimR.tsForwardKeepRecHistGen = TRUE
  )
  on.exit(options(old_opts), add = TRUE)

  st <- .ts_forward_test_setup()
  SP <- st$SP
  pop0 <- st$pop0
  chr_info <- st$chr_info

  expect_silent(AlphaSimR:::tsForwardInitOnSimParam(SP, founderPop = pop0))

  crossPlan <- matrix(c(1, 2), ncol = 2, byrow = TRUE)
  pop1 <- makeCross(pop0, crossPlan, nProgeny = 2, simParam = SP)

  recorder <- attr(SP, "tsForwardRecorder", exact = TRUE)
  expect_true(inherits(recorder, "tsForwardRecorder"))
  forward_edges <- .collect_forward_edges_for_children(recorder, chr = 1L, child_iids = pop1@iid)
  forward_edges <- .normalize_edge_df(forward_edges)
  expect_gt(nrow(forward_edges), 0L)

  bridge_edges <- .collect_bridge_edges_for_children(
    SP = SP,
    sim_output = list(pop0, pop1),
    chr_info = chr_info,
    child_iids = pop1@iid
  )

  expect_equal(nrow(forward_edges), nrow(bridge_edges))
  expect_equal(forward_edges[, c("chr", "child", "parent")],
               bridge_edges[, c("chr", "child", "parent")])
  expect_equal(forward_edges$left, bridge_edges$left, tolerance = 1e-6)
  expect_equal(forward_edges$right, bridge_edges$right, tolerance = 1e-6)
})

test_that("forward recorder uses chr_info rate map for non-constant coordinate conversion", {
  old_opts <- options(
    AlphaSimR.tsForwardKeepSeg = FALSE,
    AlphaSimR.tsForwardKeepRecHistGen = TRUE
  )
  on.exit(options(old_opts), add = TRUE)

  chr_info <- .ts_forward_test_chr_info()
  chr_info[[1L]]$breaks <- c(0, 1e6 / 3, 2e6 / 3, 1e6)
  chr_info[[1L]]$rates <- c(5e-7, 5e-6, 5e-7)
  founder <- AlphaSimR:::asMapPop(chr_info = chr_info, inbred = FALSE, ploidy = 2L)
  posMeta <- attr(founder, "tsForwardPosMeta", exact = TRUE)
  expect_equal(posMeta$breaksList[[1L]], chr_info[[1L]]$breaks)
  expect_equal(posMeta$ratesList[[1L]], chr_info[[1L]]$rates)

  SP <- SimParam$new(founder)
  SP$nThreads <- 1L
  SP$setTrackRecGen(TRUE)
  SP$quadProb <- 0
  pop0 <- newPop(founder, simParam = SP)
  AlphaSimR:::tsForwardInitOnSimParam(SP, founderPop = pop0)

  crossPlan <- matrix(c(1, 2), ncol = 2, byrow = TRUE)
  pop1 <- makeCross(pop0, crossPlan, nProgeny = 2, simParam = SP)

  recorder <- attr(SP, "tsForwardRecorder", exact = TRUE)
  forward_edges <- .collect_forward_edges_for_children(recorder, chr = 1L, child_iids = pop1@iid)
  forward_edges <- .normalize_edge_df(forward_edges)
  expect_gt(nrow(forward_edges), 0L)

  bridge_edges <- .collect_bridge_edges_for_children(
    SP = SP,
    sim_output = list(pop0, pop1),
    chr_info = chr_info,
    child_iids = pop1@iid
  )

  expect_equal(nrow(forward_edges), nrow(bridge_edges))
  expect_equal(forward_edges[, c("chr", "child", "parent")],
               bridge_edges[, c("chr", "child", "parent")])
  expect_equal(forward_edges$left, bridge_edges$left, tolerance = 1e-6)
  expect_equal(forward_edges$right, bridge_edges$right, tolerance = 1e-6)
})

test_that("forward recorder edges match bridge edges across generations with nThreads > 1", {
  old_opts <- options(
    AlphaSimR.tsForwardKeepSeg = FALSE,
    AlphaSimR.tsForwardKeepRecHistGen = TRUE
  )
  on.exit(options(old_opts), add = TRUE)

  st <- .ts_forward_test_setup()
  SP <- st$SP
  pop0 <- st$pop0
  chr_info <- st$chr_info

  SP$nThreads <- 2L
  expect_silent(AlphaSimR:::tsForwardInitOnSimParam(SP, founderPop = pop0))

  cp <- matrix(c(1, 2), ncol = 2, byrow = TRUE)
  pop1 <- makeCross(pop0, cp, nProgeny = 2, simParam = SP, nThreads = 2L)
  pop2 <- makeCross(pop1, cp, nProgeny = 2, simParam = SP, nThreads = 2L)
  pop3 <- makeCross(pop2, cp, nProgeny = 2, simParam = SP, nThreads = 2L)

  expect_true(is.null(pop1@misc$tsSegGen))
  expect_true(is.null(pop2@misc$tsSegGen))
  expect_true(is.null(pop3@misc$tsSegGen))

  recorder <- attr(SP, "tsForwardRecorder", exact = TRUE)
  expect_true(inherits(recorder, "tsForwardRecorder"))

  child_iids <- c(pop1@iid, pop2@iid, pop3@iid)
  forward_edges <- .collect_forward_edges_for_children(recorder, chr = 1L, child_iids = child_iids)
  forward_edges <- .normalize_edge_df(forward_edges)
  expect_gt(nrow(forward_edges), 0L)

  bridge_edges <- .collect_bridge_edges_for_children(
    SP = SP,
    sim_output = list(pop0, pop1, pop2, pop3),
    chr_info = chr_info,
    child_iids = child_iids
  )

  expect_equal(nrow(forward_edges), nrow(bridge_edges))
  expect_equal(forward_edges[, c("chr", "child", "parent")],
               bridge_edges[, c("chr", "child", "parent")])
  expect_equal(forward_edges$left, bridge_edges$left, tolerance = 1e-6)
  expect_equal(forward_edges$right, bridge_edges$right, tolerance = 1e-6)
})

test_that("makeCross2 + self threaded direct recorder matches bridge edges", {
  old_opts <- options(
    AlphaSimR.tsForwardKeepSeg = FALSE,
    AlphaSimR.tsForwardKeepRecHistGen = TRUE
  )
  on.exit(options(old_opts), add = TRUE)

  st <- .ts_forward_test_setup()
  SP <- st$SP
  pop0 <- st$pop0
  chr_info <- st$chr_info

  SP$nThreads <- 2L
  expect_silent(AlphaSimR:::tsForwardInitOnSimParam(SP, founderPop = pop0))

  cp <- matrix(c(1, 2), ncol = 2, byrow = TRUE)
  pop1 <- makeCross2(pop0, pop0, cp, nProgeny = 2, simParam = SP, nThreads = 2L)
  pop2 <- self(pop1, nProgeny = 1, keepParents = FALSE, simParam = SP, nThreads = 2L)

  expect_true(is.null(pop1@misc$tsSegGen))
  expect_true(is.null(pop2@misc$tsSegGen))

  recorder <- attr(SP, "tsForwardRecorder", exact = TRUE)
  expect_true(inherits(recorder, "tsForwardRecorder"))

  child_iids <- c(pop1@iid, pop2@iid)
  forward_edges <- .collect_forward_edges_for_children(recorder, chr = 1L, child_iids = child_iids)
  forward_edges <- .normalize_edge_df(forward_edges)
  expect_gt(nrow(forward_edges), 0L)

  bridge_edges <- .collect_bridge_edges_for_children(
    SP = SP,
    sim_output = list(pop0, pop1, pop2),
    chr_info = chr_info,
    child_iids = child_iids
  )

  expect_equal(nrow(forward_edges), nrow(bridge_edges))
  expect_equal(forward_edges[, c("chr", "child", "parent")],
               bridge_edges[, c("chr", "child", "parent")])
  expect_equal(forward_edges$left, bridge_edges$left, tolerance = 1e-6)
  expect_equal(forward_edges$right, bridge_edges$right, tolerance = 1e-6)
})

test_that("without recorder, makeCross suppresses tsSegGen and recHistGen by default", {
  old_opts <- options(
    AlphaSimR.tsForwardKeepSeg = FALSE,
    AlphaSimR.tsForwardKeepRecHistGen = FALSE
  )
  on.exit(options(old_opts), add = TRUE)

  st <- .ts_forward_test_setup()
  SP <- st$SP
  pop0 <- st$pop0

  expect_false(AlphaSimR:::tsForwardHasRecorder(SP))
  crossPlan <- matrix(c(1, 2), ncol = 2, byrow = TRUE)
  pop1 <- makeCross(pop0, crossPlan, nProgeny = 2, simParam = SP)

  expect_true(is.null(pop1@misc$tsSegGen))
  expect_true(length(pop1@iid) > 0L)
  expect_true(all(vapply(SP$recHistGen[pop1@iid], is.null, logical(1))))
})

test_that("without recorder, makeCross keeps recHistGen when explicitly requested", {
  old_opts <- options(
    AlphaSimR.tsForwardKeepSeg = FALSE,
    AlphaSimR.tsForwardKeepRecHistGen = TRUE
  )
  on.exit(options(old_opts), add = TRUE)

  st <- .ts_forward_test_setup()
  SP <- st$SP
  pop0 <- st$pop0

  expect_false(AlphaSimR:::tsForwardHasRecorder(SP))
  crossPlan <- matrix(c(1, 2), ncol = 2, byrow = TRUE)
  pop1 <- makeCross(pop0, crossPlan, nProgeny = 2, simParam = SP)

  expect_true(is.null(pop1@misc$tsSegGen))
  expect_true(length(pop1@iid) > 0L)
  expect_true(all(!vapply(SP$recHistGen[pop1@iid], is.null, logical(1))))
})

test_that("raw tsSegGen append matches data.frame append", {
  old_opts <- options(
    AlphaSimR.tsForwardKeepSeg = TRUE,
    AlphaSimR.tsForwardKeepRecHistGen = TRUE
  )
  on.exit(options(old_opts), add = TRUE)

  st <- .ts_forward_test_setup()
  SP <- st$SP
  pop0 <- st$pop0
  chr_info <- st$chr_info

  crossPlan <- matrix(c(1, 2), ncol = 2, byrow = TRUE)
  nProgeny <- 2L
  pop1 <- makeCross(pop0, crossPlan, nProgeny = nProgeny, simParam = SP)
  seg_df <- pop1@misc$tsSegGen
  expect_false(is.null(seg_df))
  expect_gt(nrow(seg_df), 0L)

  raw_seg <- as.matrix(seg_df[, c("childLocal", "chr", "hap", "parentSide",
                                  "parentIndex", "parentHap", "leftGen", "rightGen")])
  crossPlanExp <- cbind(
    rep(crossPlan[, 1], each = nProgeny),
    rep(crossPlan[, 2], each = nProgeny)
  )
  mother_iid <- pop0@iid[crossPlanExp[, 1]]
  father_iid <- pop0@iid[crossPlanExp[, 2]]

  rec_df <- AlphaSimR:::tsForwardInit(founderPop = pop0)
  rec_df <- AlphaSimR:::tsForwardAppendSeg(rec_df, seg_df, SP)
  edges_df <- .normalize_edge_df(.collect_forward_edges_for_children(rec_df, chr = 1L, child_iids = pop1@iid))

  rec_raw <- AlphaSimR:::tsForwardInit(founderPop = pop0)
  rec_raw <- AlphaSimR:::tsForwardAppendSegGenRaw(
    recorder = rec_raw,
    tsSegGen = raw_seg,
    childIid = pop1@iid,
    motherIid = mother_iid,
    fatherIid = father_iid,
    femaleMap = SP$femaleMap,
    maleMap = SP$maleMap,
    simParam = SP
  )
  edges_raw <- .normalize_edge_df(.collect_forward_edges_for_children(rec_raw, chr = 1L, child_iids = pop1@iid))

  expect_equal(nrow(edges_raw), nrow(edges_df))
  expect_equal(edges_raw[, c("chr", "child", "parent")], edges_df[, c("chr", "child", "parent")])
  expect_equal(edges_raw$left, edges_df$left, tolerance = 1e-6)
  expect_equal(edges_raw$right, edges_df$right, tolerance = 1e-6)
})

test_that("with recorder active and keepRecHistGen FALSE, recHistGen is suppressed", {
  old_opts <- options(
    AlphaSimR.tsForwardKeepSeg = FALSE,
    AlphaSimR.tsForwardKeepRecHistGen = FALSE
  )
  on.exit(options(old_opts), add = TRUE)

  st <- .ts_forward_test_setup()
  SP <- st$SP
  pop0 <- st$pop0
  chr_info <- st$chr_info

  AlphaSimR:::tsForwardInitOnSimParam(SP, founderPop = pop0)
  crossPlan <- matrix(c(1, 2), ncol = 2, byrow = TRUE)
  pop1 <- makeCross(pop0, crossPlan, nProgeny = 2, simParam = SP)

  expect_true(is.null(pop1@misc$tsSegGen))
  expect_true(all(vapply(SP$recHistGen[pop1@iid], is.null, logical(1))))
})

test_that("writeTreesFromSimParam defaults samples to last generation and updates sample flags", {
  old_opts <- options(
    AlphaSimR.tsForwardKeepSeg = FALSE,
    AlphaSimR.tsForwardKeepRecHistGen = FALSE
  )
  on.exit(options(old_opts), add = TRUE)

  st <- .ts_forward_test_setup()
  SP <- st$SP
  pop0 <- st$pop0
  chr_info <- st$chr_info

  AlphaSimR:::tsForwardInitOnSimParam(SP, founderPop = pop0)
  cp <- matrix(c(1, 2), ncol = 2, byrow = TRUE)
  pop1 <- makeCross(pop0, cp, nProgeny = 2, simParam = SP)
  pop2 <- makeCross(pop1, cp, nProgeny = 2, simParam = SP)
  rec <- attr(SP, "tsForwardRecorder", exact = TRUE)
  n_past_nodes <- length(AlphaSimR:::.tsForwardPastIndividualNodes(rec$chr[[1]]$tc))

  out <- AlphaSimR:::tsForwardWriteTreesFromSimParam(
    SP,
    out_dir = tempdir(),
    out_basename = paste0("ts_forward_samples_default_", Sys.getpid()),
    simplify = FALSE,
    clear = FALSE
  )
  expect_length(out, 1L)
  expect_true(file.exists(out[[1]]))

  ts <- RcppTskit::ts_load(out[[1]])
  expect_equal(length(ts$samples()), (pop0@nInd + pop2@nInd) * pop2@ploidy)
  unlink(out, force = TRUE)
})

test_that("writeTreesFromSimParam accepts manual samples and simplify", {
  old_opts <- options(
    AlphaSimR.tsForwardKeepSeg = FALSE,
    AlphaSimR.tsForwardKeepRecHistGen = FALSE
  )
  on.exit(options(old_opts), add = TRUE)

  st <- .ts_forward_test_setup()
  SP <- st$SP
  pop0 <- st$pop0
  chr_info <- st$chr_info

  AlphaSimR:::tsForwardInitOnSimParam(SP, founderPop = pop0)
  cp <- matrix(c(1, 2), ncol = 2, byrow = TRUE)
  pop1 <- makeCross(pop0, cp, nProgeny = 2, simParam = SP)

  rec <- attr(SP, "tsForwardRecorder", exact = TRUE)
  keys <- paste0(pop1@iid[[1]], "_", seq_len(pop1@ploidy))
  manual_samples <- as.integer(vapply(keys, function(k) {
    get(k, envir = rec$chr[[1]]$nodeMap, inherits = FALSE)
  }, integer(1)))

  out <- AlphaSimR:::tsForwardWriteTreesFromSimParam(
    SP,
    out_dir = tempdir(),
    out_basename = paste0("ts_forward_samples_manual_", Sys.getpid()),
    simplify = TRUE,
    samples = manual_samples,
    keep_existing_samples = FALSE,
    clear = TRUE
  )
  expect_length(out, 1L)
  expect_true(file.exists(out[[1]]))
  expect_false(AlphaSimR:::tsForwardHasRecorder(SP))

  ts <- RcppTskit::ts_load(out[[1]])
  expect_equal(length(ts$samples()), length(manual_samples))
  unlink(out, force = TRUE)
})

test_that("forward simplify update_sample_flags defaults to false", {
  expect_identical(formals(AlphaSimR:::tsForwardFinalize)$update_sample_flags, FALSE)
  expect_identical(formals(AlphaSimR:::tsForwardFinalizeFromSimParam)$update_sample_flags, FALSE)
  expect_identical(formals(AlphaSimR:::tsForwardWriteTreesFromSimParam)$update_sample_flags, FALSE)
})

test_that("default simplify samples keep imported past pedigree nodes", {
  old_opts <- options(
    AlphaSimR.tsForwardKeepSeg = FALSE,
    AlphaSimR.tsForwardKeepRecHistGen = FALSE
  )
  on.exit(options(old_opts), add = TRUE)

  st <- .ts_forward_test_setup()
  SP <- st$SP
  pop0 <- st$pop0

  AlphaSimR:::tsForwardInitOnSimParam(SP, founderPop = pop0)
  cp <- matrix(c(1, 2), ncol = 2, byrow = TRUE)
  pop1 <- makeCross(pop0, cp, nProgeny = 2, simParam = SP)
  pop2 <- makeCross(pop1, cp, nProgeny = 2, simParam = SP)
  rec <- attr(SP, "tsForwardRecorder", exact = TRUE)
  n_past_nodes <- length(AlphaSimR:::.tsForwardPastIndividualNodes(rec$chr[[1]]$tc))
  expect_gt(n_past_nodes, 0L)

  out <- AlphaSimR:::tsForwardWriteTreesFromSimParam(
    SP,
    out_dir = tempdir(),
    out_basename = paste0("ts_forward_samples_negative_time_", Sys.getpid()),
    simplify = TRUE,
    keep_unary = TRUE,
    clear = TRUE
  )
  expect_length(out, 1L)
  expect_true(file.exists(out[[1]]))

  ts <- RcppTskit::ts_load(out[[1]])
  tc <- ts$dump_tables()
  n_node <- as.integer(tc$num_nodes())
  node_rows <- lapply(seq_len(n_node) - 1L, function(nid) {
    tc$node_table_get_row(nid)
  })
  node_times <- vapply(node_rows, function(x) as.numeric(x$time), numeric(1))
  node_ind <- vapply(node_rows, function(x) as.integer(x$individual), integer(1))

  expect_gte(sum(node_times > 0 & node_ind != -1L), n_past_nodes)
  expect_gte(sum(abs(node_times + 2) < 1e-8 & node_ind != -1L), pop2@nInd * pop2@ploidy)
  unlink(out, force = TRUE)
})

test_that("simplify keeps original samples but not extra past nodes when requested", {
  old_opts <- options(
    AlphaSimR.tsForwardKeepSeg = FALSE,
    AlphaSimR.tsForwardKeepRecHistGen = FALSE
  )
  on.exit(options(old_opts), add = TRUE)

  st <- .ts_forward_test_setup()
  SP <- st$SP
  pop0 <- st$pop0

  AlphaSimR:::tsForwardInitOnSimParam(SP, founderPop = pop0)
  cp <- matrix(c(1, 2), ncol = 2, byrow = TRUE)
  pop1 <- makeCross(pop0, cp, nProgeny = 2, simParam = SP)
  pop2 <- makeCross(pop1, cp, nProgeny = 2, simParam = SP)
  rec <- attr(SP, "tsForwardRecorder", exact = TRUE)
  n_past_nodes <- length(AlphaSimR:::.tsForwardPastIndividualNodes(rec$chr[[1]]$tc))

  out <- AlphaSimR:::tsForwardWriteTreesFromSimParam(
    SP,
    out_dir = tempdir(),
    out_basename = paste0("ts_forward_samples_replace_", Sys.getpid()),
    simplify = TRUE,
    keep_unary = TRUE,
    keep_existing_samples = FALSE,
    clear = TRUE
  )
  expect_length(out, 1L)
  expect_true(file.exists(out[[1]]))

  ts <- RcppTskit::ts_load(out[[1]])
  tc <- ts$dump_tables()
  n_node <- as.integer(tc$num_nodes())
  node_rows <- lapply(seq_len(n_node) - 1L, function(nid) {
    tc$node_table_get_row(nid)
  })
  node_times <- vapply(node_rows, function(x) as.numeric(x$time), numeric(1))
  node_ind <- vapply(node_rows, function(x) as.integer(x$individual), integer(1))
  node_flags <- vapply(node_rows, function(x) as.integer(x$flags), integer(1))
  is_sample <- bitwAnd(node_flags, 1L) != 0L

  expect_lt(sum(node_times > 0 & node_ind != -1L), n_past_nodes)
  expect_gte(sum(abs(node_times + 2) < 1e-8 & node_ind != -1L), pop2@nInd * pop2@ploidy)
  expect_gte(sum(is_sample & node_times == 0 & node_ind != -1L), pop0@nInd * pop0@ploidy)
  expect_equal(sum(is_sample & node_times > 0 & node_ind != -1L), 0L)
  expect_equal(sum(is_sample & node_times < 0 & node_ind != -1L), 0L)
  unlink(out, force = TRUE)
})

test_that("build tree-sequence objects then write all chromosomes", {
  st <- .ts_forward_test_setup()
  SP <- st$SP
  pop0 <- st$pop0

  AlphaSimR:::tsForwardInitOnSimParam(SP, founderPop = pop0)
  cp <- matrix(c(1, 2), ncol = 2, byrow = TRUE)
  pop1 <- makeCross(pop0, cp, nProgeny = 2, simParam = SP)
  expect_true(is.null(pop1@misc$tsSegGen))

  rec <- attr(SP, "tsForwardRecorder", exact = TRUE)
  indTime <- AlphaSimR:::.tsForwardUpdateIndTimeCache(rec$indTime, SP$pedigree)
  rec$indTime <- indTime
  ts_list <- AlphaSimR:::tsForwardFinalize(
    recorder = rec,
    simplify = TRUE,
    keep_unary = TRUE,
    indTime = indTime,
    update_samples = TRUE
  )
  expect_length(ts_list, 1L)
  expect_gt(length(ts_list[[1]]$samples()), 0L)

  out <- AlphaSimR:::tsForwardWriteTreeSequences(
    ts_list = ts_list,
    out_dir = tempdir(),
    out_basename = paste0("ts_forward_build_write_", Sys.getpid()),
    recorder = rec
  )
  expect_length(out, 1L)
  expect_true(file.exists(out[[1]]))
  unlink(out, force = TRUE)
})
