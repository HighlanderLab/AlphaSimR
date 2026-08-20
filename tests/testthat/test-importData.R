context("importData")

test_that("importTrait",{
  # Create haplotype data
  haplo = rbind(c(1,1,0,1,0),
                c(1,1,0,1,0),
                c(0,1,1,0,0),
                c(0,1,1,0,0))
  colnames(haplo) = letters[1:5]
  
  # Create genetic map
  genMap = data.frame(markerName=letters[1:5],
                      chromosome=c(1,1,1,2,2),
                      position=c(0,0.5,1,0.15,0.4))
  
  # Create pedigree
  ped = data.frame(id=c("a","b"),
                   mother=c(0,0),
                   father=c(0,0))
  
  # Generate an external trait
  myTrait = data.frame(marker = c("a","c","d"),
                       a = c(1,-1,1))
  
  founderPop = importHaplo(haplo = haplo, 
                           genMap = genMap,
                           ploidy = 2L,
                           ped = ped)
  
  SP = SimParam$new(founderPop=founderPop)
  SP$nThreads = 1L
  
  # Import trait
  SP$importTrait(markerNames = myTrait$marker,
                 addEff = myTrait$a,
                 name = "myTrait")
  
  expect_equal(SP$traits[[1]]@addEff, myTrait$a, tolerance=1e-6)
  
  expect_equal(SP$traits[[1]]@intercept, 0, tolerance=1e-6)
  
  pop = newPop(founderPop,simParam=SP)
  
  expect_equal(unname(pop@gv[1,1]), 3, tolerance=1e-6)
  
  expect_equal(unname(pop@gv[2,1]), -3, tolerance=1e-6)
})

test_that("importVCF streams phased VCF and filters unusable sites",{
  vcfFile = tempfile(fileext=".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
          "INFO", "FORMAT", "s1", "s2", sep="\t"),
    paste("1", "10", "rs1", "A", "C", ".", "PASS", ".", "GT",
          "0|1", "1|0", sep="\t"),
    paste("1", "20", "rsMulti", "A", "C,G", ".", "PASS", ".", "GT",
          "0|1", "1|0", sep="\t"),
    paste("1", "30", "rsMissing", "A", "C", ".", "PASS", ".", "GT",
          "0|.", "1|1", sep="\t"),
    paste("1", "40", "rsBadAllele", "A", "C", ".", "PASS", ".", "GT",
          "0|2", "1|1", sep="\t"),
    paste("1", "50", "rs5", "A", "C", ".", "PASS", ".", "GT",
          "1|1", "0|0", sep="\t"),
    paste("1", "60", "rsMono", "A", "C", ".", "PASS", ".", "GT",
          "0|0", "0|0", sep="\t"),
    paste("2", "5", "rs6", "A", "C", ".", "PASS", ".", "GT",
          "0|0", "0|1", sep="\t")
  ), vcfFile)

  out = importVCF(vcfFile=vcfFile,
                  breaks=list("1"=c(0, 100), "2"=c(0, 10)),
                  rates=list("1"=c(0.01), "2"=c(0.1)),
                  returnMeta=TRUE)
  pop = out$pop

  expect_true(isNamedMapPop(pop))
  expect_equal(pop@id, c("s1", "s2"))
  expect_equal(pop@nInd, 2L)
  expect_equal(pop@ploidy, 2L)
  expect_equal(pop@nLoci, c(2L, 1L))
  expect_equal(names(pop@genMap), c("1", "2"))
  expect_equal(unname(pop@genMap[["1"]]), c(0, 0.4), tolerance=1e-12)
  expect_equal(unname(pop@genMap[["2"]]), 0, tolerance=1e-12)
  expect_equal(out$keptPos[["1"]], c(10, 50))
  expect_equal(out$keptPos[["2"]], 5)
  expect_equal(out$stats$skippedNonBiallelic, 1L)
  expect_equal(out$stats$skippedMissing, 1L)
  expect_equal(out$stats$skippedInvalidAllele, 1L)
  expect_equal(out$stats$skippedNonSegregating, 1L)

  H = pullSegSiteHaplo(pop, nThreads=1L)
  expected = rbind(c(0, 1, 0),
                   c(1, 1, 0),
                   c(1, 0, 0),
                   c(0, 0, 1))
  rownames(expected) = c("s1_1", "s1_2", "s2_1", "s2_2")
  colnames(expected) = c("rs1", "rs5", "rs6")
  expect_equal(H, expected)
})

test_that("importVCF rejects unphased genotype calls by default",{
  vcfFile = tempfile(fileext=".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
          "INFO", "FORMAT", "s1", "s2", sep="\t"),
    paste("1", "10", "rs1", "A", "C", ".", "PASS", ".", "GT",
          "0/1", "1|0", sep="\t")
  ), vcfFile)

  expect_error(
    importVCF(vcfFile=vcfFile,
              breaks=c(0, 100),
              rates=c(0.01)),
    regexp="unphased"
  )
})

test_that("importVCF reservoir samples per chromosome",{
  vcfFile = tempfile(fileext=".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
          "INFO", "FORMAT", "s1", "s2", sep="\t"),
    paste("1", "10", "rs1", "A", "C", ".", "PASS", ".", "GT",
          "0|1", "1|0", sep="\t"),
    paste("1", "20", "rs2", "A", "C", ".", "PASS", ".", "GT",
          "1|0", "0|1", sep="\t"),
    paste("2", "10", "rs3", "A", "C", ".", "PASS", ".", "GT",
          "0|1", "1|0", sep="\t"),
    paste("2", "20", "rs4", "A", "C", ".", "PASS", ".", "GT",
          "1|0", "0|1", sep="\t")
  ), vcfFile)

  out = importVCF(vcfFile=vcfFile,
                  breaks=list("1"=c(0, 100), "2"=c(0, 100)),
                  rates=list("1"=c(0.01), "2"=c(0.01)),
                  segSites=c("1"=1L, "2"=1L),
                  siteSamplingSeed=12L,
                  returnMeta=TRUE)

  expect_equal(out$pop@nLoci, c(1L, 1L))
  expect_equal(out$stats$keptByChr, c("1"=1L, "2"=1L))
})

test_that("importVCF can initialize TS forward recorder tables",{
  vcfFile = tempfile(fileext=".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
          "INFO", "FORMAT", "s1", "s2", sep="\t"),
    paste("1", "10", "rs1", "A", "C", ".", "PASS", ".", "GT",
          "0|1", "1|0", sep="\t"),
    paste("1", "50", "rs2", "A", "C", ".", "PASS", ".", "GT",
          "1|1", "0|0", sep="\t")
  ), vcfFile)

  out = importVCF(vcfFile=vcfFile,
                  breaks=list("1"=c(0, 100)),
                  rates=list("1"=c(0.01)),
                  tsRecord=TRUE,
                  returnMeta=TRUE)
  outAlias = importVCF(vcfFile=vcfFile,
                       breaks=list("1"=c(0, 100)),
                       rates=list("1"=c(0.01)),
                       tsRecorde=TRUE,
                       returnMeta=TRUE)
  outNoMut = importVCF(vcfFile=vcfFile,
                       breaks=list("1"=c(0, 100)),
                       rates=list("1"=c(0.01)),
                       tsRecord=TRUE,
                       addTsMut=FALSE,
                       returnMeta=TRUE)
  founder = out$pop

  expect_true(!is.null(attr(founder, "tsForwardSource", exact=TRUE)))
  expect_true(!is.null(attr(founder, "tsForwardPosMeta", exact=TRUE)))
  expect_length(out$tsTables, 1L)
  expect_length(outAlias$tsTables, 1L)
  expect_length(outNoMut$tsTables, 1L)

  summary = rtsk_table_collection_summary2(out$tsTables[[1L]])
  expect_equal(as.integer(summary$num_individuals), 2L)
  expect_equal(as.integer(summary$num_nodes), 4L)
  expect_equal(as.integer(summary$num_sites), 2L)
  expect_equal(as.integer(summary$num_mutations), 4L)

  summaryNoMut = rtsk_table_collection_summary2(outNoMut$tsTables[[1L]])
  expect_equal(as.integer(summaryNoMut$num_individuals), 2L)
  expect_equal(as.integer(summaryNoMut$num_nodes), 4L)
  expect_equal(as.integer(summaryNoMut$num_sites), 0L)
  expect_equal(as.integer(summaryNoMut$num_mutations), 0L)

  SP = SimParam$new(founder)
  SP$nThreads = 1L
  SP$quadProb = 0
  pop0 = newPop(founder, simParam=SP)
  expect_silent(SP$setTrackTs(TRUE, founderPop=pop0))
  expect_true(AlphaSimR:::tsForwardHasRecorder(SP))

  pop1 = makeCross(pop0, matrix(c(1, 2), ncol=2), nProgeny=1,
                   simParam=SP)
  expect_equal(pop1@nInd, 1L)

  paths = AlphaSimR:::tsForwardFinalizeFromSimParam(
    SP,
    out_dir=tempdir(),
    out_basename=paste0("vcf_ts_forward_", Sys.getpid()),
    clear=TRUE
  )
  expect_true(file.exists(paths[[1L]]))
  expect_false(AlphaSimR:::tsForwardHasRecorder(SP))
  unlink(paths, force=TRUE)

  SPnoMut = SimParam$new(outNoMut$pop)
  SPnoMut$nThreads = 1L
  SPnoMut$quadProb = 0
  pop0NoMut = newPop(outNoMut$pop, simParam=SPnoMut)
  expect_silent(SPnoMut$setTrackTs(TRUE, founderPop=pop0NoMut))
  pop1NoMut = makeCross(pop0NoMut, matrix(c(1, 2), ncol=2), nProgeny=1,
                        simParam=SPnoMut)
  expect_equal(pop1NoMut@nInd, 1L)
  pathsNoMut = AlphaSimR:::tsForwardFinalizeFromSimParam(
    SPnoMut,
    out_dir=tempdir(),
    out_basename=paste0("vcf_ts_forward_no_mut_", Sys.getpid()),
    clear=TRUE
  )
  expect_true(file.exists(pathsNoMut[[1L]]))
  unlink(pathsNoMut, force=TRUE)
})
