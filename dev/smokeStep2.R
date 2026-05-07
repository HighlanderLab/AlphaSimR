library(AlphaSimR)
use_virtualenv("~/r-reticulate-env", required = TRUE)
tskit <- import("tskit")
devtools::load_all()

# two chromosomes
L1 <- 1e6
L2 <- 2e6
chr_info <- list(
  list(ts_path="/Users/jliang2/R_scripts/AlphaSimR_test/dev/testData/msprime_chr0.trees",
       breaks=c(0, L1/2, L1), rates=c(1e-8, 2e-8), segSites=60),
  #breaks=c(0, L1/2, L1), rates=c(1e-5, 2e-5), segSites=60),
  list(ts_path="/Users/jliang2/R_scripts/AlphaSimR_test/dev/testData/msprime_chr1.trees",
       breaks=c(0, L2/3, 2*L2/3, L2), rates=c(1e-7, 1e-8, 1e-7), segSites=155)
  #breaks=c(0, L2/3, 2*L2/3, L2), rates=c(1e-4, 1e-5, 1e-4), segSites=155)
)

founderGenomes1 <- asMapPop(chr_info = chr_info, inbred=FALSE, ploidy=2L)

set.seed(42)
SP = SimParam$new(founderGenomes1)
SP$setSexes("yes_sys")
SP$addTraitA(nQtlPerChr = 5,
             mean = 500,
             var = 450)

SP$setTrackPed(TRUE)
# try the new function here, it automatically set setTrackRec also.
SP$setTrackRecGen(TRUE)
SP$recHistGen
basePop = newPop(founderGenomes1)
# the 2 objects are same now:
SP$recHistGen
SP$recHist
basePop = setPheno(basePop,
                   h2 = 0.5)

#--- n generations
nCycles<-2

# very simple container for each cycles sim output
simOutput<-list(basePop)
cycle<-1
for(cycle in 1:nCycles){
  cat(paste0(" C",cycle))
  # choose the best from last cycle
  chosenParents<- selectInd(pop=simOutput[[cycle]], nInd=6, use = "gv")
  # make crosses
  offspringPop<-randCross(pop=chosenParents, nCrosses=2, nProgeny = 5)
  # phenotype  new offspring
  offspringPop<-setPheno(pop = offspringPop, h2 = 0.5)
  # add new offspring to simOutput list
  simOutput[[cycle+1]]<-offspringPop
}

# see the difference between recHist and recHistGen
RHG <- SP$recHistGen
RH <- SP$recHist
# ind 3; chr 2; hap 1. Maybe not the same output, please check RHG and RH to find a hap with recombination
rh <- RH[[3]][[2]][[1]]
rhg <- RHG[[3]][[2]][[1]]
gm <- SP$genMap[[2]]
rh
rhg
gm[[111]]
gm[[112]]


# for RecHist
# bridgeSegDfList store the indexes of SNPs after recombination events
bridgeCollectSegFromSimOutput(SP, simOutput)
# for RecHistGen
# bridgeSegDfListGen store the positions of where recombination happen
bridgeCollectSegGenFromSimOutput(SP, simOutput)

# for RecHist
edgeDf <- bridgeAllSegToEdgeDf(chr_info)
bridgeWriteTrees(chr_info, edgeDf, SP)
# load the tree in Python...
#origin = tskit.load('/Users/jliang2/R_scripts/AlphaSimR_test/dev/testData/msprime_chr0.trees')
#marker_ts = tskit.load('/Users/jliang2/R_scripts/AlphaSimR_test/dev/testData/AlphaSimR_extended_chr0.trees')
# check the number of trees, nodes, and individual

# for RecHistGen
bridgeWriteTrees(chr_info, do.call(rbind, bridgeSegDfListGen), SP)
# real_break_ts = tskit.load('/Users/jliang2/R_scripts/AlphaSimR_test/dev/testData/AlphaSimR_extended_chr0.trees')


L1 <- 1e6
L2 <- 2e6
chr_info <- list(
  list(ts_path="/Users/jliang2/R_scripts/AlphaSimR_test/dev/testData/msprime_chr0.trees",
       #breaks=c(0, L1/2, L1), rates=c(1e-8, 2e-8), segSites=60),
  breaks=c(0, L1/2, L1), rates=c(1e-5, 2e-5), segSites=60),
  list(ts_path="/Users/jliang2/R_scripts/AlphaSimR_test/dev/testData/msprime_chr1.trees",
       #breaks=c(0, L2/3, 2*L2/3, L2), rates=c(1e-7, 1e-8, 1e-7), segSites=155)
  breaks=c(0, L2/3, 2*L2/3, L2), rates=c(1e-4, 1e-5, 1e-4), segSites=155)
)

founderGenomes2 <- asMapPop(chr_info = chr_info, inbred=FALSE, ploidy=2L)
set.seed(42)
SP2 = SimParam$new(founderGenomes2)
SP2$setSexes("yes_sys")
SP2$addTraitA(nQtlPerChr = 5,
             mean = 500,
             var = 450)

SP2$setTrackPed(TRUE)
# try the new function here, it automatically set setTrackRec also.
SP2$setTrackRecGen(TRUE)
basePop2 = newPop(founderGenomes2, simParam = SP2)
basePop2 = setPheno(basePop2,
                   h2 = 0.5,
                   simParam = SP2)

#--- n generations
nCycles<-2

# very simple container for each cycles sim output
simOutput2<-list(basePop2)
cycle<-1
for(cycle in 1:nCycles){
  cat(paste0(" C",cycle))
  # choose the best from last cycle
  chosenParents<- selectInd(pop=simOutput2[[cycle]], nInd=6, use = "gv", simParam = SP2)
  # make crosses
  offspringPop<-randCross(pop=chosenParents, nCrosses=2, nProgeny = 5, simParam = SP2)
  # phenotype  new offspring
  offspringPop<-setPheno(pop = offspringPop, h2 = 0.5, simParam = SP2)
  # add new offspring to simOutput list
  simOutput2[[cycle+1]]<-offspringPop
}

RHG <- SP2$recHistGen
RH <- SP2$recHist
gm <- SP2$genMap[[1]]

rh <- RH[[3]][[1]][[1]]
rhg <- RHG[[3]][[1]][[1]]
rh
rhg
x  <- rhg[,2]

left  <- findInterval(x, gm)
right <- pmin(left + 1, length(gm))

out <- data.frame(
  x = x,
  left_i  = left,
  left_v  = gm[left],
  right_i = right,
  right_v = gm[right]
)

out



