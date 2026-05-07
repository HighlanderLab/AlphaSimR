```
rm(list = ls())
```

0. load dependencies (already in the R scripts, but if you have different setting plz just do this in your way)
```
library(reticulate)
library(AlphaSimR)

use_virtualenv("~/r-reticulate-env", required = TRUE)
msprime <- import("msprime")
tskit <- import("tskit")
```
1. load functions
```
devtools::load_all()
#--- or ----
source("R/makeFoundersFromTs.R")
source("R/alphaSimR2Ts.R")
```

2. read the .trees files with 2 chromosomes and 2 dip individuals (from msprime)
note: chrKeptPosBpList added, so we know the position and index of sampled SNPs in the original .trees files (alphaSimR only record their index)
```
L1 <- 1e6
L2 <- 2e6

chr_info <- list(
  list(ts_path="dev/testData/msprime_chr0.trees",
       breaks=c(0, L1/2, L1), rates=c(1e-8, 2e-8), segSites=60),
  list(ts_path="dev/testData/msprime_chr1.trees",
       breaks=c(0, L2/3, 2*L2/3, L2), rates=c(1e-7, 1e-8, 1e-7), segSites=155)
)

founderGenomes <- asMapPop(chr_info = chr_info, inbred=FALSE, ploidy=2L)
```

3. run alphaSimR to set the founder genomes and parameters
```
set.seed(42)
SP = SimParam$new(founderGenomes)
SP$setSexes("yes_sys")
SP$addTraitA(nQtlPerChr = 5,
             mean = 500,
             var = 450)
SP$setTrackPed(TRUE)
SP$setTrackRec(TRUE)
basePop = newPop(founderGenomes)
basePop = setPheno(basePop,
                   h2 = 0.5)
```                   

4. run addition 2 generations
```
#--- n generations
nCycles<-2

# keep founderPop and offspringPop in SimOutput
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
```

5. Link recHist (from SP, based on sampled SNPs) with parent-child hap (in tskit positions, from chrKeptPosBpList)
```
bridgeCollectSegFromSimOutput(SP, simOutput)
```

6. make an edge table from bridgeSegDfList
```
edgeDf <- bridgeAllSegToEdgeDf(chr_info)
```

7. write tskit tables (nodes and edges)
note1: time of founder generation: 0; time of offspring: time of the youngest parent - 1
note2: check nodeIdMapByChr for ids of alphaSimR and tskit
note3: n ploidy is from asMapPop, so variable number of ploidy along generations is not allowed
note4: be careful with the metadata in the future (different behaviors between Python and R even with Reticulate)
```
bridgeWriteTrees(chr_info, edgeDf, SP)
```

8. We can see that there is no new recombination break points in chr1, let's play with chr2

```
# in Python:
import tskit
ts0 = tskit.load('.../dev/testData/msprime_chr1.trees')
ts1 = tskit.load('.../dev/testData/_AlphaSimR_extended_chr1.trees')
```

From edgeDf, there's a breakpoint at 1549443
```
ts0_1549443 = ts0.at(1549442)
ts0_1549444 = ts0.at(1549443)
ts1_1549443 = ts1.at(1549442)
ts1_1549444 = ts1.at(1549443)
```
Same tree in the original file:
```
print(ts0_1549443.draw_text())
```
```
Output:
   182 
  ┏━┻━┓
 32   ┃
 ┏┻━┓ ┃
12 1621
┏┻┓ ┃ ┃
0 2 1 3
```

```
print(ts0_1549444.draw_text())
```
```
Output:
   182 
  ┏━┻━┓
 32   ┃
 ┏┻━┓ ┃
12 1621
┏┻┓ ┃ ┃
0 2 1 3
```
Different trees in the new file
```
print(ts1_1549443.draw_text())
```
```
Output:
                                       182                                                                                         
              ┏━━━━━━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━┓                                                                 
              ┃                                                 32                                                                 
              ┃                            ┏━━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━┓                                           
             21                           16                                          12                                           
              ┃                            ┃                           ┏━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━┓                           
              3                            1                           2                               0                           
 ┏━━━┳━━━┳━━━┳┻━━━━┳━━━━━━━┓     ┏━━━┳━━━┳━┻━━━━━━━━━┓           ┏━━━┳━┻━┳━━━┓   ┏━━━┳━━━┳━━━━━━━━━┳━━━┻━━━━━━━━━┳━━━━━━━━━━━┓     
491 493 505 497   499     501   494 504 508         492         495 503 507 509 496 506 510       498           500         502    
             ┃   ┏━┻━┓   ┏━┻━┓               ┏━━━┳━━━╋━━━┳━━━┓                               ┏━━━┳━┻━┳━━━┓   ┏━━━╋━━━┓   ┏━━━╋━━━┓ 
            511 512 520 523 529             522 524 526 528 530                             513 515 517 519 514 516 518 521 525 527
```
```
print(ts1_1549444.draw_text())
```
```
Output:
                                                                                    182                                            
                                                      ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓              
                                                     32                                                             ┃              
                          ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓                                ┃              
                         12                                                       16                               21              
         ┏━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━┓                                       ┃                                ┃              
         2                                 0                                       1                                3              
 ┏━━━┳━━━╋━━━┳━━━┓   ┏━━━┳━━━┳━━━━━━━━━┳━━━┻━━━━━━━━━┳━━━━━━━━━━━┓       ┏━━━┳━━━┳━┻━━━━━━━━━┓           ┏━━━┳━━━┳━━┻━━┳━━━━━━━┓   
491 495 503 507 509 496 506 510       498           500         502     494 504 508         492         493 505 497   499     501  
                                 ┏━━━┳━┻━┳━━━┓   ┏━━━╋━━━┓   ┏━━━╋━━━┓               ┏━━━┳━━━╋━━━┳━━━┓           ┃   ┏━┻━┓   ┏━┻━┓ 
                                513 515 517 519 514 516 518 521 525 527             522 524 526 528 530         511 512 520 523 529

```
difference: one of the parent nodes of node 491 (3_1 in alphaSimR) changed from 3 (2_2) to 2 (2_1), the same as edgeDf.


New nodes look like:
```
ts1.tables.nodes[491]
```
```
Output:
NodeTableRow(flags=0, time=-1.0, population=-1, individual=-1, metadata={'alphaSimR': {'id': '3_1'}})
```
Founder nodes look like:
```
ts1.tables.nodes[0]
```
```
Output:
NodeTableRow(flags=1, time=0.0, population=0, individual=0, metadata={'alphaSimR': {'id': '1_1'}})
```
