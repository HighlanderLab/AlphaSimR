#' @title Import genetic map
#' 
#' @description
#' Formats a genetic map stored in a data.frame to
#' AlphaSimR's internal format. Map positions must be 
#' in Morgans. 
#' 
#' @param genMap genetic map as a data.frame. The first  
#' three columns must be: marker name, chromosome, and 
#' map position (Morgans). Marker name and chromosome are 
#' coerced using as.character.
#'
#' @return a list of named vectors
#' 
#' @examples 
#' genMap = data.frame(markerName=letters[1:5],
#'                     chromosome=c(1,1,1,2,2),
#'                     position=c(0,0.5,1,0.15,0.4))
#' 
#' asrMap = importGenMap(genMap=genMap)
#' 
#' str(asrMap)
#' 
#' @export
importGenMap = function(genMap){
  # Convert data type
  markerName = as.character(genMap[,1])
  chromosome = as.character(genMap[,2])
  position = as.numeric(genMap[,3])

  # Create list for map
  uniqueChr = unique(chromosome)
  genMap = vector("list", length=length(uniqueChr))
  names(genMap) = uniqueChr

  # Iterate through chromosomes
  for(i in seq_len(length(uniqueChr))){

    take = (chromosome==uniqueChr[i])
    tmpPos = position[take]
    tmpName = markerName[take]

    # Order and name
    take = order(tmpPos, decreasing=FALSE)
    tmpPos = tmpPos[take]
    names(tmpPos) = tmpName[take]

    genMap[[uniqueChr[i] ]] = tmpPos - tmpPos[1]
  }

  return(genMap)
}
 
#' @title Import inbred, diploid genotypes
#' 
#' @description
#' Formats the genotypes from inbred, diploid lines 
#' to an AlphaSimR population that can be used to 
#' initialize a simulation. An attempt is made to 
#' automatically detect 0,1,2 or -1,0,1 genotype coding. 
#' Heterozygotes or probabilistic genotypes are allowed, 
#' but will be coerced to the nearest homozygote. Pedigree 
#' information is optional and when provided will be 
#' passed to the population for easier identification 
#' in the simulation.
#' 
#' @param geno a matrix of genotypes
#' @param genMap genetic map as a data.frame. The first  
#' three columns must be: marker name, chromosome, and 
#' map position (Morgans). Marker name and chromosome are 
#' coerced using as.character. See \link{importGenMap}
#' @param ped an optional pedigree for the supplied 
#' genotypes. See details. 
#'
#' @details 
#' The optional pedigree can be a data.frame, matrix or a vector. 
#' If the object is a data.frame or matrix, the first three 
#' columns must include information in the following order: id, 
#' mother, and father. All values are coerced using 
#' as.character. If the object is a vector, it is assumed to only 
#' include the id. In this case, the mother and father will be set 
#' to "0" for all individuals.
#' 
#' @return a \code{\link{MapPop-class}} if ped is NULL,
#' otherwise a \code{\link{NamedMapPop-class}}
#' 
#' @examples 
#' geno = rbind(c(2,2,0,2,0),
#'              c(0,2,2,0,0))
#' colnames(geno) = letters[1:5]
#' 
#' genMap = data.frame(markerName=letters[1:5],
#'                     chromosome=c(1,1,1,2,2),
#'                     position=c(0,0.5,1,0.15,0.4))
#' 
#' ped = data.frame(id=c("a","b"),
#'                  mother=c(0,0),
#'                  father=c(0,0))
#' 
#' founderPop = importInbredGeno(geno=geno,
#'                               genMap=genMap,
#'                               ped=ped)
#' 
#' @export
importInbredGeno = function(geno, genMap, ped=NULL){
  # Extract pedigree, if supplied
  if(!is.null(ped)){
    if(is.vector(ped)){
      id = as.character(ped)
      stopifnot(length(id)==nrow(geno),
                !any(duplicated(id)))
      mother = father = rep("0", length(id))
    }else{
      id = as.character(ped[,1])
      stopifnot(length(id)==nrow(geno),
                !any(duplicated(id)))
      mother = as.character(ped[,2])
      father = as.character(ped[,3])
    }
  }
  
  genMap = importGenMap(genMap)
  
  # Get marker names
  if(is.data.frame(geno)){
    geno = as.matrix(geno)
  }
  markerName = colnames(geno)

  # Check marker coding and convert to haplotypes
  if(is.raw(geno)){
    geno[geno==as.raw(1)] = as.raw(0) # For consistency with round
    geno[geno==as.raw(2)] = as.raw(1)
  }else{
    minGeno = min(geno)
    maxGeno = max(geno)
    stopifnot(minGeno >= (-1-1e-8) )
    
    if(minGeno < (0-1e-8) ){
      # Suspect -1,0,1 coding
      stopifnot(maxGeno <= (1+1e-8) )
      # Converting to 0,1 haplotypes with hard thresholds
      geno = matrix(as.raw( round( (geno+1)/2 ) ),
                    ncol=ncol(geno))
    }else{
      # Suspect 0,1,2 coding
      stopifnot(maxGeno <= (2+1e-8) )
      # Converting to 0,1 haplotypes with hard thresholds
      geno = matrix(as.raw( round( geno/2 ) ),
                    ncol=ncol(geno))
    }
  }

  # Create haplotype list
  haplotypes = vector("list", length=length(genMap))

  # Order haplotypes by chromosome
  for(i in seq_len(length(genMap))){
    mapMarkers = names(genMap[[i]])
    take = match(mapMarkers, markerName)
    if(any(is.na(take))){
      genMap[[i]] = genMap[[i]][is.na(take)]
      stopifnot(length(genMap[[i]]) >= 1L)
      genMap[[i]] = genMap[[i]] - genMap[[i]]-genMap[[i]][1]
      take = na.omit(take)
    }
    haplotypes[[i]] = geno[,take]
  }

  founderPop = newMapPop(genMap=genMap,
                         haplotypes=haplotypes,
                         inbred=TRUE)

  if(!is.null(ped)){
    founderPop = new("NamedMapPop",
                     id=id,
                     mother=mother,
                     father=father,
                     founderPop)
  }

  return(founderPop)
}


#' @title Import haplotypes
#' 
#' @description
#' Formats haplotype in a matrix format to an 
#' AlphaSimR population that can be used to 
#' initialize a simulation. This function serves 
#' as wrapper for \code{\link{newMapPop}} that 
#' utilizes a more user friendly input format.
#' 
#' @param haplo a matrix of haplotypes
#' @param genMap genetic map as a data.frame. The first  
#' three columns must be: marker name, chromosome, and 
#' map position (Morgans). Marker name and chromosome are 
#' coerced using as.character. See \code{\link{importGenMap}}
#' @param ploidy ploidy level of the organism
#' @param ped an optional pedigree for the supplied 
#' genotypes. See details. 
#' 
#' @details 
#' The optional pedigree can be a data.frame, matrix or a vector. 
#' If the object is a data.frame or matrix, the first three 
#' columns must include information in the following order: id, 
#' mother, and father. All values are coerced using 
#' as.character. If the object is a vector, it is assumed to only 
#' include the id. In this case, the mother and father will be set 
#' to "0" for all individuals.
#'
#' @return a \code{\link{MapPop-class}} if ped is NULL,
#' otherwise a \code{\link{NamedMapPop-class}}
#' 
#' @examples 
#' haplo = rbind(c(1,1,0,1,0),
#'               c(1,1,0,1,0),
#'               c(0,1,1,0,0),
#'               c(0,1,1,0,0))
#' colnames(haplo) = letters[1:5]
#' 
#' genMap = data.frame(markerName=letters[1:5],
#'                     chromosome=c(1,1,1,2,2),
#'                     position=c(0,0.5,1,0.15,0.4))
#' 
#' ped = data.frame(id=c("a","b"),
#'                  mother=c(0,0),
#'                  father=c(0,0))
#' 
#' founderPop = importHaplo(haplo=haplo, 
#'                          genMap=genMap,
#'                          ploidy=2L,
#'                          ped=ped)
#' 
#' @export
importHaplo = function(haplo, genMap, ploidy=2L, ped=NULL){
  # Extract pedigree, if supplied
  if(!is.null(ped)){
    if(is.vector(ped)){
      id = as.character(ped)
      stopifnot(length(id)==(nrow(haplo)/ploidy),
                !any(duplicated(id)))
      mother = father = rep("0", length(id))
    }else{
      id = as.character(ped[,1])
      stopifnot(length(id)==(nrow(haplo)/ploidy),
                !any(duplicated(id)))
      mother = as.character(ped[,2])
      father = as.character(ped[,3])
    }
  }
  
  genMap = importGenMap(genMap)
  
  # Get marker names
  if(is.data.frame(haplo)){
    haplo = as.matrix(haplo)
  }
  markerName = colnames(haplo)
  
  # Convert haplotypes to raw
  haplo = matrix(as.raw(haplo), ncol=ncol(haplo))
  stopifnot(haplo==as.raw(0) | haplo==as.raw(1))
  
  # Create haplotype list
  haplotypes = vector("list", length=length(genMap))
  
  # Order haplotypes by chromosome
  for(i in seq_len(length(genMap))){
    mapMarkers = names(genMap[[i]])
    take = match(mapMarkers, markerName)
    if(any(is.na(take))){
      genMap[[i]] = genMap[[i]][is.na(take)]
      stopifnot(length(genMap[[i]]) >= 1L)
      genMap[[i]] = genMap[[i]] - genMap[[i]]-genMap[[i]][1]
      take = na.omit(take)
    }
    haplotypes[[i]] = haplo[,take,drop=FALSE]
  }
  
  founderPop = newMapPop(genMap=genMap,
                         haplotypes=haplotypes,
                         ploidy=ploidy)
  
  if(!is.null(ped)){
    founderPop = new("NamedMapPop",
                     id=id,
                     mother=mother,
                     father=father,
                     founderPop)
  }
  
  return(founderPop)
}

#' @title Import VCF haplotypes
#'
#' @description
#' Streams a VCF file and imports phased genotype calls as founder
#' haplotypes. Sites are filtered to biallelic segregating sites, missing
#' genotype calls are filtered or rejected, and optional reservoir sampling is
#' applied independently within each chromosome.
#'
#' @param vcfFile path to a VCF file. Plain text and gzip-compressed files are
#' supported.
#' @param breaks recombination map breakpoints. Supply either a numeric vector
#' used for every chromosome, or a list of numeric vectors with one entry per
#' chromosome. Named lists are matched to VCF chromosome names.
#' @param rates recombination rates for the intervals defined by \code{breaks}.
#' Supply either a numeric vector used for every chromosome, or a list of
#' numeric vectors with one entry per chromosome. Each rates vector must have
#' length \code{length(breaks) - 1}.
#' @param segSites optional number of segregating sites to keep per chromosome.
#' If \code{NULL}, all qualifying sites are retained. A scalar value is used for
#' every chromosome. A vector or list can be supplied per chromosome; named
#' values are matched to VCF chromosome names.
#' @param siteSamplingSeed integer seed used for reservoir sampling.
#' @param ploidy optional ploidy level. If \code{NULL}, ploidy is inferred from
#' the VCF genotype calls and checked for consistency.
#' @param ped an optional pedigree for the supplied genotypes. See details.
#' @param missing how to handle missing genotype calls. Use \code{"filter"} to
#' drop sites with missing calls or \code{"error"} to stop.
#' @param requirePhased if \code{TRUE}, genotype calls with \code{/} separators
#' are rejected.
#' @param useVCFIds if \code{TRUE} and \code{ped} is \code{NULL}, VCF sample IDs
#' are used to return a \code{\link{NamedMapPop-class}}.
#' @param tsRecord if \code{TRUE}, initialize in-memory founder tree-sequence
#' tables from the imported VCF haplotypes and attach metadata for
#' \code{SimParam$setTrackTs(TRUE, founderPop=...)}.
#' @param tsRecorde deprecated alias for \code{tsRecord}.
#' @param addTsMut if \code{TRUE} and \code{tsRecord=TRUE}, add retained VCF
#' sites and allele-1 calls as synthetic tskit site/mutation rows. If
#' \code{FALSE}, initialize only founder individuals and sample nodes, leaving
#' the TS site and mutation tables empty.
#' @param seqLen optional chromosome sequence lengths for \code{tsRecord=TRUE}.
#' Supply either a scalar, vector, or list. Named values are matched to VCF
#' chromosome names. If \code{NULL}, the maximum breakpoint for each chromosome
#' is used.
#' @param returnMeta if \code{TRUE}, return a list containing the population,
#' retained physical positions, sample IDs, and scan statistics.
#'
#' @details
#' The optional pedigree follows the same format as \code{\link{importHaplo}}.
#' If the object is a data.frame or matrix, the first three columns must include
#' id, mother, and father. If the object is a vector, it is assumed to only
#' include the id. In this case, the mother and father will be set to \code{"0"}
#' for all individuals.
#'
#' A site is retained only when it has a single ALT allele, every non-missing
#' genotype allele is coded as 0 or 1, and both alleles are observed in the
#' founder haplotypes.
#'
#' When \code{tsRecord=TRUE} and \code{addTsMut=TRUE}, duplicate retained
#' physical positions within a chromosome are not allowed because tskit requires
#' strictly increasing site positions.
#'
#' @return a \code{\link{MapPop-class}} or \code{\link{NamedMapPop-class}}. If
#' \code{returnMeta = TRUE}, a list with elements \code{pop}, \code{keptPos},
#' \code{sampleIds}, \code{stats}, and, when \code{tsRecord=TRUE},
#' \code{tsTables}.
#'
#' @examples
#' \dontrun{
#' founderPop = importVCF("founders.vcf.gz",
#'                        breaks=list(c(0, 1e8)),
#'                        rates=list(c(1e-8)),
#'                        segSites=1000)
#'
#' # Initialize founder tree-sequence tables for forward recording.
#' founderPopTs = importVCF("founders.vcf.gz",
#'                         breaks=list(c(0, 1e8)),
#'                         rates=list(c(1e-8)),
#'                         segSites=1000,
#'                         tsRecord=TRUE,
#'                         addTsMut=FALSE)
#' SP = SimParam$new(founderPopTs)
#' SP$setTrackTs(TRUE, founderPop=founderPopTs)
#' }
#'
#' @export
importVCF = function(vcfFile, breaks, rates, segSites=NULL,
                     siteSamplingSeed=42L, ploidy=NULL, ped=NULL,
                     missing=c("filter", "error"), requirePhased=TRUE,
                     useVCFIds=TRUE, tsRecord=FALSE, tsRecorde=NULL,
                     addTsMut=TRUE, seqLen=NULL,
                     returnMeta=FALSE){
  missing = match.arg(missing)
  requirePhased = isTRUE(requirePhased)
  useVCFIds = isTRUE(useVCFIds)
  if(!is.null(tsRecorde)){
    tsRecord = tsRecorde
  }
  tsRecord = isTRUE(tsRecord)
  addTsMut = isTRUE(addTsMut)
  returnMeta = isTRUE(returnMeta)
  if(!is.null(ploidy)){
    ploidy = as.integer(ploidy)
    if(length(ploidy)!=1L || is.na(ploidy) || ploidy<1L){
      stop("`ploidy` must be NULL or a single positive integer.",
           call.=FALSE)
    }
  }

  vcfData = .sampleVcfVariants(vcfFile=vcfFile,
                               segSites=segSites,
                               seed=siteSamplingSeed,
                               ploidy=ploidy,
                               requirePhased=requirePhased,
                               missing=missing)

  mapData = .vcfBuildMapData(chrData=vcfData$chrData,
                             breaks=breaks,
                             rates=rates,
                             seqLen=seqLen)
  genMap = mapData$genMap
  haplotypes = lapply(vcfData$chrData, `[[`, "H")

  founderPop = newMapPop(genMap=genMap,
                         haplotypes=haplotypes,
                         ploidy=vcfData$ploidy,
                         inbred=FALSE)

  founderPop = .vcfAddPed(founderPop=founderPop,
                          ped=ped,
                          sampleIds=vcfData$sampleIds,
                          useVCFIds=useVCFIds)

  tsTables = NULL
  if(tsRecord){
    tsTables = .vcfBuildFounderTables(chrData=vcfData$chrData,
                                      seqLenList=mapData$seqLenList,
                                      ploidy=vcfData$ploidy,
                                      addTsMut=addTsMut)
    attr(founderPop, "tsForwardSource") = lapply(tsTables, function(x){
      list(tc_xptr=x)
    })
    attr(founderPop, "tsForwardPosMeta") = list(
      posList=lapply(vcfData$chrData, `[[`, "P"),
      seqLenList=mapData$seqLenList,
      breaksList=mapData$breaksList,
      ratesList=mapData$ratesList
    )
  }

  if(returnMeta){
    out = list(pop=founderPop,
               keptPos=lapply(vcfData$chrData, `[[`, "P"),
               sampleIds=vcfData$sampleIds,
               stats=vcfData$stats)
    if(tsRecord){
      out$tsTables = tsTables
    }
    return(out)
  }
  return(founderPop)
}

#' @keywords internal
#' @noRd
.sampleVcfVariants = function(vcfFile, segSites, seed, ploidy,
                              requirePhased, missing){
  if(!file.exists(vcfFile)){
    stop("VCF file does not exist: ", vcfFile, call.=FALSE)
  }

  oldSeedExists = exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
  if(oldSeedExists){
    oldSeed = get(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
  }
  on.exit({
    if(oldSeedExists){
      assign(".Random.seed", oldSeed, envir=.GlobalEnv)
    }else if(exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)){
      rm(".Random.seed", envir=.GlobalEnv)
    }
  }, add=TRUE)
  set.seed(as.integer(seed))

  con = .vcfOpen(vcfFile)
  on.exit(close(con), add=TRUE)

  sampleIds = NULL
  chrData = list()
  chrOrder = character()
  currentPloidy = ploidy
  stats = list(records=0L,
               qualifying=0L,
               skippedNonBiallelic=0L,
               skippedMissing=0L,
               skippedInvalidAllele=0L,
               skippedNonSegregating=0L)

  repeat{
    lines = readLines(con, n=10000L, warn=FALSE)
    if(length(lines)==0L){
      break
    }
    for(line in lines){
      if(startsWith(line, "##")){
        next
      }
      if(startsWith(line, "#CHROM")){
        header = strsplit(line, "\t", fixed=TRUE)[[1]]
        if(length(header)<10L){
          stop("VCF header must contain sample genotype columns.",
               call.=FALSE)
        }
        sampleIds = header[-seq_len(9L)]
        if(any(sampleIds=="") || any(duplicated(sampleIds))){
          stop("VCF sample IDs must be non-empty and unique.",
               call.=FALSE)
        }
        next
      }
      if(is.null(sampleIds)){
        stop("VCF header line beginning with #CHROM was not found.",
             call.=FALSE)
      }
      if(startsWith(line, "#")){
        next
      }

      fields = strsplit(line, "\t", fixed=TRUE)[[1]]
      if(length(fields)<(9L+length(sampleIds))){
        stop("Malformed VCF record with too few columns.", call.=FALSE)
      }
      stats$records = stats$records + 1L

      chr = as.character(fields[[1L]])
      pos = as.numeric(fields[[2L]])
      id = fields[[3L]]
      alt = fields[[5L]]

      if(is.na(pos) || alt=="." || grepl(",", alt, fixed=TRUE)){
        stats$skippedNonBiallelic = stats$skippedNonBiallelic + 1L
        next
      }

      gt = .vcfParseGt(sampleFields=fields[-seq_len(9L)],
                       formatField=fields[[9L]],
                       sampleIds=sampleIds,
                       ploidy=currentPloidy,
                       requirePhased=requirePhased,
                       missing=missing,
                       chr=chr,
                       pos=pos)

      if(gt$status=="missing"){
        stats$skippedMissing = stats$skippedMissing + 1L
        next
      }
      if(gt$status=="invalidAllele"){
        stats$skippedInvalidAllele = stats$skippedInvalidAllele + 1L
        next
      }
      if(is.null(currentPloidy)){
        currentPloidy = gt$ploidy
      }
      if(length(unique(gt$haplo))<2L){
        stats$skippedNonSegregating = stats$skippedNonSegregating + 1L
        next
      }

      if(!chr%in%chrOrder){
        chrOrder = c(chrOrder, chr)
        target = .vcfSegSitesForChr(segSites=segSites,
                                    chr=chr,
                                    chrIndex=length(chrOrder))
        chrData[[chr]] = .vcfNewReservoir(nHaplo=length(gt$haplo),
                                          target=target)
      }

      markerName = if(id=="." || id==""){
        paste(chr, format(pos, scientific=FALSE, trim=TRUE), sep="_")
      }else{
        id
      }
      chrData[[chr]] = .vcfReservoirAdd(res=chrData[[chr]],
                                        haplo=gt$haplo,
                                        pos=pos,
                                        markerName=markerName)
      stats$qualifying = stats$qualifying + 1L
    }
  }

  if(is.null(sampleIds)){
    stop("VCF header line beginning with #CHROM was not found.",
         call.=FALSE)
  }
  if(length(chrData)==0L){
    stop("No qualifying biallelic segregating sites were found.",
         call.=FALSE)
  }
  chrData = chrData[chrOrder]
  chrData = lapply(names(chrData), function(chr){
    .vcfFinalizeReservoir(res=chrData[[chr]], chr=chr)
  })
  names(chrData) = chrOrder

  stats$kept = sum(vapply(chrData, function(x) ncol(x$H), integer(1)))
  stats$keptByChr = vapply(chrData, function(x) ncol(x$H), integer(1))

  return(list(chrData=chrData,
              sampleIds=sampleIds,
              ploidy=currentPloidy,
              stats=stats))
}

#' @keywords internal
#' @noRd
.vcfOpen = function(vcfFile){
  if(grepl("\\.gz$", vcfFile, ignore.case=TRUE)){
    con = gzfile(vcfFile, open="rt")
  }else{
    con = file(vcfFile, open="rt")
  }
  return(con)
}

#' @keywords internal
#' @noRd
.vcfParseGt = function(sampleFields, formatField, sampleIds, ploidy,
                       requirePhased, missing, chr, pos){
  formatFields = strsplit(formatField, ":", fixed=TRUE)[[1]]
  gtIndex = match("GT", formatFields)
  if(is.na(gtIndex)){
    stop("VCF FORMAT field does not contain GT.", call.=FALSE)
  }

  gt = vapply(sampleFields, function(x){
    value = strsplit(x, ":", fixed=TRUE)[[1]]
    if(length(value)<gtIndex){
      return(NA_character_)
    }
    return(value[[gtIndex]])
  }, character(1))

  if(any(is.na(gt) | gt=="") || any(grepl("\\.", gt))){
    if(missing=="error"){
      stop("Missing genotype call at ", chr, ":", pos, ".", call.=FALSE)
    }
    return(list(status="missing"))
  }

  if(requirePhased && any(grepl("/", gt, fixed=TRUE))){
    stop("VCF contains unphased genotype call at ", chr, ":", pos,
         ". Phased GT values with `|` are required.", call.=FALSE)
  }

  if(requirePhased){
    alleles = strsplit(gt, "|", fixed=TRUE)
  }else{
    alleles = strsplit(gt, "[|/]")
  }
  alleleCounts = lengths(alleles)
  if(length(unique(alleleCounts))!=1L){
    stop("Inconsistent ploidy among genotype calls at ", chr, ":", pos,
         ".", call.=FALSE)
  }
  inferredPloidy = alleleCounts[[1L]]
  if(!is.null(ploidy) && inferredPloidy!=ploidy){
    stop("VCF genotype ploidy at ", chr, ":", pos,
         " does not match expected ploidy ", ploidy, ".", call.=FALSE)
  }
  if(inferredPloidy<1L){
    stop("Invalid genotype ploidy at ", chr, ":", pos, ".", call.=FALSE)
  }

  alleles = unlist(alleles, use.names=FALSE)
  if(any(!alleles%in%c("0", "1"))){
    return(list(status="invalidAllele"))
  }

  haplo = as.integer(alleles)
  names(haplo) = paste(rep(sampleIds, each=inferredPloidy),
                       rep(seq_len(inferredPloidy), length(sampleIds)),
                       sep="_")

  return(list(status="ok",
              haplo=haplo,
              ploidy=inferredPloidy))
}

#' @keywords internal
#' @noRd
.vcfSegSitesForChr = function(segSites, chr, chrIndex){
  if(is.null(segSites)){
    return(NULL)
  }
  value = .vcfComponentForChr(x=segSites,
                              chr=chr,
                              chrIndex=chrIndex,
                              name="segSites",
                              vectorIsPerChr=TRUE)
  value = as.integer(value[[1L]])
  if(is.na(value) || value<=0L){
    return(NULL)
  }
  return(value)
}

#' @keywords internal
#' @noRd
.vcfNewReservoir = function(nHaplo, target){
  if(is.null(target)){
    H = matrix(integer(), nrow=nHaplo, ncol=0L)
    P = numeric(0L)
    markerNames = character(0L)
  }else{
    H = matrix(NA_integer_, nrow=nHaplo, ncol=target)
    P = rep(NA_real_, target)
    markerNames = rep(NA_character_, target)
  }
  return(list(H=H,
              P=P,
              markerNames=markerNames,
              target=target,
              seen=0L,
              size=0L))
}

#' @keywords internal
#' @noRd
.vcfReservoirAdd = function(res, haplo, pos, markerName){
  res$seen = res$seen + 1L
  if(is.null(res$target)){
    res$size = res$size + 1L
    res$H = cbind(res$H, haplo)
    res$P = c(res$P, pos)
    res$markerNames = c(res$markerNames, markerName)
  }else if(res$size < res$target){
    res$size = res$size + 1L
    res$H[,res$size] = haplo
    res$P[[res$size]] = pos
    res$markerNames[[res$size]] = markerName
  }else{
    j = sample.int(res$seen, 1L)
    if(j <= res$target){
      res$H[,j] = haplo
      res$P[[j]] = pos
      res$markerNames[[j]] = markerName
    }
  }
  return(res)
}

#' @keywords internal
#' @noRd
.vcfFinalizeReservoir = function(res, chr){
  if(!is.null(res$target) && res$size < res$target){
    stop("Insufficient qualifying sites on chromosome ", chr,
         " (found ", res$size, ", requested ", res$target, ").",
         call.=FALSE)
  }
  if(res$size==0L){
    stop("No qualifying sites on chromosome ", chr, ".", call.=FALSE)
  }
  take = seq_len(res$size)
  H = res$H[,take,drop=FALSE]
  P = res$P[take]
  markerNames = res$markerNames[take]

  ord = order(P)
  H = H[,ord,drop=FALSE]
  P = P[ord]
  markerNames = make.unique(markerNames[ord])

  return(list(H=H,
              P=P,
              markerNames=markerNames,
              seen=res$seen))
}

#' @keywords internal
#' @noRd
.vcfBuildMapData = function(chrData, breaks, rates, seqLen=NULL){
  chrNames = names(chrData)
  genMap = vector("list", length(chrData))
  breaksList = vector("list", length(chrData))
  ratesList = vector("list", length(chrData))
  seqLenList = vector("list", length(chrData))
  names(genMap) = chrNames
  names(breaksList) = chrNames
  names(ratesList) = chrNames
  names(seqLenList) = chrNames
  for(i in seq_along(chrData)){
    chr = chrNames[[i]]
    chrBreaks = .vcfMapComponentForChr(breaks, chr, i, "breaks")
    chrRates = .vcfMapComponentForChr(rates, chr, i, "rates")
    if(length(chrBreaks)!=(length(chrRates)+1L)){
      stop("`rates` for chromosome ", chr,
           " must have length length(breaks) - 1.", call.=FALSE)
    }
    if(any(!is.finite(chrBreaks)) || any(!is.finite(chrRates))){
      stop("`breaks` and `rates` must be finite for chromosome ", chr,
           ".", call.=FALSE)
    }
    if(any(diff(chrBreaks)<0)){
      stop("`breaks` must be sorted for chromosome ", chr, ".",
           call.=FALSE)
    }
    if(any(chrRates<0)){
      stop("`rates` must be non-negative for chromosome ", chr, ".",
           call.=FALSE)
    }
    mpos = rateMap2cumMorgan(chrData[[i]]$P, chrBreaks, chrRates)
    mpos = mpos - mpos[[1L]]
    names(mpos) = chrData[[i]]$markerNames
    genMap[[i]] = mpos
    breaksList[[i]] = chrBreaks
    ratesList[[i]] = chrRates
    seqLenList[[i]] = .vcfSeqLenForChr(seqLen=seqLen,
                                       chr=chr,
                                       chrIndex=i,
                                       breaks=chrBreaks)
  }
  return(list(genMap=genMap,
              breaksList=breaksList,
              ratesList=ratesList,
              seqLenList=seqLenList))
}

#' @keywords internal
#' @noRd
.vcfMapComponentForChr = function(x, chr, chrIndex, name){
  if(missing(x) || is.null(x)){
    stop("`", name, "` is required.", call.=FALSE)
  }
  value = .vcfComponentForChr(x=x,
                              chr=chr,
                              chrIndex=chrIndex,
                              name=name,
                              vectorIsPerChr=FALSE)
  value = as.numeric(value)
  if(length(value)==0L || any(is.na(value))){
    stop("`", name, "` contains missing values for chromosome ", chr,
         ".", call.=FALSE)
  }
  return(value)
}

#' @keywords internal
#' @noRd
.vcfSeqLenForChr = function(seqLen, chr, chrIndex, breaks){
  if(is.null(seqLen)){
    value = max(breaks)
  }else{
    value = .vcfComponentForChr(x=seqLen,
                                chr=chr,
                                chrIndex=chrIndex,
                                name="seqLen",
                                vectorIsPerChr=TRUE)
  }
  value = as.numeric(value[[1L]])
  if(length(value)!=1L || is.na(value) || !is.finite(value) || value<=0){
    stop("`seqLen` must contain finite positive values.", call.=FALSE)
  }
  return(value)
}

#' @keywords internal
#' @noRd
.vcfComponentForChr = function(x, chr, chrIndex, name, vectorIsPerChr){
  if(is.list(x)){
    if(!is.null(names(x)) && chr%in%names(x)){
      return(x[[chr]])
    }
    if(length(x)==1L){
      return(x[[1L]])
    }
    if(chrIndex<=length(x)){
      return(x[[chrIndex]])
    }
    stop("`", name, "` must have length 1 or nChr.", call.=FALSE)
  }

  if(!vectorIsPerChr){
    return(x)
  }
  if(!is.null(names(x)) && chr%in%names(x)){
    return(x[[chr]])
  }
  if(length(x)==1L){
    return(x[[1L]])
  }
  if(chrIndex<=length(x)){
    return(x[[chrIndex]])
  }
  stop("`", name, "` must have length 1 or nChr.", call.=FALSE)
}

#' @keywords internal
#' @noRd
.vcfBuildFounderTables = function(chrData, seqLenList, ploidy, addTsMut){
  out = vector("list", length(chrData))
  names(out) = names(chrData)
  for(i in seq_along(chrData)){
    pos = as.numeric(chrData[[i]]$P)
    seqLen = as.numeric(seqLenList[[i]])
    if(addTsMut && anyDuplicated(pos)){
      stop("Cannot initialize TS tables from VCF: duplicate retained ",
           "physical positions on chromosome ", names(chrData)[[i]], ".",
           call.=FALSE)
    }
    if(any(pos<0 | pos>seqLen)){
      stop("Cannot initialize TS tables from VCF: retained positions on ",
           "chromosome ", names(chrData)[[i]], " exceed `seqLen`.",
           call.=FALSE)
    }
    haplo = matrix(as.integer(chrData[[i]]$H),
                   nrow=nrow(chrData[[i]]$H))
    out[[i]] = vcfFounderTableCollection(haplo=haplo,
                                         pos=pos,
                                         seqLen=seqLen,
                                         ploidy=as.integer(ploidy),
                                         addTsMut=addTsMut)
  }
  return(out)
}

#' @keywords internal
#' @noRd
.vcfAddPed = function(founderPop, ped, sampleIds, useVCFIds){
  if(is.null(ped)){
    if(!useVCFIds){
      return(founderPop)
    }
    id = as.character(sampleIds)
    mother = father = rep("0", length(id))
  }else if(is.atomic(ped) && is.vector(ped)){
    id = as.character(ped)
    mother = father = rep("0", length(id))
  }else{
    id = as.character(ped[,1])
    mother = as.character(ped[,2])
    father = as.character(ped[,3])
  }
  stopifnot(length(id)==founderPop@nInd,
            !any(duplicated(id)))
  founderPop = new("NamedMapPop",
                   id=id,
                   mother=mother,
                   father=father,
                   founderPop)
  return(founderPop)
}
