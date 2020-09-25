


myLoadInstall <- function(packageNames) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    cat("### The 'BiocManager' package is not available on the system and will be installed!!")
    install.packages("BiocManager", repos ="https://cloud.r-project.org")
  }
  # List packages that need to be installed
  isInstalled <- packageNames %in% installed.packages()[, "Package"]
  toBeInstalled <- packageNames[!isInstalled]
  # Installing if necessary
  if(length(toBeInstalled) > 0) {
    cat("\n### The following packages are not available on", .libPaths(),
        "and will be installed using BiocManager::install() prior to loading:\n -", paste0(toBeInstalled, collapse = "\n- "),
        "\n")
    BiocManager::install(toBeInstalled, update = FALSE, ask = FALSE)
  } else {
    cat("\n### All required packages are available on your system and will be loaded\n")
  }
  # Loading the list of packages
  cat("\n### Now loading the required packages:", paste0(packageNames, collapse = " "), "\n")
  invisible(
    sapply(packageNames, FUN = function(packageName) suppressPackageStartupMessages(library(package = packageName, character.only = TRUE)))
  )
}





#### O. sativa chromosome (IRGSP1.0) Seqinfo object ####
os.chr.sizes <- c(
  "Chr1" = 43270923,
  "Chr2" = 35937250,
  "Chr3" = 36413819,
  "Chr4" = 35502694,
  "Chr5" = 29958434,
  "Chr6" = 31248787,
  "Chr7" = 29697621,
  "Chr8" = 28443022,
  "Chr9" = 23012720,
  "Chr10" = 23207287,
  "Chr11" = 29021106,
  "Chr12" = 27531856,
  "chrC" = 134525,
  "chrM" = 490520,
  "ChrUn" = 633585,
  "ChrSy" = 592136)

library(GenomicRanges)
chrominfo <- Seqinfo(seqnames = names(os.chr.sizes),
                     seqlengths=os.chr.sizes ,
                     isCircular=rep(FALSE, length(os.chr.sizes)),
                     genome="IRGSP-1.0")
isCircular(chrominfo)[c("chrC", "chrM")] <- TRUE




#####  Misc. NGS and DE analysis function definitions  #####
## Including code borrowed from package vignettes


# Constructing a custom srFilter function to keep reads with mapping qualities
# satisfying <= crit1 and >= crit2 criteria
makeMyAlQualFilter <- function(x, crit1, crit2) {
  # This function returns a srFilter function using
  # crit1 and crit2 criteria as treshold for selecting reads in object x that
  # have alignements qualities satisfying these criteria
  # Note that in principle, this should also filter out
  # reads that have no alignment quality i.e. that are unmapped.
  if (crit1 > crit2) stop("Check criteria values: crit1 > crit2")
  if (missing(crit1) | missing(crit2)) stop("Missing at least one criteria")
  myFilter <- function (x) {
    isNotNA <- !is.na(quality(alignQuality(x)))
    desiredQual <- quality(alignQuality(x)) <= crit1 |  quality(alignQuality(x)) >= crit2
    return(isNotNA & desiredQual)
  }
  return(srFilter(myFilter, name = "myWeirdAlQualFilter"))
}


makeMyReadWidthFilter <- function(x, crit1, crit2) {
  # This function returns a srFilter function using
  # crit1 and crit2 criteria as treshold for selecting reads in object x that
  # have a width included within the interval defined by crit1 and crit2
  # Note that in principle, this should also filter out
  # reads that have no width if this ever exist.
  if (crit1 > crit2) stop("Check criteria values: crit1 > crit2")
  if (crit1 < 0 | crit2 < 0) stop("At least a criterias is negative!!")
  if (missing(crit1)) {crit1 <- crit2}
  if (missing(crit2)) {crit2 <- crit1}
  myFilter <- function (x) {
    isNotNA <- !is.na(width(x))
    desiredWidht <- ((width(x) >= crit1) & (width(x) <= crit2))
    return(isNotNA & desiredWidht)
  }
  return(srFilter(myFilter, name = "MyReadWidthFilter"))
}




# Compute short reads counts for each feature and return a RangedSummarizedExperiment object
#!!! BEWARE of the which parameter !!!!
# NB: cannot use SR filters to load gappedReads fed to summarizeOverlaps. Would need to get the counts using another
# function that takes AlignedReads...
computeCountsInFeatures <- function(features = features, reads = bamList, bamParam) {
# Define fonction that perform the workflow and return a countDataSet for processing by DESeq
  library("GenomicRanges") ; library("Rsamtools")
# What are the filters we use for loading the reads
  flag <- scanBamFlag(isUnmappedQuery = FALSE, isDuplicate=NA)
  if (missing(bamParam)) {
    bamParam <- ScanBamParam(flag = flag, simpleCigar = FALSE, reverseComplement = TRUE)
  }

# Create the table of counts for genomic features of interest
# For now use summarizeOverlaps because it is convenient but beware
# of the mode. At some point countOverlaps may be less restrictive.
  message("Starting to compute read counts on query features...")
  olapCounts <- summarizeOverlaps(features = features, reads = reads,
      mode = "Union", ignore.strand = TRUE, param = bamParam, singleEnd = TRUE)
  message("Done!")

# Add info on experimental treatments. A phenoData-like DataFrame
colData(olapCounts) <- cbind(elementMetadata(reads), colData(olapCounts))
rownames(colData(olapCounts)) <- rownames(elementMetadata(reads))

return(olapCounts)
}



# DESeq2-associated functions
runDESeq2 <- function(se, varInt, batch=NULL, comps=NULL,
                      locfunc="median", fitType="parametric", sizeFactEstimMethod = "ratio",
                      pAdjustMethod="BH",
                      cooksCutoff=TRUE, independentFiltering=TRUE, alpha=0.05,
                      lfcThreshold = 0, altHypothesis = "greaterAbs", ...){
  # Wrapper to run DESeq2
  #
  # Wrapper to run DESeq2: create the \code{DESeqDataSet}, normalize data, estimate dispersions, statistical testing...
  #
  #  se \code{RangedSummarizedExperiment} object with a properle set colData slot
  #  varInt name of the factor of interest (biological condition)
  #  batch batch effect to take into account (\code{NULL} by default)
  #  comps data frame listing all the pairs of treatment levels comparisons to test (\code{NULL} by default)
  #  locfunc \code{"median"} (default) or \code{"shorth"} to estimate the size factors
  #  fitType mean-variance relationship: "parametric" (default) or "local"
  #  pAdjustMethod p-value adjustment method: \code{"BH"} (default) or \code{"BY"} for instance
  #  cooksCutoff outliers detection threshold (TRUE to let DESeq2 choosing it or FALSE to disable the outliers detection)
  #  independentFiltering \code{TRUE} or \code{FALSE} to perform the independent filtering or not
  #  alpha significance threshold to apply to the adjusted p-values
  #  ... optional arguments to be passed to \code{nbinomWaldTest()}
  # return A list containing the \code{dds} object (\code{DESeqDataSet} class), the \code{results} objects (\code{DESeqResults} class) and the vector of size factors
  # author Hugo Varet + modifications Seb Cunnac

  library("DESeq2")

  # building dds object
  dds <- DESeq2::DESeqDataSet(se=se,
                              design=formula(paste("~", ifelse(!is.null(batch), paste(batch,"+"), ""), varInt)))
  cat("Design of the statistical model:\n")
  cat(paste(as.character(design(dds)),collapse=" "),"\n")

  # normalization
  dds <- DESeq2::estimateSizeFactors(dds, type = sizeFactEstimMethod, locfunc=eval(as.name(locfunc)))
  cat("\nNormalization factors:\n")
  print(sizeFactors(dds))

  # estimating dispersions
  dds <- DESeq2::estimateDispersions(dds, fitType=fitType)

  # tests for significance of coefficients in a Negative Binomial GLM
  dds <- DESeq2::nbinomWaldTest(dds, ...)

  # statistical testing: perform all the comparisons between the levels of varInt

  if (is.null(comps)) {
    results <- list()
    for (comp in combn(nlevels(colData(dds)[,varInt]), 2, simplify=FALSE)){
      levelRef <- levels(colData(dds)[,varInt])[comp[1]]
      levelTest <- levels(colData(dds)[,varInt])[comp[2]]
      results[[paste0(levelTest,"_vs_",levelRef)]] <- results(dds, contrast=c(varInt, levelTest, levelRef),
                                                              pAdjustMethod=pAdjustMethod,
                                                              cooksCutoff=cooksCutoff,
                                                              independentFiltering=independentFiltering,
                                                              alpha=alpha,
                                                              lfcThreshold = lfcThreshold,
                                                              altHypothesis = altHypothesis,
                                                              parallel = TRUE)
      cat(paste("Comparison", levelTest, "vs", levelRef, "done\n"))
    }
  } else {
    results <- bplapply(1:nrow(comps), function(i, dds, comps, pAdjustMethod, cooksCutoff, independentFiltering, alpha){
      comp <- unlist(comps[i, ])
      comp <- c(comp, paste0(comp[1], "vs", comp[2]))
      res <- DESeq2::results(dds, contrast = c(varInt, comp[1], comp[2]),
                             pAdjustMethod=pAdjustMethod, cooksCutoff=cooksCutoff,
                             independentFiltering=independentFiltering, alpha=alpha,
                             lfcThreshold = lfcThreshold, altHypothesis = altHypothesis,
                             parallel = FALSE)
      cat(paste("Comparison", comp[1], "vs", comp[2], "done\n"))
      return(res)
    }, dds=dds, comps = comps, pAdjustMethod=pAdjustMethod, cooksCutoff=cooksCutoff, independentFiltering=independentFiltering, alpha=alpha)
    names(results) <- paste0(comparisons[[1]], "vs", comparisons[[2]])
  }
  return(list(dds=dds,results=results,sf=sizeFactors(dds)))
}


exportDESeq2Results <- function(out.DESeq2, group, alpha=0.05, export=TRUE){
  # Original function from:
  # https://github.com/PF2-pasteur-fr/SARTools/blob/master/R/exportResults.DESeq2.R
  # Modified to fit my needs
  dds <- out.DESeq2$dds
  fpmData <- DESeq2::fpm(dds, robust = TRUE)
  results <- out.DESeq2$results

  # comptages bruts et normalis?s
  counts <- data.frame(Id=rownames(counts(dds)), counts(dds), round(fpmData))
  colnames(counts) <- c("Id", colnames(counts(dds)), paste0("fpm.", colnames(counts(dds))))
  # baseMean avec identifiant
  bm <- data.frame(Id=rownames(results[[1]]),baseMean=round(results[[1]][,"baseMean"],2))
  # merge des info, comptages et baseMean selon l'Id
  base <- merge(counts, bm, by="Id", all=TRUE)
  tmp <- base[,paste("fpm", colnames(counts(dds)), sep=".")]
  for (cond in levels(group)){
    base[,cond] <- round(apply(as.data.frame(tmp[,group==cond]),1,mean),0)
  }

  complete <- list()
  for (name in names(results)){
    complete.name <- base

    # ajout d'elements depuis results
    res.name <- data.frame(Id=rownames(results[[name]]),
                           FoldChange=round(2^(results[[name]][,"log2FoldChange"]), 3),
                           log2FoldChange=round(results[[name]][,"log2FoldChange"], 3),
                           stat=round(results[[name]][,"stat"], 3),
                           pvalue=results[[name]][,"pvalue"],
                           padj=results[[name]][,"padj"])
    complete.name <- merge(complete.name, res.name, by="Id", all=TRUE)
    # ajout d'elements depuis mcols(dds)
    mcols.add <- data.frame(Id=rownames(counts(dds)),dispGeneEst=round(mcols(dds)$dispGeneEst,4),
                            dispFit=round(mcols(dds)$dispFit,4),dispMAP=round(mcols(dds)$dispMAP,4),
                            dispersion=round(mcols(dds)$dispersion,4),betaConv=mcols(dds)$betaConv,
                            maxCooks=round(mcols(dds)$maxCooks,4))
    complete.name <- merge(complete.name, mcols.add, by="Id", all=TRUE)
    complete[[name]] <- complete.name

    if (export){
      # s?lection des up et down
      up.name <- complete.name[which(complete.name$padj <= alpha & complete.name$betaConv & complete.name$log2FoldChange>=0),]
      up.name <- up.name[order(up.name$padj),]
      down.name <- complete.name[which(complete.name$padj <= alpha & complete.name$betaConv & complete.name$log2FoldChange<=0),]
      down.name <- down.name[order(down.name$padj),]

      # exports
      name <- gsub("_","",name)
      write.table(complete.name, file=paste0("tables/",name,".complete.txt"), sep="\t", row.names=FALSE, dec=".", quote=FALSE)
      write.table(up.name, file=paste0("tables/", name,".up.txt"), row.names=FALSE, sep="\t", dec=".", quote=FALSE)
      write.table(down.name, file=paste0("tables/", name,".down.txt"), row.names=FALSE, sep="\t", dec=".", quote=FALSE)
    }
  }

  return(complete)
}







# Coerce a data container from segmentSeq to a DESeq one
lociData2CountDataSet <- function(lociDataObject) {
  # Coerces a lociData object from package segmentSeq
  #  to a countDataSet object from package DESeq
  library("DESeq") ; library("segmentSeq")
  if (!class(lociDataObject) == "lociData") {stop("The provided argument is not of class 'lociData'.")}
  # Get individual elements from the lociData object
  names <-  GRanges2GBrowser(lociDataObject@coordinates)
  myCountData <- lociDataObject@data
  rownames(myCountData) <- names
  myFeatureData <- as.data.frame(lociDataObject@coordinates)
  rownames(myFeatureData) <- names
  myFeatureData <- AnnotatedDataFrame(myFeatureData)
  # Put them together in a countDataSet object
  cdset <- newCountDataSet(countData= myCountData,
      featureData = myFeatureData,
      conditions= lociDataObject@replicates)
  return(cdset)
}



# Plotting function definition from DESeq manual
# and take a look at the disp estimates...
plotDispEsts <- function( cds ) {
  # Plots the per-gene estimates against the normalized mean expressions
  # per gene, and then overlay the fitted curve in red.
  # Copied from DESq manual.
  plot(rowMeans( counts( cds, normalized=TRUE ) ),
      fitInfo(cds)$perGeneDispEsts,
      pch = '.', log="xy" )
  xg <- 10^seq( -.5, 5, length.out=300 )
  lines( xg, fitInfo(cds)$dispFun( xg ), col="red" )
}
# A MA plotting function from DESeq
plotDE <- function( res ) {
  # Plots normalised mean versus log2 fold change.
  # This plot is sometimes also called the MA-plot for the contrast untreated versus treated.
  # Copied from DESq manual.
  plot(res$baseMean,
      res$log2FoldChange,
      log="x", pch=20, cex=.3,
      col = ifelse( res$padj < .1, "red", "black" ) )
}


# Filter output for diff accum elements in DESeq
nbinomTestFilt <- function(cds, condA, condB, min.lfc = log2(2), maxAdjPval = 0.05) {
  library("DESeq")
  res <- nbinomTest(cds=cds, condA=condA, condB=condB)
  res[!is.na(res$log2FoldChange) & res$padj <= maxAdjPval & abs(res$log2FoldChange) >= min.lfc , ]
}


# Perform treatments comparisons and merge results in a data frame
mergeNBinomTestResults <- function (cds, comps, min.lfc = 2, maxAdjPval = 0.05) {
  df <- data.frame()
  for (i in 1:nrow(comps)) {
    # perform the tests for each comparison of treatment
    comp <- unlist(comps[i, ])
    comp <- c(comp, paste0(comp[1], "vs", comp[2]))
    diffExpr <- do.call(
      nbinomTestFilt,
      list(condA = comp[1], condB = comp[2], cds = cds, min.lfc = min.lfc, maxAdjPval = maxAdjPval)
    )

    message(paste("Analysis for", comp[3], "is done.",sep = " "))
    if(nrow(diffExpr) != 0) {
      diffExpr <- cbind(diffExpr, comparison = comp[3]) # appen a column with description of the comparaison to the resulting data frame
      df <- rbind(df, diffExpr) # merge this data frame to the global table of results
    } else {next}
  }
  return(df)
}



performDESeqBinCompPipeline <- function(cdsetOlap, min.lfc = 0, maxAdjPval = 1,  comparisons,
                                        method = "pooled", fitType = "local") {
  library("DESeq")
# Estimate the effective library size, i.e., correct for seq depth discrepencies
  cdsetOlap <- DESeq::estimateSizeFactors(cdsetOlap)
  print(sizeFactors(cdsetOlap))
# Estimate dispersion, i.e., biological variation
  cdsetOlap <- DESeq::estimateDispersions(cdsetOlap, method = method, sharingMode="fit-only", fitType = fitType)
  #cdsetOlap <- DESeq::estimateDispersions(cdsetOlap, method="per-condition", sharingMode="maximum")

  # save(cdsetOlap, file = "cdsetOlap.Rdata") # save an image of the CountDataSet object on disk

   print(plotDispEsts(cdsetOlap))
  ##   Calling differential expression
# Actual tests for differential expression in the case of RNA-Seq
# perform comparisons for conditions of interest.
  ## Merge together the df appending a column indicating the appropriate contraste of origin and save in DB or in a tab file

# Perform comparison and return the resulting data frame
  return(mergeNBinomTestResults(cds = cdsetOlap, comps = comparisons, min.lfc = min.lfc, maxAdjPval = maxAdjPval))
}




# Coercion utilities to switch between GRanges and genome browser location formats
# THESE are currently obsolete because the GenomicRanges package now provides
# equivalent functions.

GRanges2GBrowser <- function(x, margin = 0) {
  # construct a vector of genome positions following the synthax used in genome browsers
  # using a GRanges object (strand blind)
  paste(seqnames(x), ":", start(x)-margin, "-", end(x)+margin, sep = "")}

GBrowser2GRanges <- function (x, margin = 0) {
  # Return a GRanges when provided a vector of
  # genome positions that follow the synthax used in genome browsers: "chr:start-end"
  spl <- strsplit(x, "\\s*:\\s*", perl = TRUE) # split on ":" removing the spaces before and after
  chr <- sapply(spl, function(x) x[1])
  coord <- t(sapply(spl, function(x) as.integer(unlist(strsplit(x[2], "\\s*-\\s*", perl = TRUE)))))
  GRanges(seqnames = Rle(chr),
      ranges = IRanges(start = (coord[,1]-margin), end = (coord[,2]+margin)))
}


# TO SOME EXTENT THIS is also obsolete because several packages such as plyranges
# provide similar utilities
addRefAnnot <- function(querGR, refGR, mcolsToKeep = NULL, prefix = "OlapinRef", ignore.strand = TRUE) {
  ## Takes a query and a reference GRange and output the query GRanges with additional info from the overlapping
  ## reference GRange mcols columns as specified in the vector of column names 'mcolsToKeep'.
  if (is.null(mcolsToKeep)) {
    warning("The 'mcolsToKeep' parameter has not been specified, will include all columns in the 'refGR' object!")
  }

  requestedColsAvailable <- mcolsToKeep %in% colnames(mcols(refGR))
  if (!any(requestedColsAvailable)) {
    warning("None of the columns listed in 'mcolsToKeep' parameter are available in the 'refGR' columns! Adding only 'refGR' names info if available!")
    mcolsToKeep <- NULL
  } else if (!all(requestedColsAvailable)) {
    warning(paste("Some columns listed in 'mcolsToKeep' parameter are not available in the 'refGR' columns! Adding info only from the available ones:\n",
                  paste(mcolsToKeep, collapse = ", "))
    )
    mcolsToKeep <- mcolsToKeep[requestedColsAvailable]
  }
  # Look for overlapping GRanges
  hitsList <- findOverlaps(querGR, refGR, minoverlap = 1L,
                           ignore.strand = ignore.strand,
                           type = "any",
                           select = "first")
  refGRHitIndexes <- hitsList[!is.na(hitsList)]

  # Get mcols info of overlapping subject GRanges
  hitsInfo <- mcols(refGR)[refGRHitIndexes, mcolsToKeep, drop = FALSE]

  if (!is.null(names(refGR))) hitsInfo <- cbind(DataFrame(Names = names(refGR)[refGRHitIndexes]),
                                                hitsInfo)
  colnames(hitsInfo) <- paste0(prefix, colnames(hitsInfo))

  hitsInfo <- DataFrame(lapply(hitsInfo,  as.character)) # This is the only brutal way I found to deal with '<CharacterList>'-like objects that were messing up the merge
  hitsInfo <- as.data.frame(hitsInfo, optional = TRUE)
  # Append info on xisRNAs in the sRNALociOverlapGenes table
  info <- as.data.frame(mcols(querGR), optional = TRUE)
  cols <- colnames(hitsInfo)
  for (i in 1:length(cols)) {
    info[!is.na(hitsList), cols[i]] <- hitsInfo[,i, drop = FALSE]
  }

  if (length(querGR) != nrow(info)) warning("The length of the returned GRanges object is different from the length of querGR!")

  mcols(querGR) <- info
  return(querGR)
}


# Utility to extract DNA sequences from an had oc OsMSU6 BSgenome package
getSeqOsMSU7 <- function(gr) {
  # Returns the DNA sequence corresponding the OsMSU6 genome location specified in the GenomicRanges "gr"
  library("BSgenome.Osativa.MSU.7")
  seqlevels(gr) <- chartr("|", "_", seqlevels(gr))
  getSeq(OsMSU7, gr, as.character=TRUE)
}



# Look at the accumulation of unique read sequences in sRNA loci of interest
readSeqDistributionInGRanges <- function(gr, strd = NULL, bams) {
# Return a data frame with unique read sequences mapped to a genomic loci specified as a GRange
# with counts of occurences of each sequence broken down by libraies in bamList

# Check the value of the strand argument and assign both strands if NULL or "*"
  if (!(all(strd %in% c("+", "-", "*") | is.null(strd)))) {stop("Provided an unrecognized designation for the strand!\n")}
  if (is.null(strd) | strd %in% "*") {strd <- c("+", "-")}
# Fetch the reads corresponding to a selected location (GRange)
  readsInMyGRange <- inspectReadsInGRanges(bamfile = bams, gr= gr,
      returndf = FALSE, printPlot = FALSE, returnAlns = TRUE)
# Compute the counts of occurence for each unique sequences in the set of reads and returns a list for each library
# and assemble a data frame with a column specifying the name of libraries.
  readDistries <- data.frame()
  for (i in seq_along(readsInMyGRange)) {
    GapReads <- readsInMyGRange[[i]]
    readsDistri <- tables(qseq(GapReads[strand(GapReads) %in% strd]))[["top"]]

    if (length(readsDistri) == 0) {next} else {
      df <- cbind(
          data.frame(sequence = names(readsDistri)),
          length = nchar(names(readsDistri)),
          counts = as.numeric(readsDistri),
          library = names(readsInMyGRange[i]))
      readDistries <- rbind(readDistries, df)
    }
  }
  return(readDistries)
}



# VERY INTERESTING PLOTTING FUNCTIONS

# chartr("/", "\\", file.path(getwd(), paste(title[i], ".jpeg", sep = "")))



saveIGBview <- function(gr, width = NULL) {
  # Generate IGB genome browser image corresponding to the genome location
  # specified in the GenomicRanges gr.
  loc <- unique(gr)
  title <- chartr(":", "_", names(loc))
  if (is.null(width)) {loc <- GRanges2GBrowser(loc)} else {
    loc <- GRanges2GBrowser(resize(loc, width = (width(loc) + width), fix = "center", use.names = TRUE))}
  if (any(is.null(title))) {title <- loc}
  # Generate the script controlling display in IGB
  script <- vector()
  for (i in seq_along(loc)) {
    script <- c(script,
        paste("goto ", loc[i], sep = ""),
        "refresh",
        "sleep 200",
        paste("snapshotmainView ", file.path(getwd(), paste(title[i], ".jpeg", sep = "")), sep = "")
    )
  }
  writeLines(script, con = "IGBscript.igb")
  # Send the html request to the IGB bookmark server
  library("RCurl")
  print(try(getForm(uri = "http://localhost:7085/UnibrowControl",
              scriptfile = paste(file.path(getwd(), "IGBscript.igb"), sep = ""),
              .opts = list(connecttimeout = 1, timeout.ms = 1000))))
  #file.remove("IGBscript.igb")
}



saveIGVbatchFile <- function(gr, width = NULL, file = "IGVscript.igb") {
  # Generate IGV batch command text file for saving browser image corresponding to the genome location
  # specified in the GenomicRanges gr
  loc <- unique(gr)
  title <- chartr(":", "_", names(loc))
  if (is.null(width)) {loc <- GRanges2GBrowser(loc)} else {
    loc <- GRanges2GBrowser(resize(loc, width = (width(loc) + width), fix = "center", use.names = TRUE))}
  if (any(is.null(title))) {title <- loc}
  # Generate the script controlling display in IGV
  script <- c("setSleepInterval 500",
      paste("snapshotDirectory", getwd(), sep = " ")
  )
  for (i in seq_along(loc)) {
    script <- c(script,
        paste("goto", loc[i], sep = " "),
        paste("snapshot", paste(title[i], ".jpeg", sep = ""), sep = " ")
    )
  }
  writeLines(script, con = file)
}




myPlotGenome <- function (aD, locData, GR, margin = 0, ...) {
# A slightly modified version of plotGenome from the segmentSeq package.
# It can take a GenomicRanges object to define the limits of the region to plot.
  if (!length(GR) == 1) stop("myPlotGenome takes only GRanges objects with a single location")
  plotGenome(aD, locData,
      chr = as.character(seqnames(GR)),
      limits = c(start(GR)-margin, end(GR)+margin),
      cex.axis = 0.5, ...)
}



inspectReadsInGRanges <- function(bamfile, gr, bamFilter, plotWhat, returndf = TRUE, printPlot = TRUE, returnAlns = FALSE, ...) {
  # If used in large regions or in read dense ones, expect very poor performances because it is not optimized
  # load reads from bamfile in the regions specified by "gr" using readBamGappedReads
  # make a data frame with the information contained in the BAM files regarding the alignments
  # it is optionally returned (returndf)
  # Optionally (printPlot) return an histogram of the distribution of the counts
  # for the variable(s) in the "plotWhat" vector. Currently plotWhat handles up to two variables by plotting,
  # counts data for plotWhat[1] levels by plotWhat[2] levels, one plot per BAM file in bamfile>
  # Currently available variables are:
  # NM  Edit distance
  # X0  Number of best hits
  # XT  Type: Unique/Repeat/N/Mate-sw
  # MD  Mismatching positions/bases
  # XA  Alternative hits; format: (chr,pos,CIGAR,NM;)
  # width  read length
  # strand  genomic strand to which the read aligns
  # firstBase nucleotide at the first (5')  position of the read

  require(GenomicAlignments) ; require(BiocParallel)

  if (length(gr) > 1) message("The returned analysis will include several genomic locations!")
  if (missing(bamFilter)) {
    warning("Using built-in default ScanBamParameter object to filter reads\non provided GenomicRanges.\n")
    flag <- scanBamFlag(isPaired = NA, isProperPair = NA, isUnmappedQuery = FALSE,
        hasUnmappedMate = NA, isMinusStrand = NA, isMateMinusStrand = NA,
        isFirstMateRead = NA, isSecondMateRead = NA, isSecondaryAlignment = NA,
        isNotPassingQualityControls = NA, isDuplicate = NA)
    tag <- c("NM", "X0", "XT","MD")
    what <- character(0) # c("cigar", "rname", "strand", "pos", "qwidth")
    which <- gr
    bamFilter <- ScanBamParam(flag = flag, simpleCigar = FALSE,
        reverseComplement = FALSE, tag = tag,
        what = what, which = which)
  }

  yieldSize(bamfile) <- 100 # read everything for the time being
  alnGapList <- bplapply(bamfile, readGappedReads, param=bamFilter, BPPARAM = MulticoreParam(1))

  
  if(printPlot | returndf){
    df <- data.frame()
    for (i in 1:length(alnGapList)) {
      if (!length(alnGapList[[i]]) == 0) {
        table <- cbind(as.data.frame(mcols(alnGapList[[i]])[,tag]),
                       width = qwidth(alnGapList[[i]]),
                       #firstBase = as.character(subseq(qseq(alnGapList[[i]]), start=1, end = 1)),
                       readSeq = as.character(qseq(alnGapList[[i]])),
                       strand = as.character(strand(alnGapList[[i]])),
                       sample = names(alnGapList[i]))
        df <- rbind(df, table)
      }
    }
  }


  if (printPlot) {
    x <- as.data.frame(xtabs(
            formula = formula(paste("~", paste(plotWhat, collapse = "+"), " + sample", sep = "")),
            data = df),
        responseName = "count")

    print(barchart(formula(paste("count~", plotWhat[1], " | sample", sep = "")),
            groups = if (length(plotWhat) == 2) {eval(as.symbol(plotWhat[2]))} else {NULL},
            data = x,
            xlab = plotWhat[1],
            auto.key = list(title = plotWhat[2], cex.title = 0.8),
            as.table = TRUE,
            scales = list(y = list(relation = "free"), x = list(relation = "same")),
            main  = if (!is.null(names(gr))) {names(gr)},
            ...
        ))
  }
  if (returndf) {return(df)}
  if (returnAlns) {alnGapList}
}



plotsRNAbyFile <- function(grs, plotWhatLst = c("width", "strand", "NM", "X0", "XT", "snapshot"), ...) {
  for (i in seq_along(grs)) {
    if ("snapshot" %in% plotWhatLst) {
      ## Does not worextentedGR <- resize(grs[i], width = 200, fix = "center", use.names = TRUE)k well because Snapshot has a weird behavior, I suppose
      if (width(grs[i]) < 50) {
        extentedGR <- resize(grs[i], width = 200, fix = "center", use.names = TRUE)} else {
        extentedGR <- grs[i]}
      print(names(grs[i]))
      print(extentedGR)
      jpeg(file = paste(names(grs[i]), "snapshot.jpg", sep = "_"))
      try(
          print(Snapshot(files = bamList, range = extentedGR, ignore.strand=FALSE, currentFunction="multifine_coverage")
          )
      )
      dev.off()
      plotWhatLst <- plotWhatLst[!plotWhatLst %in% "snapshot"] # remove "snapshot" from the vector
      if (length(plotWhatLst) == 0) {break}
    }

    for (plotWhat in plotWhatLst) {
      Cairo(type="svg", file = paste(chartr(":", "_", names(grs[i])), "_", paste(plotWhat, collapse ="by"), ".svg", sep = ""),
          width = 1440, height = 760)
      print(inspectReadsInGRanges(gr = grs[i],
              plotWhat = plotWhat, returndf = F, printPlot = T, returnAlns = F,
              strip = strip.custom(par.strip.text = list(cex = 0.8)), ...))
      dev.off()
    }

  }
}



plotsRNAinTable <- function(data, plotWhat = c("width", "strand", "NM", "X0", "XT"), by = "sample", ...) {
  x <- as.data.frame(xtabs(
    formula = formula(paste("~", paste(plotWhat, collapse = "+"), " + ", by, sep = "")),
    data = data),
    responseName = "count")
  # ncol <- ifelse(nlevels(x$by) <4, ceiling(sqrt(nlevels(x$by))), 3)
  # nrow <- ceiling(nlevels(x$by)/ncol)
  # to be added as a parameter in the barchart function: layout = c(ncol,nrow)

  print(lattice::barchart(formula(paste("count~", plotWhat[1], " | ", by, sep = "")),
                 groups = if (length(plotWhat) == 2) {eval(as.symbol(plotWhat[2]))} else {NULL},
                 data = x,
                 xlab = plotWhat[1],
                 auto.key = list(title = plotWhat[2], cex.title = 0.8),
                 as.table = TRUE,
                 scales = list(y = list(relation = "free"), x = list(relation = "same")),
                 ...
        )
  )
}


colorsForFactor <-  function(var, eerFun = "paletteer_d", palFullName = "cartography::sand.pal", direct = 1) {
  # Defining a function to simplify assigning colors to annotation factor levels
  if(is.character(eerFun)){
    fn <- strsplit(eerFun, "::")[[1]]
    eerFun <- if(length(fn)==1) {
      get(fn[[1]], envir=parent.frame(), mode="function")
    } else {
      get(fn[[2]], envir=asNamespace(fn[[1]]), mode="function")
    }
  }
  pal <- as.character(do.call(eerFun, list(palette = palFullName, n = nlevels(var), direction = direct)))
  names(pal) <- levels(var)
  return(pal)
}


