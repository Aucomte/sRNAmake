#!/usr/local/R-3.5.1/bin/Rscript --vanilla

#### Loading required packages #####
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("GenomicRanges"))

#### COMMAND LINE ARGUMENTS PARSING ######
option_list <- list(
  make_option(c("-p", "--sstackDir"),
              type = "character",
              default = NULL,
              help = "Path to folder containing ShortStack output files."),
  make_option(c("-m", "--miRnaGffFile"),
              type = "character",
              default = NULL,
              help = "Path to mirBase gff file."),
  make_option(c("-o", "--outputGffFile"),
              type = "character",
              default = NULL,
              help = "Path to output gff file.")
)

##### RETRIEVEING PARAMS from optparse #####
myArgs <- parse_args(
  OptionParser(usage = "%prog", option_list = option_list,
               description = "Append information from ShortStack output 'Results.txt' file to the output ShortStack_All.gff3 file and optionally add names of overlapping mirBase pri-mirs.")
  )

###########################################
##### FOR INTERACTIVE ARGUMENTS PASSING TESTS
# testArgs  <-  c(
#   "--sstackDir=/data3/projects/pixies/TMP/shortStack",
#   "--miRnaGffFile=/data3/projects/pixies/reference_genomic_data/NB_genome/miRBase_MSU7_MIR.gff3"
# )
# 
# testArgs  <-  c(
#   "--sstackDir=/data3/projects/pixies/TMP",
#   "--miRnaGffFile=/data3/projects/pixies/reference_genomic_data/NB_genome/miRBase_MSU7_MIR.gff3"
# )
# myArgs <- parse_args(OptionParser(usage = "%prog", option_list = option_list), args = testArgs)
###########################################

# Check parameter values
invisible(dir.exists(myArgs$sstackDir) || stop(paste0("Unable to find ", myArgs$sstackDir)))

resultsFile <- file.path(myArgs$sstackDir, "Results.txt")
gffFile <- file.path(myArgs$sstackDir, "ShortStack_All.gff3")
outputGffFile <- myArgs$outputGffFile

# Load content of files
clustersGR <- rtracklayer::import.gff3(gffFile)

clustersInfo <- read.table(file = resultsFile,
                           sep = "\t",
                           header = TRUE,
                           comment.char = "",
                           row.names = NULL,
                           stringsAsFactors = FALSE)

# Verify that cluster IDs corresponds between the two files
if (!all(clustersGR$ID ==clustersInfo$Name)) stop("Cluster lists of cluster IDs in gff and Results.txt file do not match!! ")


# Populate mcols of GR with info from data.frame
mcols(clustersGR) <- mcols(clustersGR)[, c("source", "type")] # keep only source and type columns
mcols(clustersGR) <- cbind(mcols(clustersGR), clustersInfo[, -1])
mcols(clustersGR)$Strand <- NULL # Necessary otherwise gff3 import is unable to determine strand.

if (!is.null(myArgs$miRnaGffFile) && file.exists(myArgs$miRnaGffFile)) {
  # Load mirBase gff file content
  osMirAnnot <- rtracklayer::import.gff3(myArgs$miRnaGffFile)
  osMirAnnot <- osMirAnnot[values(osMirAnnot)[,"type"] == "miRNA_primary_transcript"] # Keep only MIRNA primary transcript GRanges (discard miRNA mapping info)
  
  # Find and Append information on mirBase MIRNA genes overlapping sRNA loci
  olaps <-  findOverlaps(query = clustersGR, subject = osMirAnnot,
                         type= "any", ignore.strand = TRUE, select= "first")
  mcols(clustersGR) <- cbind(mcols(clustersGR),
                             data.frame(pri_miRNA = osMirAnnot$Name[olaps], stringsAsFactors = FALSE))
} else {cat("No miRNA annotation gff file was specified or could not be found at the specified location!\n")}

# Overwrite gffFile with the version containning full info

rtracklayer::export.gff3(
  object = clustersGR,
  con = outputGffFile
)

#### Quitting #####
cat('-----------------------------------------------------------\n')
cat(date(), "\n")
cat("Analysis completed!\n")
cat('-----------------------------------------------------------\n')
print(sessionInfo()); quit(status=0, save = "no")
