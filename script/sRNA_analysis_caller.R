#!/usr/local/R-3.5.1/bin/Rscript --vanilla


#### Loading required packages #####
suppressPackageStartupMessages(library("rmarkdown"))
suppressPackageStartupMessages(library("optparse"))
options(bitmapType='cairo') # TEST IF THIS IS NECESSARY

#### Getting SGE environmental variables #####
#SGE job ID
jobid <- as.character(Sys.getenv("JOB_ID"))
#Node on which job is done
node <- Sys.getenv("HOSTNAME")
#User ID
user <- Sys.getenv("USER")

#### COMMAND LINE ARGUMENTS PARSING ######
option_list <- list(
  make_option(c("-n", "--ncores"),
              type = "integer",
              default = 1,
              help = "Number of core to use on for parallel computing by R."),
  make_option(c("-i", "--sampleInfoBamFiles"),
              type = "character",
              default = NULL,
              help = "Path to the tabulated file describing the samples under analysis."),
  make_option(c("-g", "--genomeGffFile"),
              type = "character",
              default = NULL,
              help = "Path to gff file containing reference genome annotations."),
  make_option(c("-f", "--genomeFastaFile"),
              type = "character",
              default = NULL,
              help = "Path to fasta file containing reference genome sequence."),
  make_option(c("-m", "--miRnaGffFile"),
              type = "character",
              default = NULL,
              help = "Path to gff file containing reference genome MIRNA genes annotations."),  
  make_option(c("-s", "--sstackGffFile"),
              type = "character",
              default = NULL,
              help = "Path to gff file containing inferred sRNA clusters annotations on reference genome."),
  make_option(c("-c", "--comparisonsFile"),
              type = "character",
              default = NULL,
              help = "Path to a ';' separated two-column text file specifying the treatment comparaisons that will be performed for differential accumlation testing."),
  make_option(c("-x", "--xisrnaBai3File"),
              type = "character",
              default = NULL,
              help = "Path to a ';' separated columnar text file specifying information about previously identified xisRNA using BAI3 data."),
  make_option(c("-o", "--outDir"),
              type = "character",
              default = NULL,
              help = "Path to the directory were analysis results will be ultimately written."),
  make_option(c("--scriptsDir"),
              type = "character",
              default = "/data3/projects/pixies/scripts",
              help = "Path to pipeline scripts and config files directory.")
  
)

##### RETRIEVEING PARAMS from optparse #####
myArgs <- parse_args(OptionParser(usage = "%prog", option_list = option_list))


###########################################
##### FOR INTERACTIVE ARGUMENTS PASSING TESTS
# testArgs  <-  c(
#   "--ncores=4",
#   "--sampleInfoBamFiles=/data3/projects/pixies/TMP/mappedBWA/sampleInfo.txt",
#   "--genomeGffFile=/data3/projects/pixies/reference_genomic_data/NB_genome/MSU7.gff3",
#   "--genomeFastaFile=/data3/projects/pixies/reference_genomic_data/NB_genome/NUC_MSU7.fa",
#   "--miRnaGffFile=/data3/projects/pixies/reference_genomic_data/NB_genome/miRBase_MSU7_MIR.gff3",
#   "--sstackGffFile=/data3/projects/pixies/TMP/shortStack/ShortStack_All.gff3",
#   "--comparisonsFile=/data3/projects/pixies/scripts/treatmentsComparisons.csv",
#   "--xisrnaBai3File=/data3/projects/pixies/scripts/xisRNA_BAI3.csv",
#   "--scriptsDir=/data3/projects/pixies/scripts",
#   "--outDir=/data3/projects/pixies/TMP/sRNA_diffAccum"
# )
# myArgs <- parse_args(OptionParser(usage = "%prog", option_list = option_list), args = testArgs)
# jobid <- "1"
###########################################

if(Sys.getenv("NSLOTS") != "") {myArgs$ncores <- as.integer(Sys.getenv("NSLOTS"))} # Assign the number of available cores using the SGE env variable if not empty. Else, it will be the defautl value of optparse..
invisible(dir.exists(myArgs$outDir) || dir.create(myArgs$outDir))

rmdFile <- file.path(myArgs$scriptsDir, "sRNA_analysis.Rmd")


##### Handeling qsub vs local running options #####
if (jobid == "") { # This script in NOT run via qsub (JOB_ID is empty). NOTE however that it is the case if run with qrsh
  ## Copying Rmd file to outDir
  rmdFileInOutDir <- file.path(myArgs$outDir, basename(rmdFile))
  file.copy(from = rmdFile, to = rmdFileInOutDir, overwrite = TRUE, copy.mode = TRUE) 
  
  ## Constructing list of parameters for rmarkdown::render
  renderParams <- list(input = rmdFileInOutDir, # myRmd file to render
                       clean = FALSE,
                       params = list(
                         ncores = myArgs$ncores,
                         sampleInfoBamFiles = myArgs$sampleInfoBamFiles,
                         outDirOnNode = myArgs$outDir,
                         genomeGffFile = myArgs$genomeGffFile,
                         genomeFastaFile = myArgs$genomeFastaFile,
                         miRnaGffFile = myArgs$miRnaGffFile,
                         sstackGffFile = myArgs$sstackGffFile,
                         comparisonsFile = myArgs$comparisonsFile,
                         xisrnaBai3File = myArgs$xisrnaBai3File,
                         scriptsDir = myArgs$scriptsDir
                       )
  )
  do.call(what = rmarkdown::render, args = renderParams)
  
} else { # This script is run via qsub (JOB_ID has a value)
  ## Define paths to dirs on node and create them
  dirOnNode <- file.path("/scratch", paste(user, jobid, sep = "-"))
  invisible(dir.exists(dirOnNode) || dir.create(path = dirOnNode, mode = "770"))
  inDirOnNode <- file.path(dirOnNode, "inDir")
  invisible(dir.exists(inDirOnNode) || dir.create(path = inDirOnNode, mode = "770"))
  outDirOnNode <- file.path(dirOnNode, "outDir")
  invisible(dir.exists(outDirOnNode) || dir.create(path = outDirOnNode, mode = "770"))
  
  ## Move all necessary input files to dirOnNode
  
  # rmd
  rmdFileOnNode <- file.path(outDirOnNode, basename(rmdFile))
  rsyncCmd <- paste0("rsync -avP ", user, "@nas3:", rmdFile, " ", outDirOnNode) 
  !system(rsyncCmd, intern = FALSE) || {cat("rsync unable to transfer file", rmdFile, "to cluster"); quit(status=1, save = "no")}
  
  # NCBI gff
  ####!!! VERY CLUMSY AND CAN BREAK EASILY
  NipNCBIGffFileOnNode <- file.path(dirname(myArgs$genomeGffFile), "Oryza_sativa_genomic_refseq_IRGSP_seqlevels.gff3")
  rsyncCmd <- paste0("rsync -avP ", user, "@nas3:", NipNCBIGffFileOnNode, " ", inDirOnNode) 
  !system(rsyncCmd, intern = FALSE) || {cat("rsync unable to transfer file", NipNCBIGffFileOnNode, "to cluster"); quit(status=1, save = "no")}
  
  # others
  inFiles <- c(sampleInfoBamFiles = myArgs$sampleInfoBamFiles,
               genomeGffFile = myArgs$genomeGffFile,
               genomeFastaFile = myArgs$genomeFastaFile,
               miRnaGffFile = myArgs$miRnaGffFile,
               sstackGffFile = myArgs$sstackGffFile,
               comparisonsFile = myArgs$comparisonsFile,
               xisrnaBai3File = myArgs$xisrnaBai3File
  )
  sampleInfoBamFiles <- read.delim(myArgs$sampleInfoBamFiles, as.is = 1, comment.char = "#")$FileName # Read paths from sampleInfo file
  sampleInfoBamFiles <- c(sampleInfoBamFiles, paste0(sampleInfoBamFiles, ".bai")) # also transferring the corresponding bai files
  inFilesToBeTransfered <- c(sampleInfoBamFiles, inFiles, rmdFile)
  inFilesForParams <- sapply(inFiles, function(x) file.path(inDirOnNode, basename(x)))
  
  for(i in 1:length(inFilesToBeTransfered)) {
    rsyncCmd <- paste0("rsync -avP ", user, "@nas3:", inFilesToBeTransfered[i], " ", inDirOnNode) 
    !system(rsyncCmd, intern = FALSE) || {cat("rsync unable to transfer file", inFilesToBeTransfered[i], "to cluster"); quit(status=1, save = "no")}
  }
  ## Run writeSampleInfo.R to update bam paths on node
  wInfoCmd <- paste0(file.path(myArgs$scriptsDir, "writeSampleInfo.R"),
                     " --inDir=", inDirOnNode,
                     " --infoFile=",  myArgs$sampleInfoBamFiles)
  !system(wInfoCmd, intern = FALSE) || {cat("Unable to create an updated sampleInfo.txt file on node"); quit(status=1, save = "no")}
  
  ## Constructing list of parameters for rmarkdown::render
  renderParams <- list(input = rmdFileOnNode, # myRmd file to render
                       clean = FALSE,
                       params = c(
                         list(ncores = myArgs$ncores, outDirOnNode = outDirOnNode, scriptsDir = myArgs$scriptsDir),
                         as.list(inFilesForParams)
                       )
  )
  cat('-----------------------------------------------------------\n')
  cat(date(), "\n")
  cat("Now rendering rmd file with the following parameters:\n")
  cat(str(renderParams), "\n")
  cat("On node", Sys.getenv("HOSTNAME"), "with job ID", Sys.getenv("JOB_ID"), "\n")
  cat('-----------------------------------------------------------\n')
  
  
  ## Calling rmarkdown::render
  do.call(what = rmarkdown::render, args = renderParams)
  
  
  ## Move all relevant output in outDirOnNode and subfolders to outDir
  rsyncCmd <- paste0("rsync -avP ", outDirOnNode, "/ ", user, "@nas3:", myArgs$outDir) 
  !system(rsyncCmd, intern = FALSE) || {cat("rsync unable to transfer file", toBeTransfered[i], "to cluster"); quit(status=1, save = "no")}
  

  ## Cleanup /scratch folder
  toBeDeleted <- c(dirOnNode, list.files(path = dirOnNode, all.files = TRUE, full.names = TRUE, recursive = TRUE, include.dirs = TRUE))
  (!unlink(toBeDeleted, recursive = TRUE)) || {cat("Unable to delete files on", dirOnNode, ". To be cleaned up manually!")}
  
}
#### Quitting #####
cat('-----------------------------------------------------------\n')
cat(date(), "\n")
cat("Analysis completed!\n")
cat('-----------------------------------------------------------\n')
print(sessionInfo()); quit(status=0, save = "no")
