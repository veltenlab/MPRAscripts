#!/usr/bin/env Rscript

# Process MPRA screen data
# Set the library folder
.libPaths("/path/to/snakemake/environment/MPRA_R_processing/lib/R/library")

# parse command line arguments ---------------------------------------------------------------------

library("optparse")

# create arguments list
option_list = list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Folder with the input files", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output .rds object", metavar = "character"), 
  make_option(c("-c", "--cores"), type = "integer", default = 1,
              help = "Number of cores to use in parallel", metavar = "integer"), 
  make_option(c("-b", "--bccrsass"), type = "character", default = NULL,
              help = "CRS-BC association file", metavar = "character"),
  make_option(c("-n", "--normalizetoinput"), type = "logical", default = FALSE,
              help = "Normalize to input", metavar = "logical")
)

# parse arguments
opt_parser = OptionParser(option_list = option_list, add_help_option = FALSE)

opt = parse_args(opt_parser)

# function to check for required arguments
check_required_args <- function(arg, opt, opt_parser) {
  if (is.null(opt[[arg]])) {
    print_help(opt_parser)
    stop(arg, " argument is required!", call. = FALSE) 
  } 
}
# check that all required parameters are provided
required_args <- c("input", "output", "bccrsass", "normalizetoinput")
for (i in required_args) {
  check_required_args(i, opt = opt, opt_parser = opt_parser) 
}

# Running the script
source("/path/to/project/folder/Pipelines/MPRA/3_Report_screen/scripts/get_Screen.R")
library(parallel)

outdir <- gsub("/[^/]+$", "",opt$output)
if (!dir.exists(file.path(outdir))) {dir.create(file.path(outdir), showWarnings = T, recursive = T )}

#Select folder
folders <- list.files(opt$input,full.names = T,pattern = "^State")
names(folders) <- list.files(opt$input,full.names = F,pattern = "^State")

experiments  <- mcmapply(get_Screen, folders,
                      MoreArgs = list(normalize.to.input = opt$normalizetoinput,
                                      bc_crs_ass= opt$bccrsass), 
                      SIMPLIFY = F , mc.cores = opt$cores)

saveRDS(experiments, file = opt$output )
