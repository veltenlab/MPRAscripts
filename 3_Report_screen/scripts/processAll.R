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
  make_option(c("--infile"), type = "character", default = NULL,
              help = "Input .rds object", metavar = "character"),   
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output Folder", metavar = "character"), 
  make_option(c("-c", "--cores"), type = "integer", default = 1,
              help = "Number of cores to use in parallel", metavar = "integer"),     
  make_option(c("-e", "--returnexperiments"), type = "logical", default = FALSE,
              help = "Return experiment directly, performs only 1st step of the merging", metavar = "logical"),
  make_option(c("-b", "--bccrsass"), type = "character", default = NULL,
              help = "CRS-BC association file", metavar = "character"),
  make_option(c("-q", "--metainfo"), type = "character", default = NULL,
              help = "Meta info Library", metavar = "character"),
  make_option(c("-m", "--metacontrolsinfo"), type = "character", default = NULL,
              help = "Meta controls info Library", metavar = "character"),  
  make_option(c("-n", "--normalizetoinput"), type = "logical", default = FALSE,
              help = "Normalize to input", metavar = "logical"),
  make_option(c("-p", "--pseudocount"), type = "logical", default = FALSE,
              help = "Pseudocount", metavar = "logical"),
  make_option(c("-t", "--pseudocounttype"), type = "character", default = "DNA",
              help = "Type of pseudocount: DNA, RNA or BOTH", metavar = "character"),
  make_option(c("-v", "--pseudocountvalue"), type = "integer", default = 1,
              help = "Pseudocount value", metavar = "integer"),
  make_option(c("-u", "--usecol"), type = "character", default = "DNA",
              help = "Type of pseudocount: DNA, RNA or BOTH", metavar = "character"),
  make_option(c("-x", "--mindnaumi"), type = "integer", default = 1,
              help = "Only include CRS BC covered with at least this many UMIs at DNA level", metavar = "integer"),              
  make_option(c("-y", "--minrnaumi"), type = "integer", default = 0,
              help = "Only include CRS BC covered with at least this many UMIs at RNA level", metavar = "integer"),     
  make_option(c("-d", "--dnacrsisumi"), type = "character", default = "Binarize",
              help = "Whether to treat DNA CRS observations as single molecules no matter how many UMIs there are", metavar = "character"),
  make_option(c("-z", "--rnacrsisumi"), type = "logical", default = FALSE,
              help = "Whether to treat RNA CRS observations as single molecules no matter how many UMIs there are", metavar = "logical"),
  make_option(c("-j", "--medupper"), type = "integer", default = 2,
              help = "If dnacrsisumi=Median select upper value for counts (UMI >= median(UMI))", metavar = "integer"),     
  make_option(c("-k", "--medlower"), type = "integer", default = 1,
              help = "If dnacrsisumi=Median select lower value for counts (UMI >= median(UMI))", metavar = "integer"),   
  make_option(c("-g", "--crsfilter"), type = "integer", default = 20,
              help = "min UMIs per CRS to be flagged as pass", metavar = "integer"),   
  make_option(c("-w", "--barcodeWhitelist"), type = "logical", default = FALSE,
              help = "Barcode Whitelist", metavar = "logical"),  
  make_option(c("-a", "--AssReadsFilter"), type = "integer", default = 5,
              help = "Number of minimum reads for a valid BC in the association file", metavar = "integer"),   
  make_option(c("-h", "--homopolymer"), type = "integer", default = 0,
              help = "max stretch of homopolymer allowed in barcode. If 0, no filter is applied", metavar = "integer"),   
  make_option(c("-f", "--filtercontrolGRE"), type = "logical", default = FALSE,
              help = "Filter the control barcodes out", metavar = "logical"), 
  make_option(c("-l", "--plot"), type = "logical", default = TRUE,
              help = "Make plots", metavar = "logical"),
  make_option(c("-r", "--returnreport"), type = "logical", default = TRUE,
              help = "Besides returning the CRS data frame also return a list with summary statistics", metavar = "logical"),
  make_option(c("-s", "--savereport"), type = "logical", default = TRUE,
              help = "Save the summary stats in a CSV file", metavar = "logical"),
  make_option(c("--filtercolumns"), type = "logical", default = TRUE,
              help = "Filter the data.frame for use.col to save memory in merging", metavar = "logical"),
  make_option(c("--trp53"), type = "character", default = NULL,
              help = "Trp53 Library", metavar = "character")         
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
required_args <- c("input", "usecol", "output")
for (i in required_args) {
  check_required_args(i, opt = opt, opt_parser = opt_parser) 
}

# Running the script
source("/path/to/project/folder/MPRAscripts/3_Report_screen/scripts/get_CRS_summary.R")
library(parallel)

#Select folder
folders <- list.files(opt$input,full.names = T,pattern = "^State")
names(folders) <- list.files(opt$input,full.names = F,pattern = "^State")

experiments <- readRDS(opt$infile)

all.bc  <- mcmapply(get_CRS_summary, folders, names(folders), experiments,
                      MoreArgs = list(use.col = opt$usecol, 
                                      return_report = opt$returnreport,
                                      return_experiments = opt$returnexperiments,
                                      dna.crs.is.umi = opt$dnacrsisumi, 
                                      rna.crs.is.umi= opt$rnacrsisumi,
                                      med.upper = opt$medupper,
                                      med.lower = opt$medlower,
                                      crs_filter = opt$crsfilter,
                                      filter_controlGRE = opt$filtercontrolGRE,
                                      save_report= opt$savereport,
                                      min.dna.umi = opt$mindnaumi,
                                      min.rna.umi = opt$minrnaumi,
                                      normalize.to.input = opt$normalizetoinput,
                                      pseudocount = opt$pseudocount,
                                      pseudocount_type = opt$pseudocounttype,
                                      pseudocount_value = opt$pseudocountvalue,
                                      filter_columns = opt$filtercolumns,
                                      meta_info =opt$metainfo,
                                      meta_controls_info = opt$metacontrolsinfo,
                                      bc_crs_ass= opt$bccrsass,
                                      output = opt$output,
                                      barcodeWhitelist= opt$barcodeWhitelist,
                                      AssReadsFilter = opt$AssReadsFilter,
                                      plot = opt$plot,
                                      homopolymer = opt$homopolymer,
                                      trp53 = opt$trp53), 
                      SIMPLIFY = F , mc.cores = opt$cores)


reports.umi <- lapply(names(all.bc), function(xx) if (is.list(all.bc[[xx]])) all.bc[[xx]][[2]] else NULL)
names(reports.umi) <- names(all.bc)

final.crs <- lapply(names(all.bc)[sapply(all.bc, is.list)], function(x) {
  out <- all.bc[[x]][[1]]
  out$clusterID <- x
return(out)
})

final.crs <- do.call(rbind, final.crs)
save(final.crs, reports.umi, file = paste0(opt$output, "/final.crs_report.umi.rda"))
