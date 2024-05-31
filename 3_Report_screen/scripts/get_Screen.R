
'%nin%' = Negate('%in%')

#Function
get_Screen <- function(f,
                        bc_crs_ass = NULL,
                        normalize.to.input= FALSE
                        )
{
  suppressPackageStartupMessages({
    require(plyr)
    require(dplyr)
    require(ggplot2)
    require(reshape2)
    require(readr)
    require(data.table)
  })
     
  #read in data (skipped if data is provided)
  rna <- list.files(f, "*RNA.+.map.csv.gz")
  RNA <- lapply(file.path(f, rna), function(x) read_csv(gzfile(x), col_types = cols()))
  if(normalize.to.input) {
    # Subset RNA for only BARCODE, CRS & UMI 
    RNA <- lapply(RNA, setDT)
    # Read in association
    BC_CRS <- read_csv(bc_crs_ass, col_types = cols())
    # Subset BC_CRS for only BARCODE, CRS & READS
    DNA <- BC_CRS[, c("BARCODE", "CRS", "READS")]
    DNA <- setDT(DNA)
    # Setkey for data table
    DNA <- setkey(DNA, BARCODE, CRS)
    RNA <- lapply(RNA, setkey, BARCODE, CRS)
    # Merge Samples 
    Experiments <- lapply(RNA, function(x) x[DNA, on = .(BARCODE, CRS), mult = "all"])
    Experiments <- lapply(Experiments, function(x) {colnames(x)[c(ncol(x))] <- c("READS.dna"); x})
    Experiments <- lapply(Experiments, function(x) {colnames(x)[c(3:(ncol(x)-1))] <- paste(colnames(x)[c(3:(ncol(x)-1))],".rna", sep = ""); x})
  }
  if(!normalize.to.input) {
    #Read in DNA flies
    dna <- list.files(f, "*DNA.+.map.csv.gz")
    DNA <- lapply(file.path(f, dna), function(x) read_csv(gzfile(x), col_types = cols()))
    DNA <- lapply(DNA, setDT)  
    DNA <- lapply(DNA, setkey, BARCODE, CRS)
    #Use RNA
    RNA <- lapply(RNA, setDT)     
    RNA <- lapply(RNA, setkey, BARCODE, CRS)
    # Merge Samples
    Experiments <- mapply(merge.data.table,RNA, DNA, MoreArgs = list(by = c("BARCODE","CRS"), suffixes = c(".rna",".dna"), all = T),SIMPLIFY = FALSE )
    }
  return(Experiments)
}
