#systematically create CRS level summaries from the perl pipeline output

#example 1 - process several folders in parallel:
# folders <- c(EoBaso = "/users/lvelten/lvelten/Analysis/SCG4SYN/LibA_HSC/raw/out/EoBaso/",
#              MonoPre = "/users/lvelten/lvelten/Analysis/SCG4SYN/LibA_HSC/raw/out/MonoP/",
#              EoBaso500 = "/users/lvelten/lvelten/Analysis/SCG4SYN/LibA_HSC/raw/out/EoBaso500/",
#              EoBasoDS = "/users/lvelten/lvelten/Analysis/SCG4SYN/LibA_HSC/raw/out/EoBaso_ds//")
# CRS.all <- mcmapply(get_CRS_summary,folders, names(folders), MoreArgs = list(outdir = "out",use.col = "UMI5"), SIMPLIFY = F, mc.cores = 3)

#example 2 - best conditions so far
#CRS.all <- get_CRS_summary("/users/lvelten/lvelten/Analysis/SCG4SYN/LibA_HSC/raw/out/EoBaso/", "EoBaso_conf5_dnaisumi",
#outdir = "out", use.col = "UMI5", return_report = T, 
#dna.crs.is.umi = T, crs_filter = 20)

get_CRS_summary <- function(f, run_id, Experiments = NULL,
         use.col = "UMI5", #the column in the PERL pipeline output to use (e.g. UMI5: UMIs with >=5 reads)
         plot =T, #create pllots on the way in the directory specified as outdir
         outdir = "processed", 
         min.dna.umi = 1, #to filter at the level of CRS barcodes: Only include CRS BC covered with at least this many UMIs at DNA level
         min.rna.umi = 0, #to filter at the level of CRS barcodes: Only include CRS BC covered with at least this many UMIs at RNA level
         umi.filter.connection = `&`, #DNA filter AND RNA filter or DNA filter OR RNA filter?
         dna.crs.is.umi = F, #whether to treat DNA CRS observations as single molecules no matter how many UMIs there are
         rna.crs.is.umi = F,
         bc_crs_ass = "/users/lvelten/lvelten/Analysis/SCG4SYN/LibA/CRS_BC_ass/myPipeline/outNextseq/mapped.filtered.csv",
         crs_filter =100, #min UMIs per CRS to be flagged as pass
         return_report = F, #besides returning the CRS data frame also return a list with some summary stats?
         return_experiments=F #if set to T, only runs the first steps of loading and merging data. The return value can be passed to future function calls via the Experiments parameter
)
{
  require(plyr)
  report <- list()
  dir.create(file.path(outdir, run_id), showWarnings = T, recursive = T)
  
  #read in data (skipped if data is provided)
  if (is.null(Experiments)){
    rna <- list.files(f, "RNA\\d.csv")
    dna <- list.files(f, "DNA\\d.csv")
    RNA <- lapply(file.path(f, rna), read.csv)
    DNA <- lapply(file.path(f, dna), read.csv)
    Experiments <- mapply(merge, DNA, RNA, MoreArgs = list( by = c("BARCODE","CRS"), suffixes = c(".dna",".rna"), all = T), SIMPLIFY =FALSE )
  }
  
  if (return_experiments) return(Experiments)
  
  #replace NAs by 0
  Experiments <- lapply(Experiments, function(exx) {
    for (i in 3:ncol(exx)) exx[is.na(exx[,i]),i] <- 0
    exx$norm <- exx[,paste0(use.col,".rna")] / exx[,paste0(use.col,".dna")]
    exx$USE.dna <- exx[,paste0(use.col,".dna")]
    exx$USE.rna <- exx[,paste0(use.col,".rna")]
    exx
  })
  
  #Create CRS-BC level plots
  if(plot) {
    all  <-merge(Experiments[[1]], Experiments[[2]],by = c("BARCODE","CRS"), suffixes = c(".1",".2"), all=T)
    all.filtered <- subset(all, UMI.dna.1 > 0 | UMI.dna.2 > 0 )
    toplot<-sample(1:nrow(all.filtered), 50000)
    
    png(file.path(outdir, run_id, "BC_replicates_%02d.png"),width=800,height=600,res=100)
    heatscatter(log10(all.filtered$UMI.dna.1[toplot]+1),log10( all.filtered$UMI.rna.1[toplot]+1), main = "DNA vs RNA Exp 1", cor=T, method = "pearson")
    heatscatter(log10(all.filtered$UMI.rna.1[toplot]+1),log10( all.filtered$UMI.rna.2[toplot]+1), main = "RNA vs RNA", cor=T, method = "pearson")
    heatscatter(log10(all.filtered$UMI.dna.1[toplot]+1),log10( all.filtered$UMI.dna.2[toplot]+1), main = "DNA vs DNA", cor=T, method = "pearson")
    dev.off()
  }
  
  #reporting
  report$nUMI.RNA <- sapply(Experiments, function(x) sum(x$USE.rna))
  report$nUMI.DNA <- sapply(Experiments, function(x) sum(x$USE.dna))
  report$crs_bc_covered.RNA <- sapply(Experiments, function(x) sum(x$USE.rna>0))
  report$crs_bc_covered.DNA <- sapply(Experiments, function(x) sum(x$USE.dna>0))
  
  #filtering of CRS barcodes for merging by CRS
  BC_CRS <- read.csv(bc_crs_ass)
  selected_barcodes <- subset(BC_CRS, MEANMATCHES >= 290 & DEVIANTREADS / (READS + DEVIANTREADS) < 0.2)
  
  Experiments.filtered <- lapply(Experiments, function(x) 
    subset(x, umi.filter.connection(UMI.dna >= min.dna.umi, UMI.rna >= min.rna.umi) & 
             BARCODE %in% selected_barcodes$BARCODE))
  
  report$crs_bc_filter.RNA <- sapply(Experiments.filtered, function(x) sum(x$USE.rna >0))
  report$crs_bc_filter.DNA <- sapply(Experiments.filtered, function(x) sum(x$USE.dna >0))
  
  
  #summarise CRS barcodes at the level of CRS
  CRS <- lapply(Experiments.filtered, function(exx) {
    if (!dna.crs.is.umi & !rna.crs.is.umi) crs <- ddply(exx, "CRS", summarise, RNA = sum(USE.rna), DNA = sum(USE.dna))
    if (dna.crs.is.umi & !rna.crs.is.umi) crs <- ddply(exx, "CRS", summarise, RNA = sum(USE.rna), DNA = sum(USE.dna > 0))
    if (!dna.crs.is.umi & rna.crs.is.umi) crs <- ddply(exx, "CRS", summarise, RNA = sum(USE.rna > 0), DNA = sum(USE.dna))
    if (dna.crs.is.umi & rna.crs.is.umi) crs <- ddply(exx, "CRS", summarise, RNA = sum(USE.rna > 0), DNA = sum(USE.dna > 0))
    
    crs$RNA.norm <- crs$RNA / (sum(crs$RNA) / 1e6)
    crs$DNA.norm <- crs$DNA / (sum(crs$DNA) / 1e6)
    crs$norm <- log2(crs$RNA.norm / crs$DNA.norm)
    crs
  })
  
  #merge and produce final plots
  CRS.all <- merge(CRS[[1]], CRS[[2]], by = "CRS", suffixes = c(".1",".2"), all=T)
  CRS.all$pass <- with(CRS.all, DNA.1 > crs_filter &DNA.2 > crs_filter)
  CRS.all.filtered <- subset(CRS.all, pass )
  
  if (plot) {
    png(file.path(outdir, run_id, "CRS_full_replicates_%02d.png"),width=800,height=600,res=100)
    heatscatter(CRS.all$norm.1,CRS.all$norm.2, main = "normalized", cor=T, method = "pearson")
    heatscatter(log2(CRS.all$RNA.1+1),log2( CRS.all$RNA.2+1), main = "normalized", cor=T, method = "pearson")
    heatscatter(log2(CRS.all$DNA.1+1),log2( CRS.all$DNA.2+1), main = "normalized", cor=T, method = "pearson")
    dev.off()
    
    
    png(file.path(outdir, run_id, "CRS_filtered_replicates_%02d.png"),width=800,height=600,res=100)
    heatscatter(CRS.all.filtered$norm.1,CRS.all.filtered$norm.2, main = "normalized", cor=T, method = "pearson")
    heatscatter(log2(CRS.all.filtered$RNA.1+1),log2(CRS.all.filtered$RNA.2+1), main = "normalized", cor=T, method = "pearson")
    heatscatter(log2(CRS.all.filtered$DNA.1+1),log2(CRS.all.filtered$DNA.2+1), main = "normalized", cor=T, method = "pearson")
    dev.off()
  }
  
  y <- subset(CRS.all, !is.infinite(norm.1) & !is.infinite(norm.2))
  report$norm_cor <- cor(y$norm.1,y$norm.2, use="pairwise.complete.obs")
  report$norm_n <- sum(!is.na(y$norm.1) & !is.na(y$norm.2))
  y <- subset(y, pass)
  report$norm_cor_filter <- cor(y$norm.1,y$norm.2, use="pairwise.complete.obs")
  report$norm_n_filter <- sum(!is.na(y$norm.1) & !is.na(y$norm.2))
  
  if (return_report) return(list(CRS.all, report)) else return(CRS.all)
}