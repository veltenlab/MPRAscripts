
'%nin%' = Negate('%in%')

#Function
get_CRS_summary <- function(f, run_id, Experiments = NULL,
         filter_columns = T, 
         normalize.to.input= F,
         pseudocount = T,
         pseudocount_type = "DNA", 
         pseudocount_value = 1, 
         use.col = "UMI1", #the column in the PERL pipeline output to use (e.g. UMI5: UMIs with >=5 reads)
         plot =T, #create plots on the way in the directory specified as outdir
         outdir = "~/cluster_mount/project/SCG4SYN/Pipelines/MPRA/3_Report_screen/output/4_HSC_LibraryH_minP/filtered/", 
         min.dna.umi = 1, #to filter at the level of CRS barcodes: Only include CRS BC covered with at least this many UMIs at DNA level
         min.rna.umi = 0, #to filter at the level of CRS barcodes: Only include CRS BC covered with at least this many UMIs at RNA level
         umi.filter.connection = `&`, #DNA filter AND RNA filter or DNA filter OR RNA filter?
         dna.crs.is.umi = F, #whether to treat DNA CRS observations as single molecules no matter how many UMIs there are
         rna.crs.is.umi = F,
         barcodeWhitelist = NULL,
         AssReadsFilter = 5,  # Number of minimum reads for a valid BC in the association file
         homopolymer = 0, #max stretch of homopolymer allowed in barcode. If 0, no filter is applied
         bc_crs_ass = "~/cluster_mount/project/SCG4SYN/Pipelines/MPRA/1_CRS_BC_association/output/LibraryH_minP/LibHSeqDesign_LibHminP_mapped_filtered.csv.gz",
         crs_filter =100, #min UMIs per CRS to be flagged as pass
         return_report = F, #besides returning the CRS data frame also return a list with some summary stats?
         save_report = F , # Save the summary stats in a CSV file for later comparisons
         filter_controlGRE = F, #Filter the control barcodes out
         return_experiments=F #if set to T, only runs the first steps of loading and merging data. The return value can be passed to future function calls via the Experiments parameter
)
{
  require(plyr)
  require(LSD)
  require(dplyr)
  require(ggplot2)
  require(reshape2)
  require(readr)
  require(data.table)
  report <- list()
  
  if (!dir.exists(file.path(outdir, run_id))) {dir.create(file.path(outdir, run_id), showWarnings = T, recursive = T )}
  
  #  Setup pattern for pseudocount
  if(pseudocount_type %in% "DNA") pattern_pseudo <- ".dna"
  if(pseudocount_type %in% "RNA") pattern_pseudo <- ".rna"
  if(pseudocount_type %in% "BOTH") pattern_pseudo <- c(".dna", ".rna")
  
  #read in data (skipped if data is provided)

  if (is.null(Experiments)) {
    rna <- list.files(f, "*RNA.+.map.csv.gz")
    RNA <- lapply(file.path(f, rna), function(x) read_csv(gzfile(x)))
    if(normalize.to.input) {
      if (use.col != "READS") {
        print("Select use.col = READS")
      } else {
      # Subset RNA for only BARCODE, CRS & UMI 
      # Unlike for DNA we do not select the reads, but the Umi1 column in the df
      #RNA <- lapply(RNA, function(x) {x[, c("BARCODE", "CRS", use.col)]})
      RNA <- lapply(RNA, function(x) {x[, c("BARCODE", "CRS", "UMI")]})
      RNA <- lapply(RNA, setDT)
      # Read in association
      BC_CRS <- read_csv(bc_crs_ass)
      # Subset BC_CRS for only BARCODE, CRS & READS
      DNA <- BC_CRS[, c("BARCODE", "CRS", use.col)]
      DNA <- setDT(DNA)
      # Setkey for data table
      DNA <- setkey(DNA, BARCODE, CRS)
      RNA <- lapply(RNA, setkey, BARCODE, CRS)
      # Merge Samples
      #This runs faster, but must be subsetted to use.cols beforehand
      Experiments <- lapply(RNA, function(x) x[DNA, on = .(BARCODE, CRS), mult = "all"])
      Experiments <- lapply(Experiments, function(x) {colnames(x)[c(3,4)] <- c("READS.rna", "READS.dna"); x})
      #general but slower option
      #Experiments <- lapply(RNA, function(x) {merge.data.table(x, DNA, by = c("BARCODE","CRS"), suffixes = c(".rna",".dna"), all = T)}) 
      }
      }
    if(!normalize.to.input) {
      #Read in DNA flies
      dna <- list.files(f, "*DNA.+.map.csv.gz")
      DNA <- lapply(file.path(f, dna), function(x) read_csv(gzfile(x)))
      if(filter_columns) {
      DNA <- lapply(DNA, function(x) {x[, c("BARCODE", "CRS", use.col)]})
      }
      DNA <- lapply(DNA, setDT)  
      DNA <- lapply(DNA, setkey, BARCODE, CRS)
      #Use RNA
      if(filter_columns) {
      RNA <- lapply(RNA, function(x) {x[, c("BARCODE", "CRS", use.col)]})
      }
      RNA <- lapply(RNA, setDT)      
      RNA <- lapply(RNA, setkey, BARCODE, CRS)
      # Merge Samples
      Experiments <- mapply(merge.data.table,RNA, DNA, MoreArgs = list(by = c("BARCODE","CRS"), suffixes = c(".rna",".dna"), all = T),SIMPLIFY = FALSE )
      if(pseudocount) {
        if(length(pattern_pseudo) == 1) {
          #this works
 #         Experiments2 <- Experiments
 #         for (i in seq_along(Experiments2)) {
 #           # Get the current data.table
 #           dt <- Experiments2[[i]]
 #           # Find columns matching the pattern and perform the required operations
 #           matching_cols <- grep(pattern = pattern_pseudo, x = colnames(dt), value = TRUE)
 #           for (col in matching_cols) {
 #             dt[[col]] <- ifelse(is.na(dt[[col]]), pseudocount_value, dt[[col]])
 #           }
 #           # Store the modified data.table back to the list
 #           Experiments2[[i]] <- dt
 #         }
          #Short form with lapply
          Experiments <- lapply(Experiments, function(x) {matching_cols <- grep(pattern = pattern_pseudo, x = colnames(x), value = TRUE); 
          for (col in matching_cols) { x[[col]] <- ifelse(is.na(x[[col]]), pseudocount_value, x[[col]])}; x})
          #x[[matching_cols]] <- ifelse(is.na(x[[matching_cols]]), pseudocount_value, x[[matching_cols]] ); x}) # for 1 entry in matching_cols         
          
         # Experiments <- lapply(Experiments, function(x) {x[,grepl(pattern = pattern_pseudo, x = colnames(x))] <- x[,grepl(pattern = pattern_pseudo, x = colnames(x))] %>% 
         #   mutate(any_na = rowSums(is.na(.)) > 0) %>%
         #   mutate(across(everything(), ~ ifelse(any_na, pseudocount_value, .))) %>%
         #   select(-any_na); x})
        } else {
          Experiments <- lapply(Experiments, function(x) {matching_cols <- grep(pattern = pattern_pseudo[1], x = colnames(x), value = TRUE); 
          for (col in matching_cols) { x[[col]] <- ifelse(is.na(x[[col]]), pseudocount_value, x[[col]])};
          matching_cols <- grep(pattern = pattern_pseudo[2], x = colnames(x), value = TRUE); 
          for (col in matching_cols) { x[[col]] <- ifelse(is.na(x[[col]]), pseudocount_value, x[[col]])};x})
          
    #      Experiments <- lapply(Experiments, function(x) {x[,grepl(pattern = pattern_pseudo[1], x = colnames(x))] <- x[,grepl(pattern = pattern_pseudo[1], x = colnames(x))] %>% 
    #        mutate(any_na = rowSums(is.na(.)) > 0) %>%
    #        mutate(across(everything(), ~ ifelse(any_na, pseudocount_value, .))) %>%
    #        select(-any_na); x[,grepl(pattern = pattern_pseudo[2], x = colnames(x))] <- x[,grepl(pattern = pattern_pseudo[2], x = colnames(x))] %>% 
    #          mutate(any_na = rowSums(is.na(.)) > 0) %>%
    #          mutate(across(everything(), ~ ifelse(any_na, pseudocount_value, .))) %>%
    #          select(-any_na); x})
        }
      }
    }
    rm(RNA, DNA); gc()
  }
  
  if (return_experiments) return(Experiments)

#  Experiments <- experiments[1]  

  if(filter_columns) {
    Experiments <- lapply(Experiments, function(x) {
      cbind(x[, "BARCODE", with = FALSE], x[, "CRS", with = FALSE], x[, "UMI.rna", with = FALSE], x[, paste0(use.col, ".rna"), with = FALSE], x[, "UMI.dna", with = FALSE], x[, paste0(use.col, ".dna"), with = FALSE])
      })
    gc()
  }
  
  if (!return_experiments) {
    if(dna.crs.is.umi == "Median") {
      med <- sapply(Experiments, function(x) {
        temp <- ifelse(is.na(x[[paste0(use.col, ".dna")]]), 0, x[[paste0(use.col, ".dna")]])
        temp <- temp[temp>0]
        median(temp)
        })
    }
    if(pseudocount) {
      if(length(pattern_pseudo) == 1) {
        Experiments <- lapply(Experiments, function(x) {matching_cols <- grep(pattern = pattern_pseudo, x = colnames(x), value = TRUE); 
        for (col in matching_cols) { x[[col]] <- ifelse(is.na(x[[col]]), pseudocount_value, x[[col]])}; x})
      } else {
        Experiments <- lapply(Experiments, function(x) {matching_cols <- grep(pattern = pattern_pseudo[1], x = colnames(x), value = TRUE); 
        for (col in matching_cols) { x[[col]] <- ifelse(is.na(x[[col]]), pseudocount_value, x[[col]])};
        matching_cols <- grep(pattern = pattern_pseudo[2], x = colnames(x), value = TRUE); 
        for (col in matching_cols) { x[[col]] <- ifelse(is.na(x[[col]]), pseudocount_value, x[[col]])};x})
      }
    }
  }
  
  if (homopolymer > 0) {
    filter <- sapply(c("A","C","G","T"), function(x) paste(rep(x, homopolymer), collapse = "") )
    filter <- paste(filter, collapse = "|")
    Experiments <- lapply(Experiments, function(exx) subset(exx, !grepl(filter, BARCODE)))
  }
  
  if (is.null(barcodeWhitelist)) {
    if (!exists("BC_CRS")) {
      BC_CRS <- read_csv(bc_crs_ass)
    }
    barcodeWhitelist <- subset(BC_CRS, MEANMATCHES >= 290 & DEVIANTREADS / (READS + DEVIANTREADS) < 0.2 & READS >= AssReadsFilter)$BARCODE
  }
  
  Experiments <- lapply(Experiments, function(exx) subset(exx, BARCODE %in% barcodeWhitelist))
  rm(barcodeWhitelist)
  
  #replace NAs by 0
  #
  Experiments <- lapply(Experiments, function(exx) {
    for (i in 3:ncol(exx)) { exx[[i]] <- ifelse(is.na(exx[[i]]), 0, exx[[i]]) }
    exx$norm <- exx[[paste0(use.col, ".rna")]] / exx[[paste0(use.col, ".dna")]]
    exx$USE.dna <- exx[[paste0(use.col, ".dna")]]
    exx$USE.rna <- exx[[paste0(use.col, ".rna")]]
   # for (i in 3:ncol(exx)) exx[is.na(exx[,i]),i] <- 0
   # exx$norm <- exx[,paste0(use.col,".rna")] / exx[,paste0(use.col,".dna")]
   # exx$USE.dna <- exx[,paste0(use.col,".dna")]
   # exx$USE.rna <- exx[,paste0(use.col,".rna")]
    exx
  })
  

  #Experiments <- backup
  #backup <- Experiments

  #Filter control GREs
  if(filter_controlGRE) {
    Experiments <- lapply(Experiments, function(exx){
      exx <- exx[which(!grepl(pattern = "^Ctrl_*", x = exx$CRS, ignore.case = T)),]
    })
  }
  
  #reporting
  report$nUMI.RNA <- sapply(Experiments, function(x) sum(x$USE.rna))
  report$nUMI.DNA <- sapply(Experiments, function(x) sum(x$USE.dna))
  report$crs_bc_covered.RNA <- sapply(Experiments, function(x) sum(x$USE.rna>0))
  report$crs_bc_covered.DNA <- sapply(Experiments, function(x) sum(x$USE.dna>0))
  
  #filtering of CRS barcodes for merging by CRS
  Experiments.filtered <- lapply(Experiments, function(x) 
    subset(x, umi.filter.connection(UMI.dna >= min.dna.umi, UMI.rna >= min.rna.umi)))
  
  report$crs_bc_filter.RNA <- sapply(Experiments.filtered, function(x) sum(x$USE.rna >0))
  report$crs_bc_filter.DNA <- sapply(Experiments.filtered, function(x) sum(x$USE.dna >0))
  
  #summarise CRS barcodes at the level of CRS
  CRS <- lapply(Experiments.filtered, function(exx) {
    
    #if (dna.crs.is.umi == "Median") crs <- ddply(exx, "CRS", summarise, RNA = sum(USE.rna > 0), DNA = sum(ifelse(USE.dna == 0, 0, ifelse(USE.dna > med, 2, 1)))
    # med <- median(exx$USE.dna)
    if(dna.crs.is.umi == "Median" & !rna.crs.is.umi & pseudocount) {
      return(med)
      #crs <- ddply(exx, "CRS", summarise, RNA = sum(USE.rna), DNA = sum(ifelse(USE.dna == 1, 1, ifelse(USE.dna >= med, 3, 2))))
      }
                                                      
    if (!dna.crs.is.umi & !rna.crs.is.umi) crs <- ddply(exx, "CRS", summarise, RNA = sum(USE.rna), DNA = sum(USE.dna))
    if (dna.crs.is.umi & !rna.crs.is.umi) crs <- ddply(exx, "CRS", summarise, RNA = sum(USE.rna), DNA = sum(USE.dna > 0))
    if (!dna.crs.is.umi & rna.crs.is.umi) crs <- ddply(exx, "CRS", summarise, RNA = sum(USE.rna > 0), DNA = sum(USE.dna))
    if (dna.crs.is.umi & rna.crs.is.umi) crs <- ddply(exx, "CRS", summarise, RNA = sum(USE.rna > 0), DNA = sum(USE.dna > 0))
    
    crs$RNA.norm <- crs$RNA / (sum(crs$RNA) / 1e6)
    crs$DNA.norm <- crs$DNA / (sum(crs$DNA) / 1e6)
    crs$norm <- log2((crs$RNA.norm +1)/ (crs$DNA.norm+1))
    crs
  })
  
  #remove unnecessary files
  rm(Experiments, Experiments.filtered); gc()
  
  #merge and produce final plots
  CRS.all <- merge(CRS[[1]], CRS[[2]], by = "CRS", suffixes = c(".1",".2"), all=T)
  CRS.all$pass <- with(CRS.all, DNA.1 > crs_filter &DNA.2 > crs_filter)

  CRS.all <- CRS.all %>% 
    dplyr::filter(norm.1 %nin% Inf & norm.1 %nin% -Inf & norm.1 %nin% NA & norm.1 %nin% NaN & norm.2 %nin% Inf & norm.2 %nin% -Inf & norm.2 %nin% NA & norm.2 %nin% NaN)
  
  CRS.all <- CRS.all %>%
    rowwise() %>%
    dplyr::mutate(mean.norm = mean(c(norm.1, norm.2))) %>%
    ungroup()
  
  #Create filtered version
  CRS.all.filtered <- subset(CRS.all, pass )

  if (plot) {
    png(file.path(outdir, run_id, "CRS_full_replicates_%02d.png"),width=800,height=600,res=100)
    heatscatter(CRS.all$norm.1,CRS.all$norm.2, main = "normalized", cor=T, method = "pearson")
    heatscatter(log2(CRS.all$RNA.1+1),log2( CRS.all$RNA.2+1), main = "normalized", cor=T, method = "pearson")
    heatscatter(log2(CRS.all$DNA.1+1),log2( CRS.all$DNA.2+1), main = "normalized", cor=T, method = "pearson")
    dev.off()
    
    #Filtered CRSs
    png(file.path(outdir, run_id, "CRS_filtered_replicates_%02d.png"),width=800,height=600,res=100)
    heatscatter(CRS.all.filtered$norm.1,CRS.all.filtered$norm.2, main = "normalized", cor=T, method = "pearson")
    heatscatter(log2(CRS.all.filtered$RNA.1+1),log2(CRS.all.filtered$RNA.2+1), main = "normalized", cor=T, method = "pearson")
    heatscatter(log2(CRS.all.filtered$DNA.1+1),log2(CRS.all.filtered$DNA.2+1), main = "normalized", cor=T, method = "pearson")
    dev.off()
}
  
  y <- subset(CRS.all, !is.infinite(norm.1) & !is.infinite(norm.2))
  report$norm_cor <- cor(y$norm.1,y$norm.2, use="pairwise.complete.obs")
  report$RNA_cor <- cor(log2(y$RNA.1+1),log2(y$RNA.2+1), use="pairwise.complete.obs")
  report$DNA_cor <- cor(log2(y$DNA.1+1),log2(y$DNA.2+1), use="pairwise.complete.obs")
  report$norm_n <- sum(!is.na(y$norm.1) & !is.na(y$norm.2))
  y <- subset(y, pass)
  report$norm_cor_filter <- cor(y$norm.1,y$norm.2, use="pairwise.complete.obs")
  report$norm_n_filter <- sum(!is.na(y$norm.1) & !is.na(y$norm.2))
  if(save_report) { 
    report_safe <- as.data.frame(t(as.data.frame(report)))
    colnames(report_safe) <- paste(run_id, c("__Replicate1","__Replicate2"), sep = "")
    save(report_safe, CRS.all, file = file.path(outdir, run_id, "CRS_summary_report.rda"))
    }
  if (return_report) return(list(CRS.all, report)) else return(CRS.all)
  gc()
}
