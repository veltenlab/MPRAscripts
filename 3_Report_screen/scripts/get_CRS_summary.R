
#trp53 cor all and cell states, median calc, output and wrapper


'%nin%' = Negate('%in%')

#Function
get_CRS_summary <- function(f, run_id, Experiments = NULL,
        #General
          return_experiments=F, #if set to T, only runs the first steps of loading and merging data. The return value can be passed to future function calls via the Experiments parameter
        #Additional files
          bc_crs_ass = "/path/to/project/folder/MPRAscripts/1_CRS_BC_association/output_snake/LibraryA_minPInitial/LibASeqDesign_LibAminPInitial_mapped_filtered.csv.gz",
          meta_info = NULL,
          meta_controls_info = NULL,
        #Filter section
         #Subset columns
          filter_columns = T, 
         #Input
          normalize.to.input= F,
         #Pseudocount
          pseudocount = T,
          pseudocount_type = "DNA", 
          pseudocount_value = 1, 
         #UMI behaviour
          use.col = "UMI1", #the column in the PERL pipeline output to use (e.g. UMI5: UMIs with >=5 reads)      
          min.dna.umi = 1, #to filter at the level of CRS barcodes: Only include CRS BC covered with at least this many UMIs at DNA level
          min.rna.umi = 0, #to filter at the level of CRS barcodes: Only include CRS BC covered with at least this many UMIs at RNA level
          umi.filter.connection = `&`, #DNA filter AND RNA filter or DNA filter OR RNA filter?
          dna.crs.is.umi = F, #whether to treat DNA CRS observations as single molecules no matter how many UMIs there are
          rna.crs.is.umi = F,
          med.upper = 2, 
          med.lower = 1,
          crs_filter = 40, #min UMIs per CRS to be flagged as pass
         #Association
          barcodeWhitelist = FALSE,
          AssReadsFilter = 3,  # Number of minimum reads for a valid BC in the association file
         #Other
          homopolymer = 0, #max stretch of homopolymer allowed in barcode. If 0, no filter is applied
          filter_controlGRE = F, #Filter the control barcodes out
        #output section
          output = "/path/to/project/folder/MPRAscripts/3_Report_screen/output_snake/2_HSC_LibraryA_minPInitial/filtered",         
          plot =T, #create plots on the way in the directory specified as outdir         
          return_report = F, #besides returning the CRS data frame also return a list with some summary stats?
          save_report = F, # Save the summary stats in a CSV file for later comparisons
          trp53 = "/path/to/project/folder/MPRAscripts/Library_design/meta_files/LibASeqDesign_meta_info.csv"
)
{
  suppressPackageStartupMessages({
    require(plyr)
    require(LSD)
    require(dplyr)
    require(ggplot2)
    require(reshape2)
    require(readr)
    require(data.table)
  })
  report <- list()

  #Incorporate parameter information into folder name
  outdir <- output
  report$outdir <- outdir
  
  #  Setup pattern for pseudocount
  if(pseudocount_type %in% "DNA") pattern_pseudo <- ".dna"
  if(pseudocount_type %in% "RNA") pattern_pseudo <- ".rna"
  if(pseudocount_type %in% "BOTH") pattern_pseudo <- ".dna|.rna"
   
  #Changed position, so initial folder isnt created
  if(!dir.exists(file.path(outdir, run_id))) {dir.create(file.path(outdir, run_id), showWarnings = T, recursive = T )}

  #read in data (skipped if data is provided)
  if(is.null(Experiments)) {
    rna <- list.files(f, "*RNA.+.map.csv.gz")
    RNA <- lapply(file.path(f, rna), function(x) read_csv(gzfile(x)))
    if(normalize.to.input) {
      # Subset RNA for only BARCODE, CRS & UMI 
      RNA <- lapply(RNA, function(x) {x[, c("BARCODE", "CRS", use.col)]})
      RNA <- lapply(RNA, setDT)
      # Read in association
      BC_CRS <- read_csv(bc_crs_ass)
      # Subset BC_CRS for only BARCODE, CRS & READS
      DNA <- BC_CRS[, c("BARCODE", "CRS", "READS")]
      DNA <- setDT(DNA)
      # Setkey for data table
      DNA <- setkey(DNA, BARCODE, CRS)
      RNA <- lapply(RNA, setkey, BARCODE, CRS)
      # Merge Samples - This runs faster, but must be subsetted to use.cols beforehand
      Experiments <- lapply(RNA, function(x) x[DNA, on = .(BARCODE, CRS), mult = "all"])
      Experiments <- lapply(Experiments, function(x) {colnames(x)[c(3,4)] <- c("READS.rna", "READS.dna"); x})
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
          Experiments <- lapply(Experiments, function(x) {matching_cols <- grep(pattern = pattern_pseudo, x = colnames(x), value = TRUE); 
          for (col in matching_cols) { x[[col]] <- ifelse(is.na(x[[col]]), pseudocount_value, x[[col]])}; x})
        }
      }
    }
    rm(RNA, DNA); gc()
  }
  
  if (return_experiments) return(Experiments)

#Reading in the experiment .rds and modify it accordingly
  if(!is.null(Experiments)) {
    if(normalize.to.input) {
      Experiments <- lapply(Experiments, function(x) {x$UMI.dna <- x$READS.dna; x})
      Experiments <- lapply(Experiments, function(x) {colnames(x)[c(ncol(x)-1)] <-  paste(use.col, ".dna", sep = ""); x})
    }
    if(filter_columns) {
      Experiments <- lapply(Experiments, function(x) {
        cbind(x[, "BARCODE", with = FALSE], x[, "CRS", with = FALSE], x[, "UMI.rna", with = FALSE], x[, paste0(use.col, ".rna"), with = FALSE], x[, "UMI.dna", with = FALSE], x[, paste0(use.col, ".dna"), with = FALSE])
        })
      gc()
    }
    if (!normalize.to.input) {
      if(dna.crs.is.umi == "Median") {
        med <- sapply(Experiments, function(x) {
          temp <- ifelse(is.na(x[[paste0(use.col, ".dna")]]), 0, x[[paste0(use.col, ".dna")]])
          temp <- temp[temp>0]
          median(temp)
          })}
    if(pseudocount) {
      if(length(pattern_pseudo) == 1) {
        Experiments <- lapply(Experiments, function(x) {matching_cols <- grep(pattern = pattern_pseudo, x = colnames(x), value = TRUE); 
        for (col in matching_cols) { x[[col]] <- ifelse(is.na(x[[col]]), pseudocount_value, x[[col]])}; x})
      }
    }
  }
  }

  if (homopolymer > 0) {
    filter <- sapply(c("A","C","G","T"), function(x) paste(rep(x, homopolymer), collapse = "") )
    filter <- paste(filter, collapse = "|")
    Experiments <- lapply(Experiments, function(exx) subset(exx, !grepl(filter, BARCODE)))
  }
  
  if (!barcodeWhitelist) {
    if (!exists("BC_CRS")) {
      BC_CRS <- read_csv(bc_crs_ass, col_types = cols())
    }
    barcodeWhitelist <- subset(BC_CRS, DEVIANTREADS / (READS + DEVIANTREADS) < 0.2 & READS >= AssReadsFilter)$BARCODE # MEANMATCHES >= 290 &
  }

  Experiments <- lapply(Experiments, function(exx) subset(exx, BARCODE %in% barcodeWhitelist))
  rm(barcodeWhitelist)
  
  #replace NAs by 0
  Experiments <- lapply(Experiments, function(exx) {
    for (i in 3:ncol(exx)) { exx[[i]] <- ifelse(is.na(exx[[i]]), 0, exx[[i]]) }
    exx$norm <- exx[[paste0(use.col, ".rna")]] / exx[[paste0(use.col, ".dna")]]
    exx$USE.dna <- exx[[paste0(use.col, ".dna")]]
    exx$USE.rna <- exx[[paste0(use.col, ".rna")]]
    exx
  })
  
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
  CRS <- lapply(1:length(Experiments.filtered), function(i) { 
    exx <- Experiments.filtered[[i]]
    if (dna.crs.is.umi == "Median" & !rna.crs.is.umi & pseudocount) crs <- ddply(exx, "CRS", summarise, RNA = sum(USE.rna), DNA = sum(ifelse(USE.dna == 1, 1, ifelse(USE.dna >= med[i], med.upper , med.lower))))
    if (dna.crs.is.umi == "Non-Binarize" & !rna.crs.is.umi) crs <- ddply(exx, "CRS", summarise, RNA = sum(USE.rna), DNA = sum(USE.dna))
    if (dna.crs.is.umi == "Binarize" & !rna.crs.is.umi) crs <- ddply(exx, "CRS", summarise, RNA = sum(USE.rna), DNA = sum(USE.dna > 0))
    if (dna.crs.is.umi == "Non-Binarize" & rna.crs.is.umi) crs <- ddply(exx, "CRS", summarise, RNA = sum(USE.rna > 0), DNA = sum(USE.dna))
    if (dna.crs.is.umi == "Binarize" & rna.crs.is.umi) crs <- ddply(exx, "CRS", summarise, RNA = sum(USE.rna > 0), DNA = sum(USE.dna > 0))
    
    crs$RNA.norm <- crs$RNA / (sum(crs$RNA) / 1e6)
    crs$DNA.norm <- crs$DNA / (sum(crs$DNA) / 1e6)
    crs$norm <- log2((crs$RNA.norm +1)/ (crs$DNA.norm+1))
    crs
  })
  
  #remove unnecessary files
  rm(Experiments, Experiments.filtered); gc()
  
  #merge and produce final plots
  CRS.all <- merge(CRS[[1]], CRS[[2]], by = "CRS", suffixes = c(".1",".2"), all=T)
  CRS.all$pass <- with(CRS.all, DNA.1 > crs_filter & DNA.2 > crs_filter & RNA.1 > crs_filter & RNA.2 > crs_filter ) # RNA and DNA

  CRS.all <- CRS.all %>% 
    dplyr::filter(norm.1 %nin% Inf & norm.1 %nin% -Inf & norm.1 %nin% NA & norm.1 %nin% NaN & norm.2 %nin% Inf & norm.2 %nin% -Inf & norm.2 %nin% NA & norm.2 %nin% NaN) %>%
    rowwise() %>%
    dplyr::mutate(mean.norm = mean(c(norm.1, norm.2))) %>%
    ungroup()
  
  #Create filtered version
  CRS.all.filtered <- subset(CRS.all, pass )
  
  if(!is.null(meta_controls_info)) {
    #Read in files
    LibASeqDesign <- read_csv(trp53, col_types = cols())
    Design_meta_controls <- read_csv(meta_controls_info, col_types = cols())
    #Prepare for merging
    LibASeqMerge <- LibASeqDesign[,c("SeqID", "TF.Trp53.affinity","TF.Trp53.nrepeats","TF.Trp53.orientation","Universal.Spacing")]
    colnames(LibASeqMerge)[1] <- "ControlFor"
    #Generate Trp53 data.frame
    Trp53meta <- inner_join(Design_meta_controls, LibASeqMerge, by = "ControlFor" )
    colnames(Trp53meta)[2] <- "CRS"
    Trp53meta <- Trp53meta %>%
      mutate(TF.Trp53.affinity_nrepeats = TF.Trp53.affinity* TF.Trp53.nrepeats)
    CRS.all.Trp53 <- inner_join(Trp53meta, CRS.all, by = "CRS")
    CRS.all.filtered.Trp53 <- inner_join(Trp53meta, CRS.all.filtered, by = "CRS")
    
    #Calculate cell state cor
    report$Trp53_norm_cor <- cor(x= CRS.all.Trp53$norm.1, y = CRS.all.Trp53$norm.2, method = "pearson")
    report$Trp53_cor_activity <- cor(x= CRS.all.Trp53$TF.Trp53.affinity_nrepeats, y = CRS.all.Trp53$mean.norm, method = "pearson")  
    report$Trp53_norm_cor_filter <- cor(x= CRS.all.filtered.Trp53$norm.1, y = CRS.all.filtered.Trp53$norm.2, method = "pearson")
    report$Trp53_cor_activity_filter <- cor(x= CRS.all.filtered.Trp53$TF.Trp53.affinity_nrepeats, y = CRS.all.filtered.Trp53$mean.norm, method = "pearson") 
    
    if(plot) {
      png(file.path(outdir, run_id, "Trp53_full_activity.png"),width=800,height=600,res=100)
      print(ggplot(data = CRS.all.Trp53) + geom_point(aes(x = TF.Trp53.affinity_nrepeats, y = mean.norm, color = TF.Trp53.affinity)) + theme_bw() + xlab("Affinity * nrepeats") + ylab("Mean norm GRE activity"))
      dev.off()
      
      png(file.path(outdir, run_id, "Trp53_filtered_activity.png"),width=800,height=600,res=100)
      print(ggplot(data = CRS.all.filtered.Trp53) + geom_point(aes(x = TF.Trp53.affinity_nrepeats, y = mean.norm, color = TF.Trp53.affinity)) + theme_bw() + xlab("Affinity * nrepeats") + ylab("Mean norm GRE activity"))
      dev.off()
    }
  }

  if (plot) {
    png(file.path(outdir, run_id, "CRS_full_replicates_%02d.png"),width=800,height=600,res=100)
    heatscatter(CRS.all$norm.1,CRS.all$norm.2, main = "normalized", cor=T, method = "pearson")
    heatscatter(log2(CRS.all$RNA.1+1),log2( CRS.all$RNA.2+1), main = "normalized", cor=T, method = "pearson")
    heatscatter(log2(CRS.all$DNA.1+1),log2( CRS.all$DNA.2+1), main = "normalized", cor=T, method = "pearson")
    dev.off()
    
    #Filtered CRSs
    if(length(CRS.all.filtered$norm.1) > 0 | length(CRS.all.filtered$norm.2) > 0) {
      png(file.path(outdir, run_id, "CRS_filtered_replicates_%02d.png"),width=800,height=600,res=100)
      heatscatter(CRS.all.filtered$norm.1,CRS.all.filtered$norm.2, main = "normalized", cor=T, method = "pearson")
      heatscatter(log2(CRS.all.filtered$RNA.1+1),log2(CRS.all.filtered$RNA.2+1), main = "normalized", cor=T, method = "pearson")
      heatscatter(log2(CRS.all.filtered$DNA.1+1),log2(CRS.all.filtered$DNA.2+1), main = "normalized", cor=T, method = "pearson")
      dev.off()
    }

}
  
  y <- subset(CRS.all, !is.infinite(norm.1) & !is.infinite(norm.2))
  
  report$RNA_cor <- cor(log2(y$RNA.1+1),log2(y$RNA.2+1), use="pairwise.complete.obs")
  report$DNA_cor <- cor(log2(y$DNA.1+1),log2(y$DNA.2+1), use="pairwise.complete.obs")
  ifelse(length(y$norm.1) > 0 & length(y$norm.2) > 0, report$norm_cor <- cor(y$norm.1,y$norm.2, use="pairwise.complete.obs"), report$norm_cor <- NA)
  report$norm_n <- sum(!is.na(y$norm.1) & !is.na(y$norm.2))
  y <- subset(y, pass)
  ifelse(length(y$norm.1) > 0 & length(y$norm.2) > 0, report$norm_cor_filter <- cor(y$norm.1,y$norm.2, use="pairwise.complete.obs"), report$norm_cor_filter <- NA)
  report$norm_n_filter <- sum(!is.na(y$norm.1) & !is.na(y$norm.2))

  if(save_report) { 
    report_safe <- as.data.frame(t(as.data.frame(report)))
    colnames(report_safe) <- paste(run_id, c("__Replicate1","__Replicate2"), sep = "")
    save(report_safe, CRS.all, file = file.path(outdir, run_id, "CRS_summary_report.rda"))
    }
  if (return_report) return(list(CRS.all, report)) else return(CRS.all)
}
