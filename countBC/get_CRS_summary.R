#systematically create CRS level summaries from the perl pipeline output

#example 1 - process several folders in parallel:
# folders <- c(EoBaso = "~/cluster_mount/project/SCG4SYN/Experiments/220908_RF_HSC_MPRA_screen/Analysis/count_BC/NextSeq2000/220918_EoBaso/",
#              EoBaso500 ="~/cluster_mount/project/SCG4SYN/Experiments/220908_RF_HSC_MPRA_screen/Analysis/count_BC/NextSeq500/220918_EoBaso/",
#              EoBasoMerge = "~/cluster_mount/project/SCG4SYN/Experiments/220908_RF_HSC_MPRA_screen/Analysis/count_BC/NextSeq_merge/220918_EoBaso",
#              MonoPre = "~/cluster_mount/project/SCG4SYN/Experiments/220908_RF_HSC_MPRA_screen/Analysis/count_BC/NextSeq2000/220918_MonoPre/")
              
#CRS.all_UMI5_DNAisUMI_crsfilter20 <- parallel::mcmapply(get_CRS_summary,folders, paste(names(folders), "_UMI5_DNAisUMI_crsfilter20", sep = ""), MoreArgs = list(use.col = "UMI5", return_report = T, save_report=T, dna.crs.is.umi = T, crs_filter = 20), SIMPLIFY = F, mc.cores = 4)

#example 2 - best conditions so far


#final <-  mcmapply(get_CRS_summary,folders,names(folders), replicate(7, NULL), c("UMI3","UMI3","UMI3","UMI10","UMI10","UMI10","UMI10"),
#                    MoreArgs = list(outdir = "allout/finalrerun/", return_report = T, 
#                                    dna.crs.is.umi = T, crs_filter = 10, filter_controlGRE=T,
#                                    bc_crs_ass = "/users/lvelten/project/SCG4SYN/Experiments/220908_RF_HSC_MPRA_screen/Files/220918_CRS_BC_ass/mapped.filtered.custom.csv"), 
#                    SIMPLIFY = F, mc.cores = 7, mc.preschedule = F)

##
'%nin%' = Negate('%in%')

#Function
get_CRS_summary <- function(f, run_id, Experiments = NULL,
         use.col = "UMI5", #the column in the PERL pipeline output to use (e.g. UMI5: UMIs with >=5 reads)
         plot =F, #create plots on the way in the directory specified as outdir
         outdir = "~/cluster_mount/project/SCG4SYN/Experiments/220908_RF_HSC_MPRA_screen/Analysis/count_BC/plots", 
         min.dna.umi = 1, #to filter at the level of CRS barcodes: Only include CRS BC covered with at least this many UMIs at DNA level
         min.rna.umi = 0, #to filter at the level of CRS barcodes: Only include CRS BC covered with at least this many UMIs at RNA level
         umi.filter.connection = `&`, #DNA filter AND RNA filter or DNA filter OR RNA filter?
         dna.crs.is.umi = T, #whether to treat DNA CRS observations as single molecules no matter how many UMIs there are
         rna.crs.is.umi = F,
         barcodeWhitelist = NULL, #a filtered bc crs association obejct (created from the file in bc_crs_ass if not supplied)
         homopolymer = 0, #max stretch of homopolymer allowed in barcode. If 0, no filter is applied
         bc_crs_ass = "~/cluster_mount/project/SCG4SYN/Experiments/220908_RF_HSC_MPRA_screen/Files/220918_CRS_BC_ass/mapped.filtered.custom.csv",
         crs_filter =10, #min UMIs per CRS to be flagged as pass
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
  report <- list()
  
  if (!dir.exists(file.path(outdir, run_id))) {dir.create(file.path(outdir, run_id), showWarnings = T, recursive = T )}

  #read in data (skipped if data is provided)
  if (is.null(Experiments)){
    rna <- list.files(f, "RNA\\d.csv")
    dna <- list.files(f, "DNA\\d.csv")
    RNA <- lapply(file.path(f, rna), read.csv)
    DNA <- lapply(file.path(f, dna), read.csv)
    Experiments <- mapply(merge, DNA, RNA, MoreArgs = list( by = c("BARCODE","CRS"), suffixes = c(".dna",".rna"), all = T), SIMPLIFY =FALSE )
    
  }
  
  if (return_experiments) return(Experiments)
  
  if (homopolymer > 0) {
    filter <- sapply(c("A","C","G","T"), function(x) paste(rep(x, homopolymer), collapse = "") )
    filter <- paste(filter, collapse = "|")
    Experiments <- lapply(Experiments, function(exx) subset(exx, !grepl(filter, BARCODE)))
  }
  
  if (is.null(barcodeWhitelist)) {
    BC_CRS <- read.csv(bc_crs_ass)
    barcodeWhitelist <- subset(BC_CRS, MEANMATCHES >= 290 & DEVIANTREADS / (READS + DEVIANTREADS) < 0.2)$BARCODE
  }
  Experiments <- lapply(Experiments, function(exx) subset(exx, BARCODE %in% barcodeWhitelist))
  
  #replace NAs by 0
  Experiments <- lapply(Experiments, function(exx) {
    for (i in 3:ncol(exx)) exx[is.na(exx[,i]),i] <- 0
    exx$norm <- exx[,paste0(use.col,".rna")] / exx[,paste0(use.col,".dna")]
    exx$USE.dna <- exx[,paste0(use.col,".dna")]
    exx$USE.rna <- exx[,paste0(use.col,".rna")]
    exx
  })
  
  #Filter control GREs
  if(filter_controlGRE) {
    Experiments <- lapply(Experiments, function(exx){
      exx <- exx[which(!grepl(pattern = "^Ctrl_*", x = exx$CRS, ignore.case = T)),]
    })
  }
  
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
 
  
  Experiments.filtered <- lapply(Experiments, function(x) 
    subset(x, umi.filter.connection(UMI.dna >= min.dna.umi, UMI.rna >= min.rna.umi)))
  
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
  #Add meta info
  CRS.all$TF <- gsub(".*tf_([^_]+)_.+","\\1",CRS.all$CRS)
  CRS.all$nrepeats <- gsub(".*nr_([^_]+)_.+","\\1",CRS.all$CRS)
  CRS.all$affinity <- gsub(".*aff_([^_]+)_.+","\\1",CRS.all$CRS)
  CRS.all$affinity <- gsub("NULL", "c-0-0-" , CRS.all$affinity)
  CRS.all$orientation <- gsub(".*ori_([^_]+)_.+","\\1",CRS.all$CRS)
  CRS.all$orientation <- ifelse(CRS.all$orientation %in% "fwd", "fwd", ifelse(CRS.all$orientation %in% "rev", "rev", "tandem"))
  
  CRS.all <- CRS.all %>% 
    dplyr::filter(norm.1 %nin% Inf & norm.1 %nin% -Inf & norm.1 %nin% NA & norm.1 %nin% NaN & norm.2 %nin% Inf & norm.2 %nin% -Inf & norm.2 %nin% NA & norm.2 %nin% NaN)
  
  CRS.all <- CRS.all %>%
    rowwise() %>%
    dplyr::mutate(mean.norm = mean(c(norm.1, norm.2))) %>%
    ungroup()
  
  #melt data.frame for visualization
  plf <- melt(CRS.all[,c("CRS","TF","nrepeats","affinity","orientation","norm.1","norm.2","mean.norm")], id.vars = c("CRS","TF", "nrepeats","affinity","orientation"))
  plf$nrepeats <- as.numeric(plf$nrepeats)
  
  #Create filtered version
  CRS.all.filtered <- subset(CRS.all, pass )
  #melt data.frame for visualization
  plf.filtered <- melt(CRS.all.filtered[,c("CRS","TF","nrepeats","affinity","orientation","norm.1","norm.2","mean.norm")], id.vars = c("CRS","TF", "nrepeats","affinity","orientation"))
  plf.filtered$nrepeats <- as.numeric(plf.filtered$nrepeats)
  
  if (plot) {
    png(file.path(outdir, run_id, "CRS_full_replicates_%02d.png"),width=800,height=600,res=100)
    heatscatter(CRS.all$norm.1,CRS.all$norm.2, main = "normalized", cor=T, method = "pearson")
    heatscatter(log2(CRS.all$RNA.1+1),log2( CRS.all$RNA.2+1), main = "normalized", cor=T, method = "pearson")
    heatscatter(log2(CRS.all$DNA.1+1),log2( CRS.all$DNA.2+1), main = "normalized", cor=T, method = "pearson")
    dev.off()
    
    pdf(file.path(outdir, run_id, "CRS_full_replicates_GREvisuals_%02d.pdf"), onefile = F)
    print(ggplot(aes(x = reorder(TF, value), y = value, color = variable), data=plf %>% filter(variable%in% "mean.norm")) + 
            geom_boxplot(aes(color = variable)) + 
            scale_color_discrete(name = "Replicate ", labels = c("Average")) + 
            geom_hline(aes(yintercept = 0), linetype = "dashed") + 
            ylab("Normalized RNA/DNA ratio") + xlab("Transcription factor") + 
            theme_minimal() +
            theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)))
    print(ggplot(aes(x = reorder(TF, value), y = value, color = variable), data=plf %>% filter(variable %nin% "mean.norm")) + 
            geom_boxplot(aes(color = variable)) + 
            scale_color_discrete(name = "Replicate ", labels = c("1", "2")) + 
            geom_hline(aes(yintercept = 0), linetype = "dashed") + 
            ylab("Normalized RNA/DNA ratio") + xlab("Transcription factor") + 
            theme_minimal() +
            theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)))
    print(ggplot(aes(x = reorder(TF, value), y = value, color = variable), data=plf ) + 
      geom_boxplot(aes(color = variable)) + 
      scale_color_discrete(name = "Replicate ", labels = c("1", "2", "Average")) + 
      geom_hline(aes(yintercept = 0), linetype = "dashed") + 
      ylab("Normalized RNA/DNA ratio") + xlab("Transcription factor") + 
      theme_minimal() +
      theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)))
    print(ggplot(aes(x = nrepeats, y = value, color = variable), data=plf %>% filter(variable%in% "mean.norm")) + 
      facet_wrap(~ TF) + 
      geom_smooth(aes(color = variable)) + 
      scale_color_discrete(name = "Replicate ", labels = c("Average")) + 
      geom_hline(yintercept = 0, linetype = "dashed") + 
      ylab("Normalized RNA/DNA ratio") + xlab("Number of repeats") + 
      theme_minimal()) 
    print(ggplot(aes(x = nrepeats, y = value, color = variable), data=plf %>% filter(variable%nin% "mean.norm")) + 
      facet_wrap(~ TF) + 
      geom_smooth(aes(color = variable)) + 
      scale_color_discrete(name = "Replicate ", labels = c("1","2")) + 
      geom_hline(yintercept = 0, linetype = "dashed") + 
      ylab("Normalized RNA/DNA ratio") + xlab("Number of repeats") + 
      theme_minimal()) 
    print(ggplot(aes(x = affinity, y = value), data=plf %>% filter(variable %in% "mean.norm")) + 
      facet_wrap(~ TF) + geom_boxplot(aes(color = variable)) + 
      scale_color_discrete(name = "Replicate ", labels = c("Average")) + 
      geom_hline(yintercept = 0, linetype = "dashed") + 
      ylab("Normalized RNA/DNA ratio") + xlab("Binding site affinity") + 
      theme_minimal() +
      theme(axis.text.x = element_text(size = 6, angle = 90)))
    print(ggplot(aes(x = affinity, y = value), data=plf %>% filter(variable %nin% "mean.norm")) + 
      facet_wrap(~ TF) + geom_boxplot(aes(color = variable)) + 
      scale_color_discrete(name = "Replicate ", labels = c("1", "2")) + 
      geom_hline(yintercept = 0, linetype = "dashed") + 
      ylab("Normalized RNA/DNA ratio") + xlab("Binding site affinity") + 
      theme_minimal() +
      theme(axis.text.x = element_text(size = 6, angle = 90)))
    print(ggplot(aes(x = orientation, y = value), data=plf) + 
      facet_wrap(~ TF) + 
      geom_boxplot(aes(color = variable)) + 
      scale_color_discrete(name = "Replicate ", labels = c("1", "2","Average")) + 
      geom_hline(yintercept = 0, linetype = "dashed") +
      theme_minimal() + 
      theme(axis.text.x = element_text(angle=90)) + 
      ylab("Normalized RNA/DNA ratio") + xlab("Orientation of the motif")) 
    dev.off()
    
    #Filtered CRSs
    png(file.path(outdir, run_id, "CRS_filtered_replicates_%02d.png"),width=800,height=600,res=100)
    heatscatter(CRS.all.filtered$norm.1,CRS.all.filtered$norm.2, main = "normalized", cor=T, method = "pearson")
    heatscatter(log2(CRS.all.filtered$RNA.1+1),log2(CRS.all.filtered$RNA.2+1), main = "normalized", cor=T, method = "pearson")
    heatscatter(log2(CRS.all.filtered$DNA.1+1),log2(CRS.all.filtered$DNA.2+1), main = "normalized", cor=T, method = "pearson")
    dev.off()
    

    pdf(file.path(outdir, run_id, "CRS_filtered_replicates_GREvisuals_%02d.pdf"), onefile = F)
    print(ggplot(aes(x = reorder(TF, value), y = value, color = variable), data=plf.filtered %>% filter(variable %in% "mean.norm") ) + 
            geom_boxplot(aes(color = variable)) + 
            scale_color_discrete(name = "Replicate ", labels = c("Average")) + 
            geom_hline(aes(yintercept = 0), linetype = "dashed") + 
            ylab("Normalized RNA/DNA ratio") + xlab("Transcription factor") + 
            theme_minimal() +
            theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)))
    print(ggplot(aes(x = reorder(TF, value), y = value, color = variable), data=plf.filtered %>% filter(variable %nin% "mean.norm")) + 
            geom_boxplot(aes(color = variable)) + 
            scale_color_discrete(name = "Replicate ", labels = c("1", "2")) + 
            geom_hline(aes(yintercept = 0), linetype = "dashed") + 
            ylab("Normalized RNA/DNA ratio") + xlab("Transcription factor") + 
            theme_minimal() +
            theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)))
    print(ggplot(aes(x = reorder(TF, value), y = value, color = variable), data=plf.filtered ) + 
            geom_boxplot(aes(color = variable)) + 
            scale_color_discrete(name = "Replicate ", labels = c("1", "2", "Average")) + 
            geom_hline(aes(yintercept = 0), linetype = "dashed") + 
            ylab("Normalized RNA/DNA ratio") + xlab("Transcription factor") + 
            theme_minimal() +
            theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)))
    print(ggplot(aes(x = nrepeats, y = value, color = variable), data=plf.filtered %>% filter(variable%in% "mean.norm")) + 
            facet_wrap(~ TF) + 
            geom_smooth(aes(color = variable)) + 
            scale_color_discrete(name = "Replicate ", labels = c("Average")) + 
            geom_hline(yintercept = 0, linetype = "dashed") + 
            ylab("Normalized RNA/DNA ratio") + xlab("Number of repeats") + 
            theme_minimal()) 
    print(ggplot(aes(x = nrepeats, y = value, color = variable), data=plf.filtered %>% filter(variable%nin% "mean.norm")) + 
      facet_wrap(~ TF) + 
      geom_smooth(aes(color = variable)) + 
      scale_color_discrete(name = "Replicate ", labels = c("1","2")) + 
      geom_hline(yintercept = 0, linetype = "dashed") + 
      ylab("Normalized RNA/DNA ratio") + xlab("Number of repeats") + 
      theme_minimal()) 
    print(ggplot(aes(x = affinity, y = value), data=plf.filtered %>% filter(variable %in% "mean.norm")) + 
      facet_wrap(~ TF) + geom_boxplot(aes(color = variable)) + 
      scale_color_discrete(name = "Replicate ", labels = c("Average")) + 
      geom_hline(yintercept = 0, linetype = "dashed") + 
      ylab("Normalized RNA/DNA ratio") + xlab("Binding site affinity") + 
      theme_minimal() +
      theme(axis.text.x = element_text(size = 6, angle = 90)))
    print(ggplot(aes(x = affinity, y = value), data=plf.filtered %>% filter(variable %nin% "mean.norm")) + 
      facet_wrap(~ TF) + geom_boxplot(aes(color = variable)) + 
      scale_color_discrete(name = "Replicate ", labels = c("1", "2")) + 
      geom_hline(yintercept = 0, linetype = "dashed") + 
      ylab("Normalized RNA/DNA ratio") + xlab("Binding site affinity") + 
      theme_minimal() +
      theme(axis.text.x = element_text(size = 6, angle = 90)))
    print(ggplot(aes(x = orientation, y = value), data=plf.filtered) + 
      facet_wrap(~ TF) + 
      geom_boxplot(aes(color = variable)) + 
      scale_color_discrete(name = "Replicate ", labels = c("1", "2","Average")) + 
      geom_hline(yintercept = 0, linetype = "dashed") +
      theme_minimal() + 
      theme(axis.text.x = element_text(angle=90)) + 
      ylab("Normalized RNA/DNA ratio") + xlab("Orientation of the motif")) 
    dev.off()
}
  
  y <- subset(CRS.all, !is.infinite(norm.1) & !is.infinite(norm.2))
  report$norm_cor <- cor(y$norm.1,y$norm.2, use="pairwise.complete.obs")
  report$norm_n <- sum(!is.na(y$norm.1) & !is.na(y$norm.2))
  y <- subset(y, pass)
  report$norm_cor_filter <- cor(y$norm.1,y$norm.2, use="pairwise.complete.obs")
  report$norm_n_filter <- sum(!is.na(y$norm.1) & !is.na(y$norm.2))
  if(save_report) { 
    report_safe <- as.data.frame(t(as.data.frame(report)))
    colnames(report_safe) <- paste(run_id, c("__Replicate1","__Replicate2"), sep = "")
    save(report_safe, file = file.path(outdir, run_id, "CRS_summary_report.rda"))
    }
  if (return_report) return(list(CRS.all, report)) else return(CRS.all)
}
