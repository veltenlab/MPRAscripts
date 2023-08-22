#compare seq depth per umi
require(plyr)
require(ggplot2)
require(reshape2)
require(readr)
require(dplyr)
require(stringr)
require(tibble)
require(tidyr)
umistats.H.files <- system("find ~/cluster_mount/project/SCG4SYN/Pipelines/MPRA/2_BC_count/output/4_HSC_LibraryH_minP/filtered/ -iname '*umistat.txt'", intern = T)

#umistats.A.files <- system("find ~/cluster_mount/lvelten/Analysis/SCG4SYN/LibA_HSC/raw/allout -iname '*umistat.txt'", intern = T)
umistats.A.files <- system("find ~/cluster_mount/project/SCG4SYN/Pipelines/MPRA/2_BC_count/output/2_HSC_LibraryA_minPInitial/filtered/ -iname '*umistat.txt'", intern = T)

#umistats.Arepfiles <- system("find ~/cluster_mount/lvelten/Analysis/SCG4SYN/LibA_HSC_2/raw/allout -iname '*umistat.txt'", intern = T)
umistats.Arepfiles <- system("find ~/cluster_mount/project/SCG4SYN/Pipelines/MPRA/2_BC_count/output/3_HSC_LibraryA_minPOptimized/filtered/ -iname '*umistat.txt'", intern = T)

#Read in read statistic
read_stats <- read.csv("~/cluster_mount/project/SCG4SYN/Pipelines/MPRA/2_BC_count/data/Sequencing_statistics_BS3.csv", sep =";")
read_stats$exp <- gsub(pattern = "_lib", replacement = "", x = read_stats$exp)

#Read in Umi files
umistats.A <- lapply(umistats.A.files, function(f) {
  umistat <- read.csv(f, header=F)
  umistat$library <- "A"
  umistat$exp <- gsub(".+/([^/]+).umistat.txt","\\1",f)
  umistat$DNA <- grepl("DNA|K562D",umistat$exp)
  umistat
})
umistats.A <- do.call(rbind, umistats.A)

umistats.H <- lapply(umistats.H.files, function(f) {
  umistat <- read.csv(f, header=F)
  umistat$library <- "H"
  umistat$exp <- gsub(".+/([^/]+).umistat.txt","\\1",f)
  umistat$DNA <- grepl("DNA",umistat$exp)
  umistat
})
umistats.H <- do.call(rbind, umistats.H)

umistats.Arep <- lapply(umistats.Arepfiles, function(f) {
  umistat <- read.csv(f, header=F)
  umistat$library <- "Arep"
  umistat$exp <- gsub(".+/([^/]+).umistat.txt","\\1",f)
  umistat$DNA <- grepl("DNA",umistat$exp)
  umistat
})
umistats.Arep <- do.call(rbind, umistats.Arep)

#Merging UMIs
umistats <- rbind(umistats.A, umistats.H, umistats.Arep)
umistats <- umistats %>% mutate(Cellstate = sub("(?=[DR]NA).*", "", x = exp, perl = T))

#whole overview
ggplot(data = umistats) + geom_line(aes(x = log10(V1), y = log10(V2), color = library, group = exp))  + xlab("log10(Reads per barcode)") + ylab("log10(Barcodes)")
ggplot(data = umistats) + geom_line(aes(x = log10(V1), y = log10(V2), color = library, group = exp)) + facet_wrap(~Cellstate, nrow=2)  + xlab("log10(Reads per barcode)") + ylab("log10(Barcodes)")
ggplot(data = umistats) + geom_line(aes(x = log10(V1), y = log10(V2), color = library, group = exp)) + facet_wrap(~library+DNA, nrow=3)+ xlab("log10(Reads per barcode)") + ylab("log10(Barcodes)")

#LibraryH
ggplot(data = umistats %>% filter(library %in% "H")) + geom_line(aes(x = log10(V1), y = log10(V2), color = exp, group = exp)) +ggtitle("LibraryH") + xlab("log10(Reads per barcode)") + ylab("log10(Barcodes)")

ggplot(data = umistats %>% filter(library %in% "H") %>% inner_join(., read_stats %>% select(exp, TotalReads), by = "exp")) + geom_line(aes(x = log10(V1), y = log10(V2), color = TotalReads, group = exp)) + facet_wrap(~DNA, ncol=2) + ggtitle("LibraryH DNA") + xlab("log10(Reads per barcode)") + ylab("log10(Barcodes)")
ggplot(data = umistats %>% filter(library %in% "H") %>% inner_join(., read_stats %>% select(exp, TotalReads), by = "exp")) + geom_line(aes(x = log10(V1), y = log10(V2), color = TotalReads, group = exp)) + facet_wrap(~Cellstate + DNA, ncol=7) + ggtitle("LibraryH DNA") + xlab("log10(Reads per barcode)") + ylab("log10(Barcodes)")
ggplot(data = umistats %>% filter(library %in% "H") %>% inner_join(., read_stats %>% select(exp, TotalReads), by = "exp")) + geom_line(aes(x = log10(V1), y = log10(V2), color = TotalReads, group = exp)) + facet_wrap(~Cellstate + DNA, ncol=2) + ggtitle("LibraryH DNA") + xlab("log10(Reads per barcode)") + ylab("log10(Barcodes)") 

#LibraryH with DNA or RNA
ggplot(data = umistats %>% filter(library %in% "H" & DNA %in% T)) + geom_line(aes(x = log10(V1), y = log10(V2), color = exp, group = exp)) + ggtitle("LibraryH DNA") + xlab("log10(Reads per barcode)") + ylab("log10(Barcodes)") 
ggplot(data = umistats %>% filter(library %in% "H" & DNA %in% T) %>% inner_join(., read_stats %>% select(exp, TotalReads), by = "exp")) + geom_line(aes(x = log10(V1), y = log10(V2), color = TotalReads, group = exp)) + ggtitle("LibraryH DNA") + xlab("log10(Reads per barcode)") + ylab("log10(Barcodes)")
ggplot(data = umistats %>% filter(library %in% "H" & DNA %in% F) %>% inner_join(., read_stats %>% select(exp, TotalReads), by = "exp")) + geom_line(aes(x = log10(V1), y = log10(V2), color = TotalReads, group = exp)) + ggtitle("LibraryH RNA") + xlab("log10(Reads per barcode)") + ylab("log10(Barcodes)")

ggplot(data = umistats %>% filter(library %in% "H" & DNA %in% T) ) + geom_line(aes(x = log10(V1), y = log10(V2), color = exp, group = exp))+ facet_wrap(~Cellstate, ncol=4) + ggtitle("LibraryH DNA") + xlab("log10(Reads per barcode)") + ylab("log10(Barcodes)")
ggplot(data = umistats %>% filter(library %in% "H" & DNA %in% T) %>% inner_join(., read_stats %>% select(exp, TotalReads), by = "exp") ) + geom_line(aes(x = log10(V1), y = log10(V2), color = TotalReads, group = exp))+ facet_wrap(~Cellstate, ncol=4) + ggtitle("LibraryH DNA") + xlab("log10(Reads per barcode)") + ylab("log10(Barcodes)")
ggplot(data = umistats %>% filter(library %in% "H" & DNA %in% F) %>% inner_join(., read_stats %>% select(exp, TotalReads), by = "exp") ) + geom_line(aes(x = log10(V1), y = log10(V2), color = TotalReads, group = exp))+ facet_wrap(~Cellstate, ncol=4) + ggtitle("LibraryH RNA") + xlab("log10(Reads per barcode)") + ylab("log10(Barcodes)")


#number of reads
nreads <- ddply(umistats, c("library","exp","DNA"), summarise, reads = sum(V1*V2), barcodes = sum(V2[V1 >= ifelse(library == "A", 10,1)]))
qplot(x = library, y = reads, data = nreads, geom="boxplot")

ggplot(data = nreads %>% filter(library %in% "H")%>% inner_join(., read_stats %>% select(exp, TotalReads), by = "exp")) + geom_point(aes(x= exp, y= reads/TotalReads)) + facet_wrap(~DNA, ncol=4, scales = "free_x") + theme(axis.text.x = element_text(angle= 90)) + xlab("Library H sample") + ylab("Mapped versus total reads [%]")


#how manz barcodes do we think there are per population?
#how many were there in Lib A?
#estimate from experiment vs. what the data tells us
#Probably there are 10-20 fold more barcodes, then we need to sequence 10-20 fold deeper to get
qplot(x = reads, y = barcodes, data = nreads, color = library, shape = DNA, log="xy") +xlab("Mapped reads")
qplot(x = reads, y = barcodes, data = nreads %>% filter(library %in% "H"), color = exp, shape = DNA, log="xy")   +xlab("Mapped reads")

qplot(x = Freq, y = barcodes, data = nreads %>% filter(library %in% "H")%>% inner_join(., read_stats %>% select(exp, TotalReads), by = "exp") %>% mutate(Freq=reads/TotalReads), color = exp, shape = DNA, log="y") + facet_wrap(~DNA, ncol=4, scales = "fixed") +xlab("Mapped versus total reads [%]")
#dpes this plot show that on dna level (counts not umis, but integrations, right) there are very roughly 1M integrations in the libH
#MoI was 10, how many cells were sorted?
# Also look at reads matching vs. not matching the dictionary.
ddply(umistats, c("library"), summarise, reads = sum(V1*V2))

#
read_stats_mod <- read_stats
read_stats_mod$DNA <- grepl("DNA",read_stats_mod$exp)
read_stats_mod$sample <- gsub("(\\d).+[DR]NA", "C\\1__Replicate\\2", read_stats_mod$exp)
read_stats_mod.short <- dcast(read_stats_mod,sample ~DNA, value.var = "TotalReads")

#how many DNA barcodes do we have covered in each experiment on RNA and DNA?

#umistats.H.files <- system("find /home/lars/cluster/lvelten/Analysis/SCG4SYN/LibH/processed/allout/UMI1 -iname 'CRS_summary_report.rda'", intern = T)
report_libH <- system("find ~/cluster_mount/project/SCG4SYN/Pipelines/MPRA/3_Report_screen/output/4_HSC_LibraryH_minP/filtered/Umi1_DnaIsUmiT_crsfilter20_minDnaUmi1_minRnaUmi0_Reads1_InputF_PseudoF/ -iname 'CRS_summary_report.rda'", intern = T)
"~/cluster_mount/project/SCG4SYN/Pipelines/MPRA/3_Report_screen/output/4_HSC_LibraryH_minP/filtered/Umi1_DnaIsUmiT_crsfilter20_minDnaUmi1_minRnaUmi0_Reads1_InputF_PseudoT_pseudoDna_pseudoValue1//"

#report_libH <- system("find ~/cluster_mount/project/SCG4SYN/Pipelines/MPRA/3_Report_screen/output/4_HSC_LibraryH_minP/filtered/Umi1_DnaIsUmiT_crsfilter20_minDnaUmi1_minRnaUmi0_Reads1_InputF_PseudoT_pseudoDna_pseudoValue1/ -iname 'CRS_summary_report.rda'", intern = T)


#for libH:
processing_report <- lapply(report_libH, function(x) {
  load(x)
 return(report_safe)
})
#
#names(report_libH) <- str_match(report_libH, "filtered/(.*?)/State_")[, 2]

processing_report <- as.data.frame(t(do.call(cbind, processing_report)))
processing_report$sample <- rownames(processing_report)
processing_report$sample <- gsub(pattern = "State_(\\d)(.+)__", replacement = "C\\1__", x = processing_report$sample)
processing_report.selectRNA <- processing_report
processing_report.selectRNA <- merge(processing_report.selectRNA, read_stats_mod.short)

qplot(y = crs_bc_covered.RNA, x = `FALSE`,data = processing_report.selectRNA, log="xy") + ylab("RNA Barcodes") + xlab("Reads RNA libraries")
qplot(y = crs_bc_covered.DNA, x = `TRUE`,data = processing_report.selectRNA, log="xy") + ylab("DNA Barcodes") + xlab("Reads")

qplot(x = sample, y = RNA_cor, size = crs_bc_covered.DNA,data=processing_report.selectRNA) + coord_flip()
qplot(x = sample, y = DNA_cor, size = crs_bc_covered.DNA,data=processing_report.selectRNA) + coord_flip()
qplot(x = sample, y = DNA_cor, size = crs_bc_covered.RNA,data=processing_report.selectRNA) + coord_flip()
qplot(x = sample, y = norm_cor, size = crs_bc_covered.RNA,data=processing_report.selectRNA) + coord_flip()
qplot(x = sample, y = norm_cor, size = crs_bc_covered.DNA,data=processing_report.selectRNA) + coord_flip()

processing_report_melt <- melt(subset(processing_report, select = c("crs_bc_covered.RNA" ,"crs_bc_covered.DNA", "crs_bc_filter.RNA", "crs_bc_filter.DNA", "sample")), id.vars = "sample")
ggplot(aes(x = sample, y = value, fill =variable), data = subset(processing_report_melt, variable == "crs_bc_covered.RNA"))+ geom_col(position = position_dodge()) +coord_flip()
ggplot(aes(x = sample, y = value, fill =variable), data = subset(processing_report_melt, variable %in% c("crs_bc_covered.DNA","crs_bc_filter.RNA")))+ geom_col(position = position_dodge()) +coord_flip()

# Load in all reports and compare them
reports_libH <- system("find ~/cluster_mount/project/SCG4SYN/Pipelines/MPRA/3_Report_screen/output/4_HSC_LibraryH_minP/filtered/ -iname 'CRS_summary_report.rda'", intern = T)

reports_libH <- reports_libH[!grepl(pattern = "_PseudoT_", x = reports_libH)]
#reports_libH <- reports_libH[grepl(pattern = "_PseudoT_", x = reports_libH)]

names(reports_libH) <- str_match(reports_libH, "filtered/(.*?)/State_")[, 2]

processing_reports <- lapply(reports_libH, function(x) {
  load(x)
  return(report_safe)
})

processing_reports2 <- as.data.frame(t(do.call(cbind, processing_reports)))
processing_reports2 <- processing_reports2 %>% rownames_to_column(var ="ID")
processing_reports2$sample <- gsub(pattern = ".*State_(\\d)(.+)__", replacement = "C\\1__", x = processing_reports2$ID)
processing_reports2$Run <- gsub(pattern = ".State.*", replacement = "", x = processing_reports2$ID)
processing_reports2 <- processing_reports2 %>% separate(data = ., col = "Run", into = c("use.col", "dna.crs.is.umi","crs_filter","min.dna.umi", "min.rna.umi", "AssReadsFilter","normalize.to.input", "pseudocount"), sep = "_", remove = F) 
#processing_reports2 <- processing_reports2 %>% separate(data = ., col = "Run", into = c("use.col", "dna.crs.is.umi","crs_filter","min.dna.umi", "min.rna.umi", "AssReadsFilter","normalize.to.input", "pseudocount", "pseudotype", "pseudovalue"), sep = "_", remove = F) 

processing_reports2 <- processing_reports2 %>% column_to_rownames(var = "ID")
processing_reports2.selectRNA <- processing_reports2
processing_reports2.selectRNA <- merge(processing_reports2.selectRNA, read_stats_mod.short)

#Set factors
processing_reports2$Run <- factor(processing_reports2$Run, levels = processing_reports2 %>%
                                    group_by(Run) %>%
                                    summarise(score=mean(norm_cor)) %>% arrange(desc(score)) %>% select(Run) %>% pull())

#ggplot(data = processing_reports2) +geom_point(aes(x=sample, y=RNA_cor, size = crs_bc_covered.DNA)) + facet_wrap(~Run) + theme(axis.text.x = element_text(angle = 90), strip.text = element_text(size = 6)) + ggtitle("RNA correlation") + xlab("Samples") + ylab("RNA cor")

ggplot(data = processing_reports2) +geom_point(aes(x=reorder(sample, -RNA_cor), y=RNA_cor, size = crs_bc_covered.DNA, color = Run, group = Run), alpha = 0.2) + geom_line(aes(x=reorder(sample, -RNA_cor), y=RNA_cor, color = Run, group = Run)) + theme_bw() + theme(axis.text.x = element_text(angle = 90)) + ggtitle("RNA correlation") + xlab("Samples") + ylab("RNA cor")
ggplot(data = processing_reports2) +geom_point(aes(x=reorder(sample, -DNA_cor), y=DNA_cor, size = crs_bc_covered.DNA, color = Run, group = Run), alpha = 0.2) + geom_line(aes(x=reorder(sample, -DNA_cor), y=DNA_cor, color = Run, group = Run)) + theme_bw() + theme(axis.text.x = element_text(angle = 90)) + ggtitle("DNA correlation") + xlab("Samples") + ylab("DNA cor")
ggplot(data = processing_reports2 %>% filter(use.col != "Umi3")) +geom_point(aes(x=reorder(sample, -norm_cor), y=norm_cor, size = norm_n, color = Run, group = Run), alpha = 0.2) + geom_line(aes(x=reorder(sample, -norm_cor), y=norm_cor, color = Run, group = Run)) + theme_bw() + theme(axis.text.x = element_text(angle = 90)) + ggtitle("Norm correlation") + xlab("Samples") + ylab("Norm filter")

ggplot(data = processing_reports2 %>% filter(use.col != "Umi3")) +geom_point(aes(x=reorder(sample, -norm_cor_filter), y=norm_cor_filter, size = norm_n_filter, color = Run, group = Run), alpha = 0.2) + geom_line(aes(x=reorder(sample, -norm_cor_filter), y=norm_cor_filter, color = Run, group = Run)) + theme_bw() + theme(axis.text.x = element_text(angle = 90)) + ggtitle("Norm filter correlation") + xlab("Samples") + ylab("Norm cor filter")


ggplot(data = processing_reports2 %>% 
         filter(use.col %in% "Umi1" & dna.crs.is.umi %in% "DnaIsUmiT") %>% 
         group_by(Run) %>%
         filter(all(norm_n_filter >= 10000)) %>%
         ungroup()) +
  geom_point(aes(x=reorder(sample, -norm_cor_filter), y=norm_cor_filter, size = norm_n_filter, color = Run, group = Run), alpha = 0.2) + 
  geom_line(aes(x=reorder(sample, -norm_cor_filter), y=norm_cor_filter, color = Run, group = Run), linetype = "dashed") + 
  facet_grid(~crs_filter) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90), legend.position = "bottom") + 
  guides(color = guide_legend(ncol = 2, byrow = TRUE)) +
  ggtitle("Norm filter correlation") + xlab("Samples") + ylab("Norm cor filter")

ggplot(data = processing_reports2 %>% 
         filter(use.col %in% "Umi1" & dna.crs.is.umi %in% "DnaIsUmiT") %>% 
         group_by(Run) %>%
         filter(all(norm_n_filter >= 10000 & crs_filter %in% "crsfilter20")) %>%
         ungroup()) +
  geom_point(aes(x=reorder(sample, -norm_cor_filter), y=norm_cor_filter, size = norm_n_filter, color = Run, group = Run), alpha = 0.2) + 
  geom_line(aes(x=reorder(sample, -norm_cor_filter), y=norm_cor_filter, color = Run, group = Run), linetype = "dashed") + 
  #facet_grid(~crs_filter) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90), legend.position = "bottom") + 
  guides(color = guide_legend(ncol = 2, byrow = TRUE)) +
  ggtitle("Norm filter correlation") + xlab("Samples") + ylab("Norm cor filter")

##
LibASeqDesign <- read_csv("~/cluster_mount/project/SCG4SYN/LibraryDesign/230731_MetaDesignFiles/meta_files/LibASeqDesign_meta_info.csv")
LibHSeqDesign <- read_csv("~/cluster_mount/project/SCG4SYN/LibraryDesign/230731_MetaDesignFiles/meta_files/LibHSeqDesign_meta_info.csv")
LibHSeqDesignCtrls <- read_csv("~/cluster_mount/project/SCG4SYN/LibraryDesign/230731_MetaDesignFiles/meta_files/LibHSeqDesign_controls_meta_info.csv")

LibASeqMerge <- LibASeqDesign[,c("SeqID", "TF.Trp53.affinity","TF.Trp53.nrepeats","TF.Trp53.orientation","Universal.Spacing")]
colnames(LibASeqMerge)[1] <- "ControlFor"


Trp53meta <- inner_join(LibHSeqDesignCtrls, LibASeqMerge, by = "ControlFor" )
colnames(Trp53meta)[2] <- "CRS"
Trp53meta <- Trp53meta %>%
  mutate(TF.Trp53.affinity_nrepeats = TF.Trp53.affinity* TF.Trp53.nrepeats)

#1
report_libH_ps <- system("find ~/cluster_mount/project/SCG4SYN/Pipelines/MPRA/3_Report_screen/output/4_HSC_LibraryH_minP/filtered/Umi1_DnaIsUmiT_crsfilter20_minDnaUmi1_minRnaUmi0_Reads1_InputF_PseudoT_pseudoDna_pseudoValue1/ -iname 'final.crs_report.umi.rda'", intern = T)
#2
report_libH_ps2 <- system("find ~/cluster_mount/project/SCG4SYN/Pipelines/MPRA/3_Report_screen/output/4_HSC_LibraryH_minP/filtered/Umi1_DnaIsUmiT_crsfilter20_minDnaUmi1_minRnaUmi0_Reads1_InputF_PseudoF/ -iname 'final.crs_report.umi.rda'", intern = T)

#for libH:
processing_report_ps <- lapply(report_libH_ps, function(x) {
  load(x)
  return(final.crs)
})

processing_report_ps2 <- lapply(report_libH_ps2, function(x) {
  load(x)
  return(final.crs)
})


processing_report_ps <- processing_report_ps[[1]]
processing_report_ps2 <- processing_report_ps2[[1]]
Trp53_df <- inner_join(Trp53meta, processing_report_ps, by = "CRS")
Trp53_df2 <- inner_join(Trp53meta, processing_report_ps2, by = "CRS")

ddply(Trp53_df, "clusterID", summarise, R2 = cor(TF.Trp53.affinity* TF.Trp53.nrepeats, mean.norm)^2)
ddply(Trp53_df2, "clusterID", summarise, R2 = cor(TF.Trp53.affinity* TF.Trp53.nrepeats, mean.norm)^2)

qplot(x = TF.Trp53.affinity_nrepeats, y = mean.norm,data = Trp53_df, color = TF.Trp53.affinity) + facet_wrap(~clusterID, scales="free_y") + ggtitle("PseudoT")
qplot(x = TF.Trp53.affinity_nrepeats, y = mean.norm,data = Trp53_df2, color = TF.Trp53.affinity) + facet_wrap(~clusterID, scales="free_y") + ggtitle("PseudoF")
#cor(Trp53_df2$TF.Trp53.affinity_nrepeats,Trp53_df2$mean.norm)

Trp53_df3 <- inner_join(Trp53_df, Trp53_df2, by = c("CRS", "clusterID"))
#cor(Trp53_df3$mean.norm.x,Trp53_df3$mean.norm.y)

qplot(x = mean.norm.x, y = mean.norm.y, data = Trp53_df3, color = TF.Trp53.affinity.x) + facet_wrap(~clusterID, scales="fixed") + xlab("Mean norm PseudoT") + ylab("Mean norm PseudoF")
qplot(x = mean.norm.x, y = mean.norm.y, data = Trp53_df3, color = TF.Trp53.affinity.x) + facet_wrap(~clusterID, scales="free") + xlab("Mean norm PseudoT") + ylab("Mean norm PseudoF")

ggplot(data = Trp53_df3) + geom_point(aes(x = mean.norm.x, y = mean.norm.y, color = TF.Trp53.affinity.x, shape = as.factor(Universal.Spacing.x))) + facet_wrap(~clusterID, scales="free") + theme_bw() + xlab("Mean norm PseudoT") + ylab("Mean norm PseudoF")




###

report_libH_ps3 <- system("find ~/cluster_mount/project/SCG4SYN/Pipelines/MPRA/3_Report_screen/output/4_HSC_LibraryH_minP/filtered/Umi1_DnaIsUmiT_crsfilter20_minDnaUmi1_minRnaUmi0_Reads1_InputF_PseudoT_pseudoDna_pseudoValue1/ -iname 'CRS_summary_report.rda'", intern = T)
report_libH_ps3 <- system("find ~/cluster_mount/project/SCG4SYN/Pipelines/MPRA/3_Report_screen/output/4_HSC_LibraryH_minP/filtered/Umi1_DnaIsUmiT_crsfilter20_minDnaUmi1_minRnaUmi0_Reads1_InputF_PseudoF/ -iname 'CRS_summary_report.rda'", intern = T)


processing_report_ps3 <- lapply(report_libH_ps3, function(x) {
  load(x)
  return(CRS.all)
})

processing_report_ps3 <- lapply(processing_report_ps3, function(x) { x[x$pass %in% T,]})

heatscatter(processing_report_ps3[[1]]$norm.2, processing_report_ps3[[1]]$norm.1, main = "normalized", cor=T, method = "pearson")
heatscatter(log10(processing_report_ps3[[1]]$DNA.norm.1), processing_report_ps3[[1]]$norm.1 , main = "normalized", cor=T, method = "pearson")
heatscatter(log10(processing_report_ps3[[1]]$RNA.norm.1), processing_report_ps3[[1]]$norm.1 , main = "normalized", cor=T, method = "pearson")

heatscatter(processing_report_ps3[[1]]$mean.norm, processing_report_ps3[[1]]$DNA.norm.1, main = "normalized", cor=T, method = "pearson")
heatscatter(processing_report_ps3[[1]]$mean.norm, processing_report_ps3[[1]]$RNA.norm.1, main = "normalized", cor=T, method = "pearson")


heatscatter(CRS.all.filtered$norm.1,CRS.all.filtered$norm.2, main = "normalized", cor=T, method = "pearson")
heatscatter(log2(CRS.all.filtered$RNA.1+1),log2(CRS.all.filtered$RNA.2+1), main = "normalized", cor=T, method = "pearson")
