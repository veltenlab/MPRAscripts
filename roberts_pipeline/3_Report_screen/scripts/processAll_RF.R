#.libPaths("/home/lars/RLibs/R4.2.1/")

setwd("~/cluster_mount/project/SCG4SYN/Pipelines/MPRA/3_Report_screen/scripts/")
require(ggplot2)
require(parallel)
require(LSD)
require(reshape2)
require(ggrepel)
#require(readr)
#require(dplyr)
source("get_CRS_summary_RF.R")

#read non-filtered BC ass
#new_bc_crs_ass <- read.csv(gzfile("/users/lvelten/project/SCG4SYN/Experiments/221118_LibraryA_recloned_association/Analysis/mapped_filter.csv.gz"))

main.folder <- "~/cluster_mount/project/SCG4SYN/Pipelines/MPRA/2_BC_count/output/4_HSC_LibraryH_minP/filtered"
folders <- list.files(main.folder,full.names = T,pattern = "^State")
names(folders) <- list.files(main.folder,full.names = F,pattern = "^State")

#folders <- folders[c(1)]
#folders <- folders[c(5)]
#f <- folders[1]
#folders <- folders[c(5,6)]

### standard filters ####

# BC_CRS <- read.csv(gzfile("/users/lvelten/project/SCG4SYN/Pipelines/MPRA/1_CRS_BC_association/output/LibraryH_minP/LibHSeqDesign_LibHminP_mapped_filtered.csv.gz"))
# barcodeWhitelist <- subset(BC_CRS, MEANMATCHES >= 290 & MEANMATCHES <= 292 & DEVIANTREADS / (READS + DEVIANTREADS) < 0.2)$BARCODE
#barcodeWhiteList <- readRDS("~/fastdata/LibH/barcodeWhitelist.rds")


#Generate experiment for the cell state
experiments  <- mcmapply(get_CRS_summary,folders,names(folders),
                         MoreArgs = list(return_experiments = T,
                                         filter_columns = F, 
                                         normalize.to.input= T,
                                         use.col = "READS",
                                        # pseudocount_type = "DNA",
                                         #pseudocount_value = 1, 
                                         pseudocount = F), 
                         SIMPLIFY = F, mc.cores = 2, mc.preschedule = F)

#saveRDS(experiments, "~/cluster_mount/project/SCG4SYN/Pipelines/MPRA/3_Report_screen/output/4_HSC_LibraryH_minP/filtered/experiments_InputT_PseudoF.rds")
#experiments <- readRDS("~/cluster_mount/project/SCG4SYN/Pipelines/MPRA/3_Report_screen/output/4_HSC_LibraryH_minP/filtered/experiments_InputF_PseudoF.rds")


#Start debugging
#debug(get_CRS_summary)

#Input experiment for this
OUTDIR = "~/cluster_mount/project/SCG4SYN/Pipelines/MPRA/3_Report_screen/output/4_HSC_LibraryH_minP/filtered/Umi1_DnaIsUmiT_crsfilter20_minUmi1"
OUTDIR = "~/cluster_mount/project/SCG4SYN/Pipelines/MPRA/3_Report_screen/output/4_HSC_LibraryH_minP/filtered/Umi1_DnaIsUmiT_crsfilter20_minUmi3"
OUTDIR = "~/cluster_mount/project/SCG4SYN/Pipelines/MPRA/3_Report_screen/output/4_HSC_LibraryH_minP/filtered/Umi1_DnaIsUmiT_crsfilter30_minUmi3"
OUTDIR = "~/cluster_mount/project/SCG4SYN/Pipelines/MPRA/3_Report_screen/output/4_HSC_LibraryH_minP/filtered/Umi3_DnaIsUmiT_crsfilter10_minUmi1"
OUTDIR = "~/cluster_mount/project/SCG4SYN/Pipelines/MPRA/3_Report_screen/output/4_HSC_LibraryH_minP/filtered/Umi3_DnaIsUmiT_crsfilter20_minUmi1"
OUTDIR = "~/cluster_mount/project/SCG4SYN/Pipelines/MPRA/3_Report_screen/output/4_HSC_LibraryH_minP/filtered/Umi1_DnaIsUmiT_crsfilter10_minUmi1"
OUTDIR = "~/cluster_mount/project/SCG4SYN/Pipelines/MPRA/3_Report_screen/output/4_HSC_LibraryH_minP/filtered/Umi1_DnaIsUmiF_crsfilter10_minUmi1"
OUTDIR = "~/cluster_mount/project/SCG4SYN/Pipelines/MPRA/3_Report_screen/output/4_HSC_LibraryH_minP/filtered/Umi1_DnaIsUmiT_crsfilter20_minUmi1_minUmi1"
OUTDIR = "~/cluster_mount/project/SCG4SYN/Pipelines/MPRA/3_Report_screen/output/4_HSC_LibraryH_minP/filtered/Umi1_DnaIsUmiT_crsfilter10_minUmi1_minUmi1"
OUTDIR = "~/cluster_mount/project/SCG4SYN/Pipelines/MPRA/3_Report_screen/output/4_HSC_LibraryH_minP/filtered/Umi1_DnaIsUmiT_crsfilter20_minUmi1_Reads3"
OUTDIR = "~/cluster_mount/project/SCG4SYN/Pipelines/MPRA/3_Report_screen/output/4_HSC_LibraryH_minP/filtered/Umi1_DnaIsUmiT_crsfilter20_minUmi1_Reads10"
OUTDIR = "~/cluster_mount/project/SCG4SYN/Pipelines/MPRA/3_Report_screen/output/4_HSC_LibraryH_minP/filtered/READS_DnaIsUmiT_crsfilter20_minUmi1_Reads3"

#next one
OUTDIR = "~/cluster_mount/project/SCG4SYN/Pipelines/MPRA/3_Report_screen/output/4_HSC_LibraryH_minP/filtered/Umi1_DnaIsUmiT_crsfilter20_minDnaUmi1_minRnaUmi0_Reads1_InputF_PseudoT_pseudoDna_Median"   

source("get_CRS_summary_RF.R")
all.bc  <- mcmapply(get_CRS_summary,folders,names(folders) , experiments,
       MoreArgs = list(outdir = OUTDIR, use.col = "UMI1", return_report = T, return_experiments = F, 
                       dna.crs.is.umi = "Median", rna.crs.is.umi = F, crs_filter = 20, filter_controlGRE=F,save_report=T,min.dna.umi=1,min.rna.umi=0,
                       normalize.to.input= F, pseudocount = T,  pseudocount_type = "DNA", pseudocount_value = 1, filter_columns = T, 
                       bc_crs_ass = "/home/robert/cluster_mount/project/SCG4SYN/Pipelines/MPRA/1_CRS_BC_association/output/LibraryH_minP/LibHSeqDesign_LibHminP_mapped_filtered.csv.gz",
                       barcodeWhitelist = NULL, AssReadsFilter = 1, plot = T), 
       SIMPLIFY = F, mc.cores = 1, mc.preschedule = F)


reports.umi <- lapply(names(all.bc), function(xx) if (is.list(all.bc[[xx]])) all.bc[[xx]][[2]] else NULL)
names(reports.umi) <- names(all.bc)

final.crs <- lapply(names(all.bc)[sapply(all.bc, is.list)], function(x) {
  out <- all.bc[[x]][[1]]
  out$clusterID <- x
  return(out)
})

final.crs <- do.call(rbind, final.crs)
save(final.crs, reports.umi, file = paste0(OUTDIR, "/final.crs_report.umi.rda"))

#saveRDS(experiments, "~/cluster_mount/project/SCG4SYN/Pipelines/MPRA/3_Report_screen/output/4_HSC_LibraryH_minP/filtered/experiments.rds")
#experiments <- readRDS("~/cluster_mount/project/SCG4SYN/Pipelines/MPRA/3_Report_screen/output/4_HSC_LibraryH_minP/filtered/experiments.rds")
#rm(all.bc, final.crs, reports.umi); gc()

reports.umiReads10 <- reports.umi
reports.umiReads3 <- reports.umi



meta_info <- read_csv("~/cluster_mount/project/SCG4SYN/LibraryDesign/230731_MetaDesignFiles/meta_files/LibHSeqDesign_meta_info.csv")
#meta_info
experiments_subset <- experiments$State_5M[[1]]
colnames(experiments_subset)[2] <- "SeqID"
experiments_subset <- inner_join(experiments_subset, meta_info[,c(1:5,138:140)], by= "SeqID")

View(experiments_subset[experiments_subset$UMI11.rna %nin% NA & experiments_subset$UMI11.rna >0, ])
table(experiments_subset$UMI11.rna[experiments_subset$UMI11.rna %nin% NA & experiments_subset$UMI11.rna >0]  )
  rm(experiments_subset)
gc()

sum(experiments$State_1M[[1]][,25] %nin% NA)

m <- lm(log2(CRS.all$RNA.norm.1) ~ log2(CRS.all$DNA.norm.1))
r <- residuals(m)

heatscatter(log2(CRS.all$DNA.norm.1),abs(r), main = "normalized", cor=T, method = "pearson")


lapply(experiments[[1]], function(x) {x[,grepl(pattern = ".[dr]na", x = colnames(x))]})

##
test_df <- Experiments[[1]] %>% group_by(CRS) %>% count(RNA = UMI1.rna >0 , DNA = UMI1.dna >0)
summary(test_df$n)
ggplot(test_df %>% filter(DNA %in% T)) +geom_violin(aes(x = CRS, y = n)) +facet_grid(~RNA+DNA)

test_df %>% group_by(DNA, RNA)  %>% summarize(
  MeanValue = mean(n, na.rm = TRUE),
  MedianValue = median(n, na.rm = TRUE),
  Quant25 = quantile(n, probs = 0.25),
  Quant75 = quantile(n, probs = 0.75),
  MaxValue = max(n, na.rm = TRUE)
)
#
Experiments2 <- lapply(Experiments, function(x) {x[,grepl(pattern = ".dna", x = colnames(x))] <- x[,grepl(pattern = ".dna", x = colnames(x))] %>% 
  mutate(any_na = rowSums(is.na(.)) > 0) %>%
  mutate(across(everything(), ~ ifelse(any_na, 0.001, .))) %>%
  select(-any_na); x})

test2_df <- Experiments2[[1]] %>% group_by(CRS) %>% count(RNA = UMI1.rna >=1 , DNA = UMI1.dna >=1)

test2_df %>% group_by(DNA, RNA)  %>% summarize(
  MeanValue = mean(n, na.rm = TRUE),
  MedianValue = median(n, na.rm = TRUE),
  Quant25 = quantile(n, probs = 0.25),
  Quant75 = quantile(n, probs = 0.75),
  MaxValue = max(n, na.rm = TRUE)
)

#

test_df2 <- melt(test_df, id.vars = c("CRS","n"), variable.name = "Nucleo")
test_df2$value <- ifelse(test_df2$value %in% NA, FALSE, TRUE)  
ggplot(test_df2[c(1:10000, 65323:75323),]) +geom_point(aes(x=CRS, y= n)) +facet_grid(~Nucleo+value)

##
experiments <- experiments[c(1:2)]

for(i in 1:length(reports.umi)) {
  print(reports.umi[[1]][i])
  print(test[[1]][i])
}

summary(final.crs$DNA.1)
