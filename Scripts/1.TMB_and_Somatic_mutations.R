options(stringsAsFactors=F)

library(ggstatsplot)
library(ggpubr)
library(ggridges)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(plyr)

LOCAL_DIR <- "" # Set your working directory

setwd(LOCAL_DIR)

# 1.0 DNA somatic analysis ====

db <- readRDS("Rdata/metadata.rds")
mut <- readRDS("Rdata/somatic_mutations.rds")

# VAF correction based on the purity

purity = readRDS('Rdata/Purity_ploidy_sequenza_selected.rds')
mut_purity = mut
mut_purity$FREQ_corrected = mut_purity$FREQ/purity$cellularity[1]
mut_purity$FREQ_corrected[which(mut_purity$FREQ_corrected>1)] = 1

saveRDS(mut_purity, "Rdata/somatic_mutations_purity_correction.rds") 

# 1.1 Tumor mutational burden (TMB)  ====

target  = fread('Capture_kits/Agilent_Exome_V6_plus_COSMIC/S07604715_on_target.bed', skip = 1, sep='\t')

exome_target_size_Mb = sum(target[,3]-target[,2])/10^6

stat = ddply(mut, .(sample), summarise
             ,  n_ont=sum(on_target), n_ont_sel=sum(selected[on_target])
             ,  n=length(sample), n_sel=sum(selected)
)
stat$TMB_selected     = stat$n_sel/exome_target_size_Mb

stat <- stat %>%
  separate(sample, into = c("barcode", "date"), sep = "\\.")

stat <- stat %>%
  left_join(
    db %>%
      dplyr::select(bridge_barcode, disease, tissue, tissue_general,
                    age_years, sex, key, site, organ),
    by = c("barcode" = "bridge_barcode")
  )

saveRDS(stat, "Rdata/Tumor_mutational_burden.rds")

####

TMB <- readRDS("Rdata/Tumor_mutational_burden.rds")

selected_subtypes <- c("Osteosarcoma", "Ewing's Sarcoma","Rhabdomyosarcoma", "Synovial Sarcoma")

TMB_selected <- TMB %>%
  filter(disease %in% selected_subtypes)

TMB_summary <- TMB_selected %>%
  filter(disease %in% selected_subtypes) %>%
  ungroup() %>% 
  summarise(
    mean_TMB = round(mean(TMB_selected, na.rm = TRUE), 2),
    median_TMB = round(median(TMB_selected, na.rm = TRUE), 2)
  )

my_comparisons <- combn(selected_subtypes, 2, simplify = FALSE)  # all pairwise comparisons
TMB_selected$disease <- factor(TMB_selected$disease, levels = selected_subtypes)

##### Figure 2A ====

pdf("Figures/Box_plot_TMB_diseases.pdf", useDingbats = F, width=unit(14, "cm"), height=unit(7, "cm"))
ggplot(TMB_selected, aes(x = disease, y = TMB_selected)) +
  geom_boxplot(alpha = 0.7, notch = TRUE, outlier.shape = NA) +
  geom_jitter(
    fill = "black",
    color = "black",
    size = 1.5,
    stroke = 0.5,
    width = 0.2,
    shape = 21
  ) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = my_comparisons,
    step.increase = 0.1,
    hide.ns = TRUE
  ) +
  scale_y_log10() +
  labs(
    title = "TMB by Sarcoma Subtype",
    x = NULL,
    y = "TMB selected"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(hjust = 1)
  )
dev.off()

##### Figure 2B ====

TMB_ewing <- TMB %>% filter(disease == "Ewing's Sarcoma")
TMB_ewing$tissue_general <- factor(TMB_ewing$tissue_general, levels = c("PRIMARY", "RECURRENT", "METASTASIS"))

pdf("Figures/Box_plot_TMB_EW.pdf", useDingbats = F, width=unit(14, "cm"), height=unit(7, "cm"))
ggplot(TMB_ewing, aes(x = tissue_general, y = TMB_selected)) +
  geom_boxplot(alpha = 0.7, notch = TRUE, outlier.shape = NA) +
  
  geom_jitter(
    position = position_jitter(width = 0.2, height = 0),
    shape = 21,
    fill = "black",
    color = "black",
    size = 2,
  ) +
  
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(
                       c("PRIMARY", "RECURRENT"),
                       c("PRIMARY", "METASTASIS"),
                       c("RECURRENT", "METASTASIS")
                     ),
                     hide.ns = TRUE) +
  scale_y_log10() +
  labs(
    title = "TMB in Ewing Sarcoma by Tissue Type",
    x = NULL,
    y = "Tumor Mutational Burden"
  ) +
  theme_minimal()
dev.off()

# 1.2 Maftools ====

## Create the summary of the mutations using maftools

library(maftools)
source("Scripts/config/mutations.R")

mut_os=subset(mut, disease=="Osteosarcoma")

FP = unique(c(subset(mut_os, false_cancer_gene)$symbol, 'TTN','MUC4','MUC5B','MUC12'))

filtered_mut_os <- mut_os[!mut_os$Gene.refGene %in% FP, ]

maf = ToMAF(
  subset(filtered_mut_os, selected)
  , refBuild ="hg19"
  , MAFobj = T
  , basename = 'Rdata/Maf_OS_selected'
)

pdf("Figures/summary_maftools.pdf", useDingbats = F, width=unit(14, "cm"), height=unit(7, "cm"))
plotmafSummary(maf = maf
               , rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = F, showBarcodes = T
               , top = 25)
dev.off

