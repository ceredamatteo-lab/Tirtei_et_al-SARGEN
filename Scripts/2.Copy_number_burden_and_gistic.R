options(stringsAsFactors=F)

library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(ggstatsplot)
library(ggpubr)


LOCAL_DIR <- "" # Set your working directory

setwd(LOCAL_DIR)

# Copy number variation (CNV) ==== 

db <- readRDS("Rdata/metadata.rds")
cnv <- readRDS("Rdata/somatic_cnvs_included.rds") 
cnB <- readRDS("Rdata/somatic_cnv_burden.rds") 

cnB_summary <- cnB %>%
  filter(disease %in% selected_subtypes) %>%
  ungroup() %>% 
  summarise(
    mean_cnB = round(mean(cna_burden, na.rm = TRUE), 2),
    median_cnB = round(median(cna_burden, na.rm = TRUE), 2)
  )

# Figure 3 Box plot

cnv_sel <- cnB %>%
  filter(disease %in% c("Rhabdomyosarcoma", "Osteosarcoma","Ewing's Sarcoma","Synovial Sarcoma")) 
cnv_sel$disease <- factor(cnv_sel$disease, levels = c("Synovial Sarcoma", "Rhabdomyosarcoma", "Osteosarcoma","Ewing's Sarcoma"))

selected_subtypes <- c("Synovial Sarcoma","Rhabdomyosarcoma","Osteosarcoma","Ewing's Sarcoma")
my_comparisons <- combn(selected_subtypes, 2, simplify = FALSE)  # all pairwise comparisons
cnv_sel$disease <- factor(cnv_sel$disease, levels = selected_subtypes)

pdf("Figures/Box_plot_cnB.pdf", useDingbats = F, width=unit(14, "cm"), height=unit(7, "cm"))
ggplot(cnv_sel, aes(x = disease, y = cna_burden)) +
  geom_boxplot(alpha = 0.7, notch = TRUE, outlier.shape = NA) +
  geom_jitter(
    fill = "black",
    color = "black",
    size = 2,
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
  scale_x_discrete(limits = selected_subtypes) +
  scale_y_continuous(limits = c(-0.1, 1.6), expand = c(0, 0)) +  # ← Set scale up to 2
  labs(
    title = "CNA burden by Sarcoma Subtype",
    x = NULL,
    y = "CNA burden"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(hjust = 1)
  )
dev.off()


## GISTIC ====

library(maftools)

all.lesions = 'gistic/OS_P/all_lesions.conf_95.txt'
amp.genes   = 'gistic/OS_P/amp_genes.conf_95.txt'
del.genes   = 'gistic/OS_P/del_genes.conf_95.txt'
scores.gis  = 'gistic/OS_P/scores.gistic'

#If it doesn't work, restart RStudio
gistic = readGistic(gisticAllLesionsFile = all.lesions,
                    gisticAmpGenesFile = amp.genes,
                    gisticDelGenesFile = del.genes,
                    gisticScoresFile = scores.gis,
                    isTCGA = F )

genes <- gistic@data %>%
  dplyr:: mutate(key = paste0(Variant_Classification, "_", Cytoband))

cito <- gistic@cytoband.summary %>%
  dplyr::mutate(key = paste0(Variant_Classification, "_", Unique_Name)) %>%
  filter(qvalues <= 0.1) 

# Summarize genes per sample
genes2 <- genes %>%
  dplyr::group_by(Hugo_Symbol, key, Variant_Classification) %>%
  dplyr::summarise(
    N = n_distinct(Tumor_Sample_Barcode),
    samples = paste(unique(Tumor_Sample_Barcode), collapse = ","),
    .groups = "drop"
  ) %>%
  left_join(dplyr::select(cito, key, qvalues, nSamples), by = "key") %>%
  filter(qvalues <= 0.1) %>%
  dplyr::mutate(
    cito = sub(".*:", "", key),
    total = gistic@summary$summary[gistic@summary$ID == "Samples"],
    Fract = nSamples / total
  )

# Filter genes based on cito$key presence
genes = subset(genes, key %in% cito$key)
genes$cito = sapply(strsplit(genes$Cytoband, "\\:"), "[[",2)
genes$N_samples = cito$nSamples[match(genes$key, cito$key)]
genes$total = gistic@summary$summary[which(gistic@summary$ID == "Samples")]
genes$Fract = genes$N_samples/genes$total

# Save results
write.csv(cito, "Results/OS_P_citobands_Gistic.csv", row.names = FALSE)
write.csv(genes, "Results/OS_P_genes_Gistic.csv", row.names = FALSE)
write.csv(genes2, "Results/OS_P_summary_Gistic.csv", row.names = FALSE)

#### Chromplot

pdf("Figures/Gistic_chromPlot_OS_P.pdf", useDingbats = F, width=unit(14, "cm"), height=unit(7, "cm"))
gisticChromPlot(
  gistic = gistic,
  fdrCutOff = 0.1,
  markBands =genes$cito ,
  color = NULL,
  ref.build = "hg19",
  cytobandOffset = 0.1,
  txtSize = 0.8,
  cytobandTxtSize = 0.6,
  maf = NULL,
  mutGenes = NULL,
  y_lims = NULL,
  mutGenesTxtSize = 0.6
)
dev.off()
