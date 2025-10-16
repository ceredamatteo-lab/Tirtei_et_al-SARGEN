options(stringsAsFactors=F)

library(dplyr)
library(readr) #recommended v  ≥ 2.0.0
library(stringr)
library(purrr)
library(ggstatsplot)
library(ggpubr)

LOCAL_DIR <- "" # Set your working directory

setwd(LOCAL_DIR)

db <- readRDS( "Rdata/metadata.rds")

# Intratumoral heterogeneity (ITH) ==============

data_dir <- "Pyclone/"

files <- list.files(data_dir, pattern = "_cell_prevalence\\.tsv$", full.names = TRUE)

read_file_with_sample <- function(file_path) {
  df <- read_tsv(file_path, show_col_types = FALSE)
  sample_name <- str_extract(basename(file_path), "^[^_]+")  
  df$sample_id <- sample_name
  return(df)
}

pyclone <- map_dfr(files, read_file_with_sample)

# Hierarchical clustering for merging the singletons

merge_singletons_hclust <- function(df, max_dist = 0.05) {
  cluster_summary <- df %>%
    group_by(cluster_id) %>%
    summarise(
      mean_cp = mean(cellular_prevalence, na.rm = TRUE),
      mean_vaf = mean(variant_allele_frequency, na.rm = TRUE),
      n_mutations = n(),
      .groups = "drop"
    )
  
  singletons <- cluster_summary %>% filter(n_mutations == 1)
  
  df <- df %>% mutate(new_cluster_id = as.character(cluster_id))
  
  if (nrow(singletons) == 0) {
    return(df)
  }
  
  if (nrow(singletons) > 1) {
    mat <- singletons %>% select(mean_cp, mean_vaf)
    dist_matrix <- dist(mat)
    hc <- hclust(dist_matrix, method = "average")
    singleton_groups <- cutree(hc, h = max_dist)
    singletons$new_cluster_id <- paste0("singleton_group_", singleton_groups)
  } else {
    singletons$new_cluster_id <- "singleton_group_1"
  }
  
  df <- df %>%
    left_join(singletons %>% select(cluster_id, new_cluster_id),
              by = "cluster_id",
              suffix = c("", "_singleton")) %>%
    mutate(
      new_cluster_id = ifelse(!is.na(new_cluster_id_singleton),
                              new_cluster_id_singleton,
                              new_cluster_id)
    ) %>%
    select(-new_cluster_id_singleton)
  
  return(df)
}

cluster_singletons_merged <- pyclone %>%
  group_by(sample_id) %>%
  group_modify(~merge_singletons_hclust(.x)) %>%
  ungroup()

# Wilcoxon-based merging
threshold_p <- 0.05

merge_clusters <- function(df) {
  clusters <- unique(df$new_cluster_id)  # use new_cluster_id
  merged <- clusters
  
  for (i in 1:(length(clusters)-1)) {
    for (j in (i+1):length(clusters)) {
      c1 <- clusters[i]
      c2 <- clusters[j]
      
      cp1 <- df %>% filter(new_cluster_id == c1) %>% pull(cellular_prevalence)
      cp2 <- df %>% filter(new_cluster_id == c2) %>% pull(cellular_prevalence)
      vaf1 <- df %>% filter(new_cluster_id == c1) %>% pull(variant_allele_frequency)
      vaf2 <- df %>% filter(new_cluster_id == c2) %>% pull(variant_allele_frequency)
      
      if (length(cp1) > 0 & length(cp2) > 0 &
          length(vaf1) > 0 & length(vaf2) > 0) {
        
        p_cp  <- wilcox.test(cp1, cp2, exact = FALSE)$p.value
        p_vaf <- wilcox.test(vaf1, vaf2, exact = FALSE)$p.value
        
        if (!is.na(p_cp) && !is.na(p_vaf) &&
            p_cp > threshold_p && p_vaf > threshold_p) {
          merged[merged == c2] <- c1
        }
      }
    }
  }
  
  df$new_cluster_id <- merged[match(df$new_cluster_id, clusters)]  # update new_cluster_id
  df
}

# Apply Wilcoxon merging after hierarchical merging

cluster_final <- cluster_singletons_merged %>%
  group_by(sample_id) %>%
  group_modify(~merge_clusters(.x)) %>%
  ungroup()


#### ITH score calculation shannon index 

ith_scores <- cluster_final %>%
  group_by(sample_id, new_cluster_id) %>%
  summarise(mean_cp = mean(cellular_prevalence), .groups = "drop") %>%
  group_by(sample_id) %>%
  mutate(p = mean_cp / sum(mean_cp)) %>%
  summarise(
    ITH_shannon = -sum(p * log(p)),
    n_clusters = n()
  )

ith_scores$disease=db$disease[match(ith_scores$sample_id, db$bridge_barcode)]
ith_scores$tissue_general=db$tissue_general[match(ith_scores$sample_id, db$bridge_barcode)]
ith_scores$location=db$organ[match(ith_scores$sample_id, db$bridge_barcode)]

selected_subtypes <- c("Osteosarcoma", "Ewing's Sarcoma","Rhabdomyosarcoma", "Synovial Sarcoma")

ith_scores <- ith_scores %>%
  filter(disease %in% selected_subtypes)

saveRDS(ith_scores,"Rdata/ITH_scores.rds")

##### Fig 4A ====

ith_scores <- readRDS("Rdata/ITH_scores.rds")

df <- ith_scores[ith_scores$disease == "Osteosarcoma", ]

df <- df %>%
  group_by(location) %>%
  filter(n() >= 3) %>%
  ungroup()

selected_locations <- unique(df$location)
my_comparisons <- combn(selected_locations, 2, simplify = FALSE)

pdf("Figures/ITH_location_OS.pdf", width = 8, height = 6)
ggplot(df, aes(x = location, y = ITH_shannon)) +
  geom_boxplot(alpha = 0.7, notch = TRUE, outlier.shape = NA) +
  geom_jitter(
    fill = "black",
    color = "black",
    shape = 21,
    size = 2.5,
    stroke = 0.5,
    width = 0.2
  ) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = my_comparisons,
                     hide.ns = TRUE) +
  labs(
    x = NULL,
    y = "ITH"
  ) +
  theme_minimal()

dev.off()

##### Fig 4B ====

df <- ith_scores[ith_scores$disease == "Ewing's Sarcoma", ]
df$tissue_general <- factor(df$tissue_general, levels = c("PRIMARY", "RECURRENT", "METASTASIS"))

pdf("Figures/ITH_EW_tissue.pdf", width = 8, height = 6)
ggplot(df, aes(x = tissue_general, y = ITH_shannon)) +
  geom_boxplot(alpha = 0.7, notch = TRUE, outlier.shape = NA) +
  
  geom_jitter(
    fill = "black",   # OUTSIDE of aes()
    color = "black",    # optional: defines border color
    size = 1.5,
    stroke = 0.5,
    width = 0.2,
    shape = 21
  ) +
  scale_y_continuous(limits = c(-0.5, 2.2), expand = c(0, 0)) +  # ← Set scale up to 2
  
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(
                       c("PRIMARY", "RECURRENT"),
                       c("PRIMARY", "METASTASIS"),
                       c("RECURRENT", "METASTASIS")
                     ),
                     hide.ns = TRUE) +
  labs(
    x = NULL,
    y = "ITH"
  ) +
  theme_minimal()
dev.off()

