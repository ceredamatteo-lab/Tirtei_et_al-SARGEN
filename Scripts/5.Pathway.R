options(stringsAsFactors=F)

library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)
library(stringr)
library(readr)

LOCAL_DIR <- "" # Set your working directory

setwd(LOCAL_DIR)

### The Genes_pathways.csv was dowloaded from Oncogenic Signaling Pathways in The Cancer Genome Atlas. 
# Sanchez-Vega, FranciscoCaesar-Johnson, Samantha J. et al. 
# Cell, Volume 173, Issue 2, 321 - 337.e10

path_genes=read.csv("Tables/Genes_pathways.csv",sep =";")

#### Mutation pathways ====

mut <- readRDS("Rdata/somatic_mutations.rds")
db=readRDS("Rdata/metadata.rds")

mut$tissue_general=db$tissue_general[match(mut$barcode, db$bridge_barcode)]

# Filter mutations for selected diseases
mut_selected <- subset(mut, disease %in% c("Osteosarcoma", "Ewing's Sarcoma", 
                                           "Rhabdomyosarcoma", "Synovial Sarcoma"))

#Count mutations per gene × disease × tissue
gene_counts <- mut_selected %>%
  group_by(Gene.refGene, disease, tissue_general) %>%
  dplyr::summarise(MutationCount = n(), .groups = "drop")

gene_pathway_counts <- gene_counts %>%
  inner_join(path_genes, by = c("Gene.refGene" = "Gene"), relationship = "many-to-many")

annotated_table <- gene_pathway_counts %>%
  inner_join(mut_selected, 
             by = c("Gene.refGene", "disease", "tissue_general"), relationship = "many-to-many")

final_selected <- annotated_table %>%
  dplyr::select(Gene.refGene, disease, tissue_general, key, id, Pathway,ExonicFunc.refGene)

ordered_table <- final_selected %>%
  mutate(
    disease = factor(disease, levels = c("Osteosarcoma", "Ewing's Sarcoma", "Rhabdomyosarcoma", "Synovial Sarcoma")),
    tissue_general = factor(tissue_general, levels = c("PRIMARY", "RECURRENT", "METASTASIS"))
  ) %>%
  arrange(disease, tissue_general, Pathway)

write.csv(ordered_table, "Tables/pathway_snvs_results.csv", row.names = FALSE)

# Sum mutation counts per pathway × disease × tissue
pathway_counts <- gene_pathway_counts %>%
  group_by(Pathway, disease, tissue_general) %>%
  dplyr::summarise(MutationCount = sum(MutationCount), .groups = "drop") %>%
  mutate(Disease_Tissue = paste(disease, tissue_general, sep = "_")) %>%
  dplyr::select(Pathway, Disease_Tissue, MutationCount)

disease_levels <- c("Osteosarcoma", "Ewing's Sarcoma", "Rhabdomyosarcoma", "Synovial Sarcoma")
tissue_levels <- c("PRIMARY", "RECURRENT", "METASTASIS")

combined_levels <- as.vector(sapply(disease_levels, function(d) paste(d, tissue_levels, sep = "_")))

pathway_counts$Disease_Tissue <- factor(pathway_counts$Disease_Tissue, levels = combined_levels)

pathway_matrix_df <- pathway_counts %>%
  pivot_wider(names_from = Disease_Tissue, values_from = MutationCount, values_fill = 0)

pathway_matrix <- as.matrix(pathway_matrix_df[, -1])
rownames(pathway_matrix) <- pathway_matrix_df$Pathway

# Create top and left barplot data 
top_contrib <- pathway_matrix
top_contrib_prop <- sweep(top_contrib, 2, colSums(top_contrib), FUN = "/")

mut_pathway <- mut_selected %>%
  inner_join(path_genes, by = c("Gene.refGene" = "Gene"), relationship = "many-to-many")

mutation_types_pathway <- mut_pathway %>%
  dplyr::group_by(Pathway, ExonicFunc.refGene) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = ExonicFunc.refGene, values_from = n, values_fill = 0)

desired_order <- as.vector(sapply(disease_levels, function(d) paste(d, tissue_levels, sep = "_")))
desired_order <- desired_order[desired_order != "Synovial Sarcoma_METASTASIS"]

pathway_matrix <- pathway_matrix[, desired_order]
rownames(pathway_matrix) <- pathway_matrix_df$Pathway

left_bar_matrix <- as.matrix(mutation_types_pathway[, -1])
rownames(left_bar_matrix) <- mutation_types_pathway$Pathway

left_bar_matrix <- left_bar_matrix[rownames(pathway_matrix), , drop = FALSE]

# Calculate proportions for left barplot
left_bar_prop <- sweep(left_bar_matrix, 1, rowSums(left_bar_matrix), "/")

# Define colors
# Left bar (mutation types)
n_colors_left <- max(length(colnames(left_bar_prop)), 3)
left_colors <- brewer.pal(min(n_colors_left, 8), "Set2")[1:length(colnames(left_bar_prop))]
bar_colors <- setNames(left_colors, colnames(left_bar_prop))

# Top bar (disease_tissue)
n_colors_top <- length(colnames(top_contrib_prop))
top_colors <- brewer.pal(min(n_colors_top, 12), "Set3")[1:n_colors_top]
bar_colors_1 <- setNames(top_colors, colnames(top_contrib_prop))

# Annotations
top_bar <- HeatmapAnnotation(
  Contribution = anno_barplot(t(top_contrib_prop),
                              gp = gpar(fill = bar_colors_1),
                              border = FALSE,
                              height = unit(2, "cm")),
  annotation_name_side = "left"
)

left_bar <- rowAnnotation(
  MutationTypes = anno_barplot(left_bar_prop,
                               gp = gpar(fill = bar_colors),
                               border = FALSE,
                               width = unit(2, "cm")),
  annotation_name_side = "top"
)

# Color scale for heatmap
col_fun <- colorRamp2(c(0, max(pathway_matrix)), c("white", "darkred"))

##### Figure 5E ====
pdf("Figures/Heatmap_pathways_snvs.pdf",useDingbats = F, width = unit(6,"cm"), height = unit(8,"cm"))
Heatmap(
  pathway_matrix,
  name = "Mutation Count",
  col = col_fun,
  #top_annotation = top_bar,
  left_annotation = left_bar,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width, height, gp = gpar(fill = fill, col = NA))
    grid.text(pathway_matrix[i, j], x, y, gp = gpar(fontsize = 8))
  }
)
dev.off()

#### CNA GISTIC =====

extract_info <- function(filename) {
  disease <- case_when(
    str_detect(filename, "EW") ~ "Ewing's Sarcoma",
    str_detect(filename, "OS") ~ "Osteosarcoma",
    str_detect(filename, "RH") ~ "Rhabdomyosarcoma",
    str_detect(filename, "SS") ~ "Synovial Sarcoma",
    TRUE ~ NA_character_
  )
  
  tissue_general <- case_when(
    str_detect(filename, "_P") ~ "PRIMARY",
    str_detect(filename, "_R") ~ "RECURRENT",
    str_detect(filename, "_M") ~ "METASTASIS",
    TRUE ~ NA_character_
  )
  
  return(list(disease = disease, tissue_general = tissue_general))
}

setwd("Results/")

files <- list.files(pattern = "_summary_Gistic\\.csv$")

all_data <- list()

for (f in files) {
  df <- read.csv(f, stringsAsFactors = FALSE)
  
  info <- extract_info(f)
  df$disease <- info$disease
  df$tissue_general <- info$tissue_general
  
  all_data[[f]] <- df
}

combined_df <- bind_rows(all_data)

combined_df <- combined_df %>%
  distinct(Hugo_Symbol, samples, Fract, disease, .keep_all = TRUE)

names(combined_df)[names(combined_df) == "Hugo_Symbol"] <- "symbol"

write_csv(combined_df,"Supplementary_table_gistic.csv")

#Count mutations per gene per disease and tissue
gene_counts_cnv <- combined_df %>%
  group_by(symbol, disease, tissue_general,N) %>%
  dplyr::summarise(MutationCount = N, .groups = "drop")

gene_pathway_counts_cnv <- gene_counts_cnv %>%
  inner_join(path_genes, by = c("symbol" = "Gene"), relationship = "many-to-many")

annotated_table <- gene_pathway_counts_cnv %>%
  inner_join(combined_df, 
             by = c("symbol", "disease", "tissue_general"), relationship = "many-to-many")

final_selected <- annotated_table %>%
  dplyr::select(symbol, disease, tissue_general, samples, cito, Pathway,Variant_Classification)

ordered_table <- final_selected %>%
  mutate(
    disease = factor(disease, levels = c("Osteosarcoma", "Ewing's Sarcoma", "Rhabdomyosarcoma", "Synovial Sarcoma")),
    tissue_general = factor(tissue_general, levels = c("PRIMARY", "RECURRENT", "METASTASIS"))
  ) %>%
  arrange(disease, tissue_general, Pathway)

write.csv(ordered_table, "Tables/pathway_cnv_results.csv", row.names = FALSE)

pathway_counts_cnv <- gene_pathway_counts_cnv %>%
  dplyr::group_by(Pathway, disease, tissue_general) %>%
  dplyr::summarise(MutationCount = sum(N), .groups = "drop") %>%
  dplyr::mutate(Disease_Tissue = paste(disease, tissue_general, sep = "_")) %>%
  dplyr::select(Pathway, Disease_Tissue, MutationCount)

disease_levels <- c("Osteosarcoma", "Ewing's Sarcoma")
tissue_levels <- c("PRIMARY", "RECURRENT", "METASTASIS")

combined_levels <- as.vector(sapply(disease_levels, function(d) paste(d, tissue_levels, sep = "_")))

pathway_counts_cnv$Disease_Tissue <- factor(pathway_counts_cnv$Disease_Tissue, levels = combined_levels)

pathway_matrix_df_cnv <- pathway_counts_cnv %>%
  pivot_wider(names_from = Disease_Tissue, values_from = MutationCount, values_fill = 0)

pathway_matrix_cnv<- as.matrix(pathway_matrix_df_cnv[, -1])
rownames(pathway_matrix_cnv) <- pathway_matrix_df_cnv$Pathway
existing_cols <- intersect(combined_levels, colnames(pathway_matrix_cnv))

pathway_matrix_cnv <- pathway_matrix_cnv[, existing_cols]

# Create top and left barplot data
top_contrib <- pathway_matrix_cnv
top_contrib_prop <- sweep(top_contrib, 2, colSums(top_contrib), FUN = "/")

cnv_pathway <- combined_df %>%
  inner_join(path_genes, by = c("symbol" = "Gene"), relationship = "many-to-many")

cnv_types_pathway <- cnv_pathway %>%
  group_by(Pathway, Variant_Classification) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = Variant_Classification, values_from = n, values_fill = 0)

desired_order <- as.vector(sapply(disease_levels, function(d) paste(d, tissue_levels, sep = "_")))
desired_order <- desired_order[desired_order != "Osteosarcoma_METASTASIS"]
desired_order <- desired_order[desired_order != "Ewing's Sarcoma_RECURRENT"]

pathway_matrix_cnv <- pathway_matrix_cnv[, desired_order]

left_bar_matrix <- as.matrix(cnv_types_pathway[, -1])
rownames(left_bar_matrix) <- cnv_types_pathway$Pathway

left_bar_matrix <- left_bar_matrix[rownames(pathway_matrix_cnv), , drop = FALSE]

left_bar_prop <- sweep(left_bar_matrix, 1, rowSums(left_bar_matrix), "/")

# Define colors
# Left bar (mutation types)
n_colors_left <- max(length(colnames(left_bar_prop)), 3)
left_colors <- brewer.pal(min(n_colors_left, 8), "Set2")[1:length(colnames(left_bar_prop))]
bar_colors <- setNames(left_colors, colnames(left_bar_prop))

bar_colors["Amp"] <- "#E41A1C"   # red
bar_colors["Del"] <- "#377EB8"   # blue

# Top bar (disease_tissue)
n_colors_top <- length(colnames(top_contrib_prop))
top_colors <- brewer.pal(min(n_colors_top, 12), "Set3")[1:n_colors_top]
bar_colors_1 <- setNames(top_colors, colnames(top_contrib_prop))

# Annotations
top_bar <- HeatmapAnnotation(
  Contribution = anno_barplot(t(top_contrib_prop),
                              gp = gpar(fill = bar_colors_1),
                              border = FALSE,
                              height = unit(2, "cm")),
  annotation_name_side = "left"
)

left_bar <- rowAnnotation(
  MutationTypes = anno_barplot(left_bar_prop,
                               gp = gpar(fill = bar_colors),
                               border = FALSE,
                               width = unit(2, "cm")),
  annotation_name_side = "top"
)

col_fun <- colorRamp2(c(0, max(pathway_matrix_cnv)), c("white", "darkred"))

###### Figure 5F =====
pdf("Figures/Heatmap_pathways_cnvs.pdf",useDingbats = F, width = unit(6,"cm"), height = unit(8,"cm"))
Heatmap(
  pathway_matrix_cnv,
  name = "Mutation Count",
  col = col_fun,
  #top_annotation = top_bar,
  left_annotation = left_bar,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width, height, gp = gpar(fill = fill, col = NA))
    grid.text(pathway_matrix_cnv[i, j], x, y, gp = gpar(fontsize = 8))
  }
)
dev.off()

### Mirrored circular barplot ==== 

# Install if needed
# install.packages("geomtextpath")

library(geomtextpath)
library(ggplot2)
library(dplyr)
library(tidyr)

setwd(LOCAL_DIR)

cnv_path=read_csv("Tables/pathway_cnv_results.csv")
mut_path=read_csv("Tables/pathway_snvs_results.csv")
db=readRDS("Rdata/metadata.rds")

names(mut_path)[names(mut_path) == "ExonicFunc.refGene"] <- "mutation"
names(mut_path)[names(mut_path) == "Gene.refGene"] <- "Hugo_Symbol"

setDT(mut_path)

map_vals <- c(
  "Osteosarcoma" = 53,
  "Ewing's Sarcoma" = 39,
  "Rhabdomyosarcoma" = 13,
  "Synovial Sarcoma" = 5
)

mut_path[, N_tot_samples := map_vals[disease]]

mut_counts <- mut_path %>%
  group_by(disease, Pathway, N_tot_samples) %>%
  dplyr::summarise(
    n_samples = n_distinct(key),
    n_genes = n_distinct(Hugo_Symbol),
    genes_list = paste(unique(Hugo_Symbol), collapse = ", "),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    percentage = round(100 * n_samples / N_tot_samples, 2)
  )

df_long <- cnv_path %>%
  separate_rows(samples, sep = ",\\s*")

df_long$key=db$key[match(df_long$samples, db$bridge_barcode)]
df_long$mutation <- "cnv"
df_long <- df_long %>% dplyr::rename(id = Variant_Classification)
df_long <- df_long %>% dplyr::select(-samples,-cito)
names(df_long)[names(df_long) == "symbol"] <- "Hugo_Symbol"

setDT(df_long)

map_vals <- c(
  "Osteosarcoma" = 53,
  "Ewing's Sarcoma" = 38,
  "Rhabdomyosarcoma" = 12,
  "Synovial Sarcoma" = 5
)

df_long[, N_tot_samples := map_vals[disease]]

cnv_counts <- df_long %>%
  group_by(disease, Pathway, N_tot_samples) %>%
  dplyr::summarise(
    n_samples = n_distinct(key),
    n_genes = n_distinct(Hugo_Symbol), # count unique sample ids
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    percentage = round(100 * n_samples / N_tot_samples, 2)
  )

combined_df <- mut_counts %>%
  dplyr::select(disease, Pathway, percentage_snv = percentage) %>%
  full_join(
    cnv_counts %>%
      dplyr::select(disease, Pathway, percentage_cnv = percentage),
    by = c("disease", "Pathway")
  ) %>%
  dplyr::mutate(
    percentage_snv = ifelse(is.na(percentage_snv), 0, percentage_snv),
    percentage_cnv = ifelse(is.na(percentage_cnv), 0, percentage_cnv)
  )

df_long <- combined_df %>%
  pivot_longer(cols = c(percentage_snv, percentage_cnv),
               names_to = "Type", values_to = "Percentage") %>%
  dplyr::mutate(
    Percentage = ifelse(Type == "Percentage_cnv", -Percentage, Percentage),
    Type = recode(Type,
                  "percentage_snv" = "SNV",
                  "percentage_cnv" = "CNV")
  )

df <- df_long %>%
  filter(disease == "Osteosarcoma") #Change disease subtype based on what you would like to plot

# Colors
colors <- c("SNV" = "#E69F00", "CNV" = "#56B4E9")

# Circle levels
circle_levels <- c(25, 50, 75, 100)

df_label <- df %>%
  distinct(Pathway) %>%
  mutate(
    x = seq_along(Pathway),
    label_radius = 110  
  )

df <- df %>%
  mutate(Percentage = round(Percentage, 1))

df_plot <- df %>%
  mutate(
    Percentage = ifelse(Type %in% c("CNV", "Del", "Deletion"), -Percentage, Percentage)
  )

###### Figure 5A-D ====

pdf("Figures/pathway_circos_OS.pdf", useDingbats = F, width = unit(10, "cm"), height = unit(10, "cm"))
ggplot(df_plot, aes(x = Pathway, y = Percentage, fill = Type)) +
  # Main bars
  geom_col(width = 0.9, color = "white") +
  geom_text(
    aes(
      label = paste0(abs(Percentage), "%"),
      y = ifelse(Type == "SNV", Percentage + 5, Percentage - 5)
    ),
    color = "black",
    size = 3,
    fontface = "bold",
    position = position_dodge(width = 0.9)
  ) +
  coord_polar(clip = "off") +
  scale_fill_manual(values = colors) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
  geom_hline(yintercept = circle_levels, linetype = "dotted", color = "grey70") +
  geom_hline(yintercept = -circle_levels, linetype = "dotted", color = "grey70") +
  geom_text(data = data.frame(y = circle_levels),
            aes(x = 0, y = y, label = paste0(y, "%")),
            inherit.aes = FALSE, size = 3, color = "grey30") +
  geom_text(data = data.frame(y = -circle_levels),
            aes(x = 0, y = y, label = paste0(abs(y), "%")),
            inherit.aes = FALSE, size = 3, color = "grey30") +
  scale_y_continuous(limits = c(-100, 100)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(
    title = "Osteosarcoma – SNV vs CNV by Pathway (Mirrored)",
    fill = ""
  )+
  geom_textpath(
    data = df_label,
    aes(x = seq_along(Pathway), y = label_radius, label = Pathway),
    inherit.aes = FALSE,
    size = 3, color = "black",
    vjust = 0
  )
dev.off()