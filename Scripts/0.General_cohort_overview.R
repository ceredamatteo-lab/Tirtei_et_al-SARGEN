options(stringsAsFactors=F)

library(dplyr)
library(ggplot2)
library(forcats)

LOCAL_DIR <- "" # Set your working directory

setwd(LOCAL_DIR)

# 1.0 General cohort overview ====

db <- readRDS("Rdata/metadata.rds")
purity = readRDS('Rdata/Purity_ploidy_sequenza_selected.rds')

db$purity=purity$cellularity[match(db$bridge_barcode, purity$barcode)]
db$purity <- as.numeric(db$purity)

##### Figure 1B ====

# Create bar plot of tissue_general by subtype

plot_data <- db %>%
  dplyr::count(disease, tissue_general) %>%
  group_by(disease) %>%
  dplyr::mutate(percent = n / sum(n) * 100,
         total = sum(n)) %>%
  ungroup() %>%
  mutate(
    tissue_general = factor(tissue_general, levels = c("PRIMARY", "RECURRENT", "METASTASIS")),
    disease = fct_reorder(disease, total, .desc = TRUE)  # Order by total count
  )

plot_data$disease <- fct_rev(plot_data$disease)

pdf("Figures/Fig1_tissue_general.pdf",useDingbats = F, width = unit(6,"cm"), height = unit(8,"cm"))
ggplot(plot_data, aes(x = disease, y = percent, fill = tissue_general)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("PRIMARY" = "#1f78b4",
                               "RECURRENT" = "#1b9e77",
                               "METASTASIS" = "#984ea3")) +
  labs(x = "Subtype", y = "Percent (%)", fill = "Tissue Type") +
  theme_minimal(base_size = 12)
dev.off()

# Create bar plot of tissue_type by subtype

plot_data_tissue <- db %>%
  dplyr::count(disease, tissue_type) %>%
  group_by(disease) %>%
  dplyr::mutate(percent = n / sum(n) * 100,
         total = sum(n)) %>%
  ungroup() %>%
  mutate(
    tissue_type = factor(tissue_type, levels = c("FRESH", "FFPE")),
    disease = fct_reorder(disease, total, .desc = TRUE),
    disease = fct_rev(disease)
  )

pdf("Figures/Fig1_tissue_type.pdf",useDingbats = F, width = unit(6,"cm"), height = unit(8,"cm"))
ggplot(plot_data_tissue, aes(x = disease, y = percent, fill = tissue_type)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c(
    "FRESH" = "#fdaaaa",  
    "FFPE"  = "#fde0e0"   
  )) +
  labs(x = "Subtype", y = "Percent (%)", fill = "Tissue Type") +
  theme_minimal(base_size = 12)
dev.off()

# Create bar plot of purity by subtype

db <- db %>%
  mutate(
    purity_percent = purity * 100,  # Convert to percentage
    purity_bin = case_when(
      purity_percent < 25                    ~ "<25%",
      purity_percent >= 25 & purity_percent < 50 ~ "25-49%",
      purity_percent >= 50 & purity_percent < 75 ~ "50-74%",
      purity_percent >= 75                   ~ ">=75%",
      TRUE ~ NA_character_
    ),
    purity_bin = factor(purity_bin, levels = c("<25%", "25-49%", "50-74%", ">=75%"))
  )

plot_data_purity <- db %>%
  dplyr::count(disease, purity_bin) %>%
  group_by(disease) %>%
  dplyr::mutate(
    percent = n / sum(n) * 100,
    total = sum(n)
  ) %>%
  ungroup() %>%
  mutate(
    disease = fct_reorder(disease, total, .desc = TRUE),
    disease = fct_rev(disease)
  )

pdf("Figures/Fig1_tumor_purity.pdf",useDingbats = F, width = unit(6,"cm"), height = unit(8,"cm"))
ggplot(plot_data_purity, aes(x = disease, y = percent, fill = purity_bin)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c(
    "<25%"   = "#fde0dd",  # Light pink
    "25-49%" = "#fa9fb5",  # Rose
    "50-74%" = "#c51b8a",  # Deep pink
    ">=75%"  = "#7a0177",   # Dark purple
    "NA" = "grey"
  )) +
  labs(x = "Subtype", y = "Percent (%)", fill = "Tumor Purity") +
  theme_minimal(base_size = 12)
dev.off()

# Create bar plot of sex by subtype

db <- db %>%
  mutate(
    sex_clean = factor(sex, levels = c("M", "F"))
  )

plot_data_sex <- db %>%
  dplyr::count(disease, sex_clean) %>%
  group_by(disease) %>%
  dplyr::mutate(
    percent = n / sum(n) * 100,
    total = sum(n)
  ) %>%
  ungroup() %>%
  mutate(
    disease = fct_reorder(disease, total, .desc = TRUE),
    disease = fct_rev(disease)
  )

pdf("Figures/Fig1_sex.pdf",useDingbats = F, width = unit(6,"cm"), height = unit(8,"cm"))
ggplot(plot_data_sex, aes(x = disease, y = percent, fill = sex_clean)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c(
    "M" = "#a6cee3",  # Pastel blue
    "F" = "#fb9a99"   # Pastel red/pink
  )) +
  labs(x = "Subtype", y = "Percent (%)", fill = "Sex") +
  theme_minimal(base_size = 12)
dev.off()

# Create box plot of age distrubution by subtype

desired_order <- c(
  "Osteosarcoma",
  "Ewing's Sarcoma",
  "Rhabdomyosarcoma",
  "Synovial Sarcoma",
  "Alveolar Soft Part Sarcoma",
  "Liposarcoma",
  "Dedifferentiated Liposarcoma",
  "Epithelioid Sarcoma",
  "Infantile Fibrosarcoma",
  "Leiomyosarcoma",
  "Malignant fibrous histiocytoma of bone",
  "Primitive Mesenchimal Myxoid Tumor"
)

db <- db %>%
  mutate(
    disease = factor(disease, levels = rev(desired_order))
  )

pdf("Figures/Fig1_age.pdf",useDingbats = F, width = unit(6,"cm"), height = unit(8,"cm"))
ggplot(db, aes(x = disease, y = age_years, fill = disease)) +
  geom_boxplot() +
  coord_flip() +
  labs(
    x = "Disease Subtype",
    y = "Age",
    title = "Age Distribution by Disease Subtype"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")
dev.off()

