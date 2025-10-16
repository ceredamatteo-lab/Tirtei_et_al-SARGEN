options(stringsAsFactors=F)

# Set up the libraries and environment ----------

library(ggplot2)
library(ggpubr)
library(forcats)
library(maftools)
library(dplyr)
library(survival)
library(survminer)
library(broom)
library(tidyr)

# set-up the environment

LOCAL_DIR <- "" # Set your working directory

setwd(LOCAL_DIR)

db <- readRDS("Rdata/metadata.rds")

# Survival cytobands GISTIC =====

# Extract GISTIC Data for Osteosarcoma primary samples
all.lesions = 'gistic/OS_P/all_lesions.conf_95.txt'
amp.genes   = 'gistic/OS_P/amp_genes.conf_95.txt'
del.genes   = 'gistic/OS_P/del_genes.conf_95.txt'
scores.gis  = 'gistic/OS_P/scores.gistic'

gistic = readGistic(gisticAllLesionsFile = all.lesions,
                    gisticAmpGenesFile = amp.genes,
                    gisticDelGenesFile = del.genes,
                    gisticScoresFile = scores.gis,
                    isTCGA = F )

# Extract GISTIC Data
genes <- gistic@data %>%
  mutate(key = paste0(Variant_Classification, "_", Cytoband))

# Extract Cytoband Summary
cito <- gistic@cytoband.summary %>%
  mutate(key = paste0(Variant_Classification, "_", Unique_Name)) %>%
  filter(qvalues <= 0.1) 

genes = subset(genes, key %in% cito$key)
genes$cito = sapply(strsplit(genes$Cytoband, "\\:"), "[[",2)
genes$N_samples = cito$nSamples[match(genes$key, cito$key)]
genes$total = gistic@summary$summary[which(gistic@summary$ID == "Samples")]
genes$Fract = genes$N_samples/genes$total

# Copy genes dataset
df <- genes

# Keep only rows where cytoband exists in cito$key
cito_cyto <- sapply(strsplit(cito$key, ":"), "[[", 2)
df <- df[df$cito %in% cito_cyto, ]

# Extract cytoband information
df$cito = sapply(strsplit(df$key, "\\:"), "[[",2)

# Map patient IDs and tissue information from db
df$patient_key <- db$key[match(df$Tumor_Sample_Barcode, db$bridge_barcode)]
df$tissue <- db$tissue[match(df$Tumor_Sample_Barcode, db$bridge_barcode)]
#df$tissue_number <- db$tissue_number[match(df$Tumor_Sample_Barcode, db$bridge_barcode)]
#df$patient_key <- paste0(df$patient_id, "-", df$tissue_number)

# Extract cytoband from key
df$cito <- sapply(strsplit(df$key, ":"), "[[", 2)

df_short <- df[, c("Cytoband", "Hugo_Symbol", "patient_key", "cito", "Variant_Classification")]

df_short <- df_short %>%
  # mutate(value = TRUE) %>%
  pivot_wider(names_from = Cytoband, values_from = Variant_Classification, values_fill = list(value = FALSE))  

df_short <- df_short[, c(2:ncol(df_short))] %>%
  group_by(patient_key) %>%
  dplyr::summarise(across(everything(), ~ if (all(is.na(.))) "None" else dplyr::first(na.omit(.)))
  )

names(df_short)[names(df_short) == "patient_key"] <- "key"
names(df_short)[3:length(names(df_short))] <- paste0("OS_cyto_", names(df_short)[3:length(names(df_short))])

surv_data <- left_join(df_short,db, by = "key")
surv_data$survival_from_enrollment_months <- as.numeric(surv_data$survival_from_enrollment_months)

# Ensure Status is numeric (1 = dead, 0 = alive)
surv_data$statusFUP <- as.numeric(surv_data$statusFUP)

# Set up the survial at 24 months (2 years)
cutoff <- 24

# Apply the cutoff to survival time
surv_data$survival_time_censored <- pmin(surv_data$survival_from_enrollment_months, cutoff)  # Cap at 24 months

# Censor patients who survived beyond 24 months
surv_data$survival_status_censored <- ifelse(surv_data$survival_from_enrollment_months > cutoff, 0, surv_data$statusFUP)

surv_data <- surv_data %>%
  mutate(across(starts_with("OS_cyto_"), as.factor))

surv_data <- surv_data %>%
  mutate(across(starts_with("OS_cyto_"), fct_rev))

saveRDS(surv_data, "Rdata/Survival_data_OS_P_GISTIC.rds")

##### Univariate analysis =====

surv_data <- readRDS("Rdata/Survival_data_OS_P_GISTIC.rds")

columns_to_plot <- colnames(surv_data)[grep("OS_cyto_", colnames(surv_data))]

pdf(file='Figures/Survival_K-M_DNA_GISTIC_amp-del_regions_OS_P.pdf', width = 9, height = 9, useDingbats = F)
for (column_name in columns_to_plot){
  data = surv_data
  df_sub <- data[!is.na(data$survival_time_censored) & !is.na(data$survival_status_censored), ]  
  df_sub <- data[!is.na(data$cito), ]
  
  surv_object <- Surv(df_sub$survival_time_censored, df_sub$survival_status_censored)
  fit <- survfit(surv_object ~ df_sub[[column_name]])
  
  cox <- coxph(surv_object ~ df_sub[[column_name]])
  summary_cox <- summary(cox)
  HR <- summary_cox$coefficients[1, "exp(coef)"]
  
  names(fit$strata) <- sapply(strsplit(names(fit$strata), "\\="), "[[", 2)
  
  p <- ggsurvplot(fit, data = df_sub, risk.table = TRUE, 
                  pval = TRUE, 
                  pval.method = TRUE,
                  conf.int = TRUE,
                  title = paste0("Kaplan-Meier Curve for: ", sapply(strsplit(column_name, "_cyto_"), "[[", 1), " - " , sapply(strsplit(column_name, "_cyto_"), "[[", 2)),
                  legend.title = column_name,
                  palette = c("Amp" = "red", "Del" = "blue", "None" = "darkgoldenrod2"),
                  ggtheme = theme_bw(),
                  break.time.by = 4
  ) 
  print(p)
}

dev.off()

### Multivariate survival ======

surv_data <- readRDS("Rdata/Survival_data_OS_P_GISTIC.rds")

multiv_data <- surv_data
names(multiv_data) <- gsub(":", "__", names(multiv_data))
columns_to_plot <- colnames(multiv_data)[grep("OS_cyto", colnames(multiv_data))]

names(multiv_data)[names(multiv_data) %in% columns_to_plot] <- sapply(strsplit(columns_to_plot, "_cyto_"), "[[", 2)
columns_to_plot <- sapply(strsplit(columns_to_plot, "_cyto_"), "[[", 2)

multiv_res <- data.frame()
all_models <- list()

for (i in columns_to_plot){
  
  data = multiv_data
  df_sub <- data[!is.na(data$survival_time_censored) & !is.na(data$survival_status_censored), ]  # Filter out rows with NA in surv_months or status_FUP
  
  df_sub <- data[!is.na(data$cito), ]
  
  surv_object <- Surv(df_sub$survival_time_censored, df_sub$survival_status_censored)
  
  formula <- as.formula(paste("Surv(survival_time_censored, survival_status_censored) ~", i, " + sex + age_group  + stadiation"))
  big_model <- coxph( formula = formula , data = as.data.frame(df_sub))
  tmp <- tidy(big_model, exponentiate = TRUE, conf.int = TRUE)
  
  tmp$Variant <- i
  
  # Adding info on region classification by GISTIC
  tmp$Variant_Classification <- sapply(strsplit(tmp$term[1], i), "[[", 2)
  tmp$n_variant <- sum(df_sub[[i]] %in% c("Amp", "Del"), na.rm = TRUE)
  tmp$tot_samples <- sum(!is.na(df_sub[[i]]))
  
  tmp$Variant_FDR <- p.adjust(tmp$p.value, method = "fdr")
  
  multiv_res <- rbind(multiv_res, tmp)
  
  all_models[[i]] <- big_model
  
}

multiv_res$FDR <- p.adjust(multiv_res$p.value, method = "fdr")
names(multiv_res)[names(multiv_res) == "estimate"] <- "hazard_ratio"
write.csv(multiv_res, file = "Tables/Survival_cox_OS_R_GISTIC_SurvInfo_MULTIV_summary_table.csv")

pdf(file='Figures/Survival_Multivariate_Forestplot_DNA_Osteo_P_GISTIC_SurvInfo_regions.pdf', width = 9, height = 9, useDingbats = F)

j = 1
i = names(all_models)[j]

for (i in names(all_models)){
  data = multiv_data
  column_name = columns_to_plot[19]
  df_sub <- data[!is.na(data$survival_time_censored) & !is.na(data$survival_status_censored), ]  # Filter out rows with NA in surv_months or status_FUP
  
  df_sub <- data[!is.na(data$cito), ]
  
  
  p <- ggforest(all_models[[i]], data = as.data.frame(df_sub), 
                fontsize = 1.2, 
                main = "Forest plot of Hazard Ratios for multivariate Cox model")
  
  print(p)
}

dev.off()