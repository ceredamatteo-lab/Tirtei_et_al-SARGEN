options(stringsAsFactors=F)

library(dplyr)
library(ggplot2)

LOCAL_DIR <- "" # Set your working directory

setwd(LOCAL_DIR)

# Signatures (CNV) ==== 

db <- readRDS( "Rdata/metadata.rds")

tissue_annotation <- db[, c("bridge_barcode", "tissue_general")]

process_tumor_data <- function(file_path, tumor_name, signature_type) {
  # Read the data
  data <- read.delim2(file_path)
  
  # Calculate summary
  summary <- cbind.data.frame(colSums(data[, 2:ncol(data)] > 0))
  colnames(summary) <- "N"
  summary$total <- nrow(data)
  summary$Percentage_of_samples <- summary$N / summary$total
  summary$Tumor <- tumor_name
  summary$Signature <- rownames(summary)
  
  # Optionally remove certain signatures (like 'SBS96E' for some cases)
  if (signature_type == "EW") {
    summary <- subset(summary, !rownames(summary) == "SBS96E")
  }
  
  return(summary)
}

#### Figure 2G and Figure 3E ====

tumor_types <- c("OS_P", "OS_R", "OS_M","EW_P", "EW_R", "EW_M","RMS_P", "RMS_R", "RMS_M", "SS_R")

# Initialize empty list to store results
final_results <- list()

# Process SBS and ID data
for (tumor in tumor_types) {
  # SBS files (handle by subfolder)
  sbs_file_path <- paste0("Signatures/", tumor, "/COSMIC_SBS96_Activities.txt")
  id_file_path <- paste0("Signatures/", tumor, "/COSMIC_ID83_Activities.txt")
  cnv_file_path <- paste0("Signatures/", tumor, "/COSMIC_CNV48_Activities.txt")
  
  # Process SBS and ID data and store the results
  if (file.exists(sbs_file_path)) {
    
    sbs_summary <- process_tumor_data(sbs_file_path, tumor, "SBS")
    sbs_summary$Cat <- "SBS"  # Add category column for SBS
    final_results[[paste0("SBS_", tumor)]] <- sbs_summary
  }
  
  if (file.exists(id_file_path)) {
    id_summary <- process_tumor_data(id_file_path, tumor, "ID")
    id_summary$Cat <- "ID"  # Add category column for ID
    final_results[[paste0("ID_", tumor)]] <- id_summary
  }
  
  if (file.exists(cnv_file_path)) {
    cnv_summary <- process_tumor_data(cnv_file_path, tumor, "CNV")
    cnv_summary$Cat <- "CNV"  # Add category column for CNV
    final_results[[paste0("CNV_", tumor)]] <- cnv_summary
  }
}

# Combine all tumor summaries into a single data frame
combined_data <- do.call(rbind.data.frame, final_results)

#The signature annotation file was downloaded from https://cancer.sanger.ac.uk/signatures/
anno_signature <- read.csv("Tables/Signature_annotation.csv", sep=",")

f <- merge(combined_data,anno_signature, by = "Signature")

f$Tumor <- factor(f$Tumor, levels = c("OS_P", "OS_R", "OS_M","EW_P", "EW_R", "EW_M","RMS_P", "RMS_R", "RMS_M", "SS_R"))
f$Subclass <- factor(f$Subclass, levels = c("Clock-like signature", "DNA repair deficiency", "Tobacco Exposure","Chemical exposure", "UV Exposure", "Chemotherapy treatment","Deaminase and Editing Enzyme Activity", "Focal LOH", "Chromosomal LOH","Chromosomal alteration", "Structural chromosomal alterations","Ploidy changes","Unknown","Sequencing Artefact"))

# Filter only SBS 

f_sbs <- f %>% filter(Cat %in% c("SBS"))

f_sbs <- f_sbs %>%
  mutate(Tumor_group = sub("_.*", "", Tumor))  # keep only part before '_'

pdf("Figures/Signatures_SBS.pdf", useDingbats = F, width=unit(14, "cm"), height=unit(7, "cm"))
ggplot(f_sbs, aes(x = Tumor, 
                  y = reorder(Signature, Percentage_of_samples), 
                  size = Percentage_of_samples, 
                  fill = N)) +
  scale_fill_gradient(low = "#FFF4BC", high = "#E0115F") +
  geom_point(color = "black", shape = 21) +
  facet_grid(rows = vars(Subclass), cols = vars(Tumor_group), 
             scales = "free", space = "free") +
  theme_minimal() +
  theme(legend.position = "right", 
        legend.direction = "vertical", 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        text = element_text(size = 8),
        strip.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1),
        strip.text.y = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
  ylab("Signatures") +
  xlab("Tumor type")
dev.off()

# Filter only ID 
library(dplyr)

f_id <- f %>% filter(Cat %in% c("ID"))

f_id <- f_id %>%
  mutate(Tumor_group = sub("_.*", "", Tumor))  # keep only part before '_'

f_id <- f_id %>%
  mutate(Tumor_group = sub("_.*", "", Tumor))  # keep only part before '_'

pdf("Figures/Signatures_ID.pdf", useDingbats = F, width=unit(14, "cm"), height=unit(7, "cm"))
ggplot(f_id, aes(x = Tumor, 
                  y = reorder(Signature, Percentage_of_samples), 
                  size = Percentage_of_samples, 
                  fill = N)) +
  scale_fill_gradient(low = "#FFF4BC", high = "#E0115F") +
  geom_point(color = "black", shape = 21) +
  facet_grid(rows = vars(Subclass), cols = vars(Tumor_group), 
             scales = "free", space = "free") +
  theme_minimal() +
  theme(legend.position = "right", 
        legend.direction = "vertical", 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        text = element_text(size = 8),
        strip.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1),
        strip.text.y = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
  ylab("Signatures") +
  xlab("Tumor type")
dev.off()

# Filter only CNV rows

f_cnv <- f %>% filter(Cat == "CNV")

f_cnv <- f_cnv %>%
  mutate(Tumor_group = sub("_.*", "", Tumor))  # keep only part before '_'

pdf("Figures/Signatures_CNV.pdf", useDingbats = F, width=unit(14, "cm"), height=unit(7, "cm"))
ggplot(f_cnv, aes(x = Tumor, 
                  y = reorder(Signature, Percentage_of_samples), 
                  size = Percentage_of_samples, 
                  fill = N)) +
  scale_fill_gradient(low = "#FFF4BC", high = "#E0115F") +
  geom_point(color = "black", shape = 21) +
  facet_grid(rows = vars(Subclass), cols = vars(Tumor_group), 
             scales = "free", space = "free") +
  theme_minimal() +
  theme(legend.position = "right", 
        legend.direction = "vertical", 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        text = element_text(size = 8),
        strip.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1),
        strip.text.y = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
  ylab("Signatures") +
  xlab("Tumor type")
dev.off()

