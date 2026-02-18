###############################################################################
# Script: RNAseq_PCA.R
# Description:
#   Perform PCA on RNA-seq samples using TPM values.
#   Samples are colored by strain and shaped by timepoint.
###############################################################################

library(tidyverse)
library(FactoMineR)
library(factoextra)

# Load TPM table
data <- read.delim("/path/to/your/TPM_table.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Set gene names as rownames and convert to numeric
data <- column_to_rownames(data, var = colnames(data)[1]) %>%
        mutate(across(everything(), as.numeric))

# Replace NA with 0 and remove genes with no variance
data[is.na(data)] <- 0
data <- data[apply(data, 1, var) > 0, ]

# Log-transform and transpose (samples as rows)
data_t <- as.data.frame(t(log1p(data)))

# Metadata for coloring and shapes
samples <- rownames(data_t)
metadata <- data.frame(
  sample = samples,
  strain = case_when(
    str_detect(samples, "1237") ~ "I1237",
    str_detect(samples, "3041") ~ "P3041",
    str_detect(samples, "3080") ~ "P3080",
    str_detect(samples, "3116") ~ "P3116",
    str_detect(samples, "1295") ~ "MMS1295",
    str_detect(samples, "1508") ~ "N1508",
    TRUE ~ "Other"
  ),
  timepoint = case_when(
    str_detect(samples, "_a_") ~ "after",
    str_detect(samples, "_b_") ~ "before",
    TRUE ~ "unknown"
  )
)
rownames(metadata) <- metadata$sample
metadata$strain <- factor(metadata$strain, levels = c("I1237","P3041","P3080","P3116","MMS1295","N1508"))
metadata$timepoint <- factor(metadata$timepoint, levels = c("before","after"))

# PCA
res_pca <- PCA(data_t, scale.unit = TRUE, graph = FALSE)

# Prepare coordinates
ind_coord <- as.data.frame(res_pca$ind$coord)
ind_coord$sample <- rownames(ind_coord)
ind_coord$strain <- metadata$strain
ind_coord$timepoint <- metadata$timepoint

# Plot PCA
ggplot(ind_coord, aes(x = Dim.1, y = Dim.2, color = strain, shape = timepoint)) +
  geom_point(size = 4) +
  stat_ellipse(aes(group = strain), linetype = 2, alpha = 0.3) +
  scale_color_brewer(palette = "Set2") +
  theme_minimal() +
  labs(x = paste0("Dim 1 (", round(res_pca$eig[1,2],1), "%)"),
       y = paste0("Dim 2 (", round(res_pca$eig[2,2],1), "%)"),
       color = "Strains",
       shape = "Timepoint") +
  theme(text = element_text(size = 12))
