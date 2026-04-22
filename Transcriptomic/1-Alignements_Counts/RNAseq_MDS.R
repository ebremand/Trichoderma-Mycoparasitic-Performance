###############################################################################
# Script: RNAseq_MDS.R
# Description:
#   Perform MDS (PCoA) analysis on RNA-seq samples using log-transformed
#   expression values and Bray-Curtis distance.
#   Samples are colored by strain and shaped by timepoint.
###############################################################################

library(tidyverse)
library(vegan)
library(ape)

# --- 1) Load data ---
data <- read.delim(
  "Y:/irhs-FungiSem/Projet TRICHOSEED/5- RNAseq/5.3- Analyse RNAseq 2/2- Alignements/ACP/N1508_net.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# --- 2) Set gene names as rownames ---
data <- column_to_rownames(data, var = colnames(data)[1])

# --- 3) Convert all values to numeric ---
data <- data %>%
  mutate(across(everything(), as.numeric))

# --- 4) Replace NA with 0 ---
data[is.na(data)] <- 0

# --- 5) Remove genes with no expression ---
data <- data[rowSums(data) > 0, ]

# --- 6) Log transformation ---
data_log <- log1p(data)

# --- 7) Transpose (samples as rows) ---
data_t <- as.data.frame(t(data_log))

# --- 8) Compute Bray-Curtis distance matrix ---
dist_matrix <- vegdist(data_t, method = "bray")

# --- 9) Perform MDS (PCoA) ---
pcoa_res <- cmdscale(dist_matrix, k = 2, eig = TRUE)

# Variance explained
pcoa_eig <- round(100 * pcoa_res$eig / sum(pcoa_res$eig), 1)

# --- 10) Build coordinates dataframe ---
mds_coord <- as.data.frame(pcoa_res$points)
colnames(mds_coord) <- c("MDS1", "MDS2")
mds_coord$sample <- rownames(mds_coord)

# --- 11) Metadata annotation ---
mds_coord <- mds_coord %>%
  mutate(
    strain = case_when(
      str_detect(sample, "1237") ~ "I1237",
      str_detect(sample, "3041") ~ "P3041",
      str_detect(sample, "3080") ~ "P3080",
      str_detect(sample, "3116") ~ "P3116",
      str_detect(sample, "1295") ~ "MMS1295",
      str_detect(sample, "1508") ~ "N1508",
      TRUE ~ "Other"
    ),
    timepoint = case_when(
      str_detect(sample, "_a_") ~ "after",
      str_detect(sample, "_b_") ~ "before",
      TRUE ~ "unknown"
    )
  )

# --- 12) Set factor order ---
mds_coord$strain <- factor(
  mds_coord$strain,
  levels = c("I1237", "P3041", "P3080", "P3116", "MMS1295", "N1508")
)

mds_coord$timepoint <- factor(
  mds_coord$timepoint,
  levels = c("before", "after")
)

# --- 13) Plot MDS ---
ggplot(mds_coord, aes(x = MDS1, y = MDS2, color = strain, shape = timepoint)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(group = strain), linetype = 2, alpha = 0.3) +
  scale_color_brewer(palette = "Set2") +
  theme_minimal() +
  labs(
    x = paste0("Axis 1 (", pcoa_eig[1], "%)"),
    y = paste0("Axis 2 (", pcoa_eig[2], "%)"),
    color = "Strains",
    shape = "Timepoint"
  ) +
  theme(
    text = element_text(size = 12),
    panel.grid.minor = element_blank()
  )
