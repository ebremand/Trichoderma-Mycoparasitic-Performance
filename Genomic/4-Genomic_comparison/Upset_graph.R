#!/usr/bin/env Rscript
###############################################################################
# Script: UpSet_Orthogroups.R
# Description:
#   Generate an UpSet plot from OrthoFinder orthogroups to visualize shared
#   orthologues between strains.
#
#   To use:
#     - Edit the variables below to point to your Orthogroups file and the 
#       list of strains you want to include.
#     - Adjust the UpSet plot parameters as needed.
###############################################################################

# -------------------------
# === USER VARIABLES ===
# -------------------------
# Path to the Orthogroups file from OrthoFinder
ORTHO_FILE <- "/path/to/your/Orthogroups.txt"

# List of strains to include in the UpSet plot
STRAINS <- c("Strain1", "Strain2", "Strain3", "Strain4", "Strain5")

# Maximum number of intersections to display in the plot
N_INTERSECTS <- 21

# -------------------------
# === LOAD LIBRARIES ===
# -------------------------
library(UpSetR)
library(tidyverse)

# -------------------------
# === LOAD DATA ===
# -------------------------
df <- read.delim(ORTHO_FILE, check.names = FALSE, header = TRUE)

# Keep only the Orthogroup column + the selected strains
binary_df <- df[, c("Orthogroup", STRAINS)]

# Convert to binary: 1 if present, 0 if absent
binary_df <- binary_df %>%
  mutate(across(all_of(STRAINS), ~ ifelse(. > 0, 1, 0)))

# Remove Orthogroup column for UpSetR input
upset_input <- binary_df[, STRAINS]

# -------------------------
# === GENERATE UPSET PLOT ===
# -------------------------
upset(
  upset_input,
  sets = rev(STRAINS),           # reverse order for plotting
  nsets = length(STRAINS),
  nintersects = N_INTERSECTS,
  order.by = "freq",
  keep.order = TRUE,             # respect the reversed order
  main.bar.color = "#3498db",    # main bar color (blue)
  sets.bar.color = "grey",       # sets bar color
  text.scale = c(1.2, 1.1, 1, 1, 1.5, 1.5), # adjust text size
  mainbar.y.label = "Number of shared orthogroups",
  sets.x.label = "Number of orthogroups per strain"
)
