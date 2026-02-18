#!/bin/bash
###############################################################################
# Script: OrthoFinder.sh
# Description:
#   Run OrthoFinder to infer orthologous groups from multiple proteomes.
#
#   To use:
#     - Edit the variables below to point to your proteome folder and output directory.
###############################################################################

# -------------------------
# === USER VARIABLES ===
# -------------------------
# Folder containing all proteomes in .faa format
PROTEOMES="/path/to/your/proteomes_folder"

# Output directory for OrthoFinder results
OUTPUT_DIR="/path/to/your/orthofinder_results"

# Number of threads to use
THREADS=10

# -------------------------
# Run OrthoFinder
# -------------------------
echo "=== Running OrthoFinder on proteomes in $PROTEOMES ==="

# Remove previous output folder to start fresh
rm -rf "$OUTPUT_DIR/Orthofinder"

orthofinder -f "$PROTEOMES" \
            -t "$THREADS" \
            -a "$THREADS" \
            -S diamond \
            -M msa \
            -o "$OUTPUT_DIR/Orthofinder"

echo "OrthoFinder analysis completed. Results are in $OUTPUT_DIR/Orthofinder"
