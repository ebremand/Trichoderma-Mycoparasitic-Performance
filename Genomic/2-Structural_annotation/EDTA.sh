#!/bin/bash
###############################################################################
# Script: EDTA.sh
# Description:
#   Run EDTA to annotate transposable elements in a genome assembly.
#
#   To use:
#     - Edit the variables below to point to your genome assembly and output directory.
###############################################################################

# -------------------------
# === USER VARIABLES ===
# -------------------------

# Sample name
SAMPLE="sample_name"

# Path to genome assembly (FASTA)
GENOME="/path/to/your/assembly.fasta"

# Output directory
OUTPUT_DIR="./EDTA_results"

# Number of threads to use
THREADS=10

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

# -------------------------
# === RUN EDTA ===
# -------------------------

echo "=== Running EDTA on $SAMPLE ==="

EDTA.pl \
  --genome "$GENOME" \
  --species others \
  --step all \
  --anno 1 \
  --threads "$THREADS" \
  --sensitive 1 \
  --evaluate 1 \
  --overwrite 0 \
  --force 1

echo "EDTA annotation completed."
echo "Results are available in: $OUTPUT_DIR"
