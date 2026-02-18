#!/bin/bash
###############################################################################
# Script: Flye.sh
# Description:
#   Run Flye to assemble a genome from raw reads. 
#   Provides options for genome size, number of threads, and coverage.
#
#   To use:
#     - Edit the variables below to point to your reads, output directory, and genome size.
###############################################################################

# -------------------------
# === USER VARIABLES ===
# -------------------------
# Path to your raw reads (FASTQ)
READS="/path/to/your/reads.fastq"

# Output directory for Flye assembly results
OUTPUT_DIR="./flye_results"

# Number of threads to use
THREADS=10

# Estimated genome size (e.g., 40m for 40 Mb)
GENOME_SIZE="40m"

# Optional: maximum assembly coverage (Flye will downsample if needed)
ASM_COVERAGE=30

# Create output directory
mkdir -p "$OUTPUT_DIR"

# -------------------------
# Run Flye
# -------------------------
echo "=== Running Flye assembly ==="
cd "$OUTPUT_DIR"

flye --nano-raw "$READS" \
     --out-dir "$OUTPUT_DIR/assembly_flye_unfiltered" \
     --genome-size "$GENOME_SIZE" \
     --threads "$THREADS" \
     --asm-coverage "$ASM_COVERAGE"

echo "Flye assembly completed. Results are in $OUTPUT_DIR/assembly_flye_unfiltered"
