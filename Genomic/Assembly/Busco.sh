#!/bin/bash
###############################################################################
# Script: Busco.sh
# Description:
#   Run BUSCO to evaluate the completeness of a genome assembly.
#
#   To use:
#     - Edit the variables below to point to your assembly, output directory, lineage, and sample name.
###############################################################################

# -------------------------
# === USER VARIABLES ===
# -------------------------
# Sample name
SAMPLE="sample_name"

# Path to your genome assembly
ASSEMBLY="/path/to/your/assembly.fasta"

# Output directory for BUSCO results
OUTPUT_DIR="./busco_results"

# Path to BUSCO lineage dataset
LINEAGE="/path/to/busco/lineages/fungi_odb10"

# Number of threads to use
THREADS=10

# -------------------------
# Run BUSCO
# -------------------------
echo "=== Running BUSCO on $SAMPLE ==="

busco -i "$ASSEMBLY" \
      -o "$SAMPLE" \
      -l "$LINEAGE" \
      -m genome \
      --cpu "$THREADS" \
      -f \
      --out "$OUTPUT_DIR"

echo "BUSCO analysis completed. Results are in $OUTPUT_DIR/$SAMPLE"
