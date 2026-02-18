#!/bin/bash
###############################################################################
# Script: run_quast.sh
# Description:
#   Run QUAST (Quality Assessment Tool for Genome Assemblies) to evaluate
#   the quality of genome assemblies. Provides metrics such as N50, misassemblies,
#   genome fraction, and GC content.
#
#   To use:
#     - Edit the variables below to point to your assembly, reference, and output directory.
###############################################################################

# -------------------------
# === USER VARIABLES ===
# -------------------------
# Path to your genome assembly file
ASSEMBLY="/path/to/your/assembly.fasta"

# Path to your reference genome file (optional, can be same as assembly or leave empty)
REFERENCE="/path/to/your/reference.fasta"

# Output directory for QUAST results
OUTPUT_DIR="./quast_results"

# Number of threads to use
THREADS=10

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# -------------------------
# Run QUAST
# -------------------------
echo "Running QUAST..."
if [ -z "$REFERENCE" ]; then
    quast.py "$ASSEMBLY" -o "$OUTPUT_DIR" --threads "$THREADS"
else
    quast.py "$ASSEMBLY" -o "$OUTPUT_DIR" -r "$REFERENCE" --threads "$THREADS"
fi

echo "QUAST analysis completed. Results are in $OUTPUT_DIR"
