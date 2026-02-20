#!/bin/bash
###############################################################################
# Script: gffread.sh
# Description:
#   Extract CDS, transcripts and/or protein sequences from a genome annotation
#   using gffread.
#
#   To use:
#     - Edit the variables below.
#     - Set extraction options to true or false depending on what you want.
###############################################################################

# -------------------------
# === USER VARIABLES ===
# -------------------------

# Sample prefix for output files
PREFIX="sample_name"

# Genome FASTA file
FASTA="/path/to/your/genome.fasta"

# Annotation file (GFF/GFF3)
GFF="/path/to/your/annotation.gff3"

# Output directory
OUTDIR="./gffread_results"

# Set to true or false
EXTRACT_CDS=true
EXTRACT_TRANSCRIPTS=true
EXTRACT_PROTEINS=true

mkdir -p "$OUTDIR"

# -------------------------
# === RUN GFFREAD ===
# -------------------------

echo "=== Running gffread extraction ==="

# Base gffread command
GFFREAD_CMD="gffread \"$GFF\" -g \"$FASTA\""

# Conditional options
if [ "$EXTRACT_CDS" = true ]; then
    GFFREAD_CMD+=" -x \"$OUTDIR/${PREFIX}_cds.fa\""
    echo "[INFO] CDS extraction enabled"
fi

if [ "$EXTRACT_TRANSCRIPTS" = true ]; then
    GFFREAD_CMD+=" -w \"$OUTDIR/${PREFIX}_transcripts.fa\""
    echo "[INFO] Transcript extraction enabled"
fi

if [ "$EXTRACT_PROTEINS" = true ]; then
    GFFREAD_CMD+=" -y \"$OUTDIR/${PREFIX}_proteins.fa\""
    echo "[INFO] Protein extraction enabled"
fi

# Check that at least one output is requested
if [ "$EXTRACT_CDS" = false ] && \
   [ "$EXTRACT_TRANSCRIPTS" = false ] && \
   [ "$EXTRACT_PROTEINS" = false ]; then
    echo "[WARNING] No extraction selected. Nothing to do."
    exit 1
fi


echo ""
echo "[COMMAND]"
echo "$GFFREAD_CMD"
echo ""

eval "$GFFREAD_CMD"

echo "Extraction completed."
echo "Results are available in: $OUTDIR"
