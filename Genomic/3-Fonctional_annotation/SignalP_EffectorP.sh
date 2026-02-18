#!/bin/bash
###############################################################################
# Script: SignalP_EffectorP.sh
# Description:
#   Run SignalP and EffectorP to predict secreted effectors from a proteome.
#   Generates a table with protein ID, length, and cysteine count.
#
#   To use:
#     - Edit the variables below to point to your proteome, output directory, and sample name.
###############################################################################

# -------------------------
# === USER VARIABLES ===
# -------------------------
# Sample name
SAMPLE="sample_name"

# Path to your proteome (FASTA format)
PROTEOME="/path/to/"$SAMPLE"_proteins.fa"

# Output directory for results
OUTPUT_DIR="./Effectors_results"

# Path to EffectorP installation
EFFECTORP_DIR="/path/to/EffectorP-3.0"

mkdir -p "$OUTPUT_DIR"

# -------------------------
# Run SignalP to detect secreted proteins
# -------------------------
echo "=== Running SignalP on $SAMPLE ==="
signalp -fasta "$PROTEOME" \
        -org euk \
        -format short \
        -prefix "$OUTPUT_DIR/Fungi_"

SUMMARY_FILE="$OUTPUT_DIR/Fungi__summary.signalp5"

awk '$2=="SP(Sec/SPI)" {print $1}' "$SUMMARY_FILE" > "$OUTPUT_DIR/secreted_ids.txt"

# -------------------------
# Run EffectorP on secreted proteins
# -------------------------
echo "=== Running EffectorP ==="
seqtk subseq "$PROTEOME" "$OUTPUT_DIR/secreted_ids.txt" > "$OUTPUT_DIR/secreted_proteins.fasta"

cd "$EFFECTORP_DIR"

python3 EffectorP.py -f -i "$OUTPUT_DIR/secreted_proteins.fasta" \
    -o "$OUTPUT_DIR/EffectorP_results.txt" \
    -E "$OUTPUT_DIR/Effectors.fasta" \
    -N "$OUTPUT_DIR/NonEffectors.fasta"

# -------------------------
# Generate summary table: Protein ID | Length | Cysteine count
# -------------------------
echo "=== Generating summary table ==="
awk '/^>/{id=substr($0,2); getline seq; cys=gsub(/C/, "", seq); print id"\t"length(seq)"\t"cys}' \
    "$OUTPUT_DIR/Effectors.fasta" > "$OUTPUT_DIR/Effectors_table.tsv"

echo "Effectors prediction completed. Results are in $OUTPUT_DIR"
