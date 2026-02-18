#!/bin/bash
###############################################################################
# Script: Mash_Tree.sh
# Description:
#   Build a genome distance matrix using Mash, convert it to PHYLIP format with Python,
#   and generate a phylogenetic tree with FastME.
#
#   To use:
#     - Edit the variables below to point to your genome folder and output directory.
###############################################################################

## Variables
GENOMES="/path/to/your/genomes"       # folder containing genome FASTA files
OUTPUT="./mash_results"                # output folder
THREADS=10                             # number of threads to use

# -------------------------
# Preparation
# -------------------------
mkdir -p "$OUTPUT"
cd "$OUTPUT"

# Get all genome assemblies
ASSEMBLIES=(${GENOMES}/*.fna ${GENOMES}/*.fasta)

echo "[INFO] Detected genomes:"
printf '%s\n' "${ASSEMBLIES[@]}"

# -------------------------
# Mash sketches
# -------------------------
echo "[INFO] Creating Mash sketches..."
mash sketch -o "${OUTPUT}/all_genomes" -k 21 -s 10000 "${ASSEMBLIES[@]}"

# -------------------------
# Mash distance matrix
# -------------------------
MASH_DIST_FILE="${OUTPUT}/mash_distances.tab"
echo "[INFO] Calculating Mash distances (full matrix)..."
mash dist "${OUTPUT}/all_genomes.msh" "${OUTPUT}/all_genomes.msh" > "$MASH_DIST_FILE"

# Optional triangular matrix
echo "[INFO] Creating triangular matrix for reference..."
mash triangle "${OUTPUT}/all_genomes.msh" > "${OUTPUT}/mash_triangle.tab"

# -------------------------
# Convert to PHYLIP and run FastME
# -------------------------
echo "[INFO] Converting Mash distances to PHYLIP and generating tree..."

FASTME_INPUT="${OUTPUT}/fastme_input.phy"
PYTHON_SCRIPT="${OUTPUT}/phylip_converter.py"
TREE_OUTPUT="${OUTPUT}/mash_phylogeny.nwk"

cat << EOF > "$PYTHON_SCRIPT"
import sys
import os

lines = [l.strip().split('\t') for l in sys.stdin if l.strip()]

names = sorted(list(set(l[0] for l in lines)))
name_to_index = {name: i for i, name in enumerate(names)}
num_genomes = len(names)

matrix = [[0.0] * num_genomes for _ in range(num_genomes)]

for g1, g2, dist_str, _, _ in lines:
    try:
        i = name_to_index[g1]
        j = name_to_index[g2]
        dist = float(dist_str)
        matrix[i][j] = dist
    except (ValueError, KeyError):
        continue

print(num_genomes)
for i in range(num_genomes):
    genome_name = os.path.basename(names[i]).split('.')[0].replace('_', '')[:10]
    dist_line = " ".join([f"{d:.6f}" for d in matrix[i]])
    print(f"{genome_name:<10} {dist_line}")
EOF

python "$PYTHON_SCRIPT" < "$MASH_DIST_FILE" > "$FASTME_INPUT"
rm -f "$PYTHON_SCRIPT"

echo "[INFO] Running FastME..."
fastme -i "$FASTME_INPUT" -o "$TREE_OUTPUT" -m N

echo "[SUCCESS] Newick tree generated: ${TREE_OUTPUT}"
