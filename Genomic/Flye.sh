



reads="/scratch/ebremand/Scratch/Job/Genomes_Pu_Rs/Advanced_Filtering_Rs/Rhizoctonia_solani_PURE_reads.fastq"
output="/scratch/ebremand/Scratch/Job/Genome_Pu_Rs_2/Rhizoctonia_pure/Flye"
threads=10
genome_size="40m" 

mkdir -p "$output"


echo "=== 1. Lancement de l'assemblage Flye (Reads bruts) ==="
cd "$output"

# Lancement de Flye en utilisant directement le fichier $reads
# Le paramètre --asm-coverage 30 gère la profondeur de données utilisées
flye --nano-raw "$reads" --out-dir "$output/assembly_flye_unfiltered" \
     --genome-size $genome_size \
     --threads $threads \
     --asm-coverage 30

