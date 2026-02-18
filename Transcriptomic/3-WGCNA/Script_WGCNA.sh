
###############################################################################
## Desciption du script

# Permet d'aligner des données de RNAseq sur un génome ou un transcriptome
# Utilise le logiciel Salmon

###############################################################################

# Variables

fastq="/scratch/ebremand/Scratch/Data/RNAseq_2/*a_Pu*"

transcriptome="/scratch/ebremand/Scratch/Data/References/Transcriptomes/N1508_Pu19_transcripts.fa" # Transcriptome de reference

output="/scratch/ebremand/Scratch/Job/Genomes_T.atroviride_2/4-RNAseq/N1508_Pu19" 
mkdir -p $output

###############################################################################

salmon --version

echo "########## Running salmon index with Transcriptome ##########"

salmon index -t $transcriptome -i $output/Reference_index

mkdir $output/quant


echo "########## Début Salmon quant ##########"

for r1 in $fastq/*_1.fq.gz; do
    prefix=$(basename $r1 _1.fq.gz)
    salmon quant -i $output/Reference_index --libType A -1 $fastq/${prefix}_1.fq.gz -2 $fastq/${prefix}_2.fq.gz -p 20 --seqBias --useVBOpt --validateMappings -o $output/quant/${prefix}
done

#for r1 in $fastq/*R1.fastq; do
#    prefix=$(basename $r1 R1.fastq)
#    salmon quant -i $output/Reference_index --libType A -1 $fastq/${prefix}R1.fastq -2 $fastq/${prefix}R2.fastq -p 20 --seqBias --useVBOpt --validateMappings -o $output/quant/${prefix}
#done



echo "Fin de Salmon quant"
cd $output/quant/
for i in *; do  cp "$i/quant.sf" "$i.sf" ; done # Renomme les fichiers quant en .sf

echo "Création des fichiers TPM et count"
# Créer une table de comptage en count et en TPM
file1=`ls -1 *.sf | head -1` # On récupère le 1er fichier pour extraire la colonne correspondant au nom des gènes
awk -F" " '{print $1}' $file1 > geneID # Extraction du nom des gènes


for i in *.sf # Création du tableau TPM en reprenant les valeurs des colonnes 4 de tous les fichiers
do
  colName=$(echo $(basename $i .${i##*.}))
  awk -F" " '{print $4}' $i > $i.temp
  sed -i 's/TPM/'$colName'/g' $i.temp  
done
paste geneID *.temp > ../TPM.txt # Le fichier TPM.txt est formé avec les TPM et les noms de gènes


for i in *.sf # Création du tableau count en reprenant les valeurs des colonnes 5 de tous les fichiers
do
  colName=$(echo $(basename $i .${i##*.}))
  awk -F" " '{print $5}' $i > $i.temp
  sed -i 's/NumReads/'$colName'/g' $i.temp  
done
paste geneID *.temp > ../counts.txt # Le fichier count.txt est formé avec les count et les noms de gènes

rm *.temp # Suppression des fichiers temporaires
rm geneID

echo "Lancement de multiqc"
cd ../
multiqc .

