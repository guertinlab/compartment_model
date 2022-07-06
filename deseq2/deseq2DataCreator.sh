#!/bin/bash
#SBATCH --job-name=deseq2DataCreator
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=60G
#SBATCH --mail-user=rmukherjee@uchc.edu
#SBATCH -o deseq2DataCreator_%j.out
#SBATCH -e deseq2DataCreator_%j.err

set -e
module load bedtools # sprcific for Xanadu cluster

release=105
#annotation_prefix=/home/FCAM/rmukherjee/Spring2022Rotation/genome/annotations/Mus_musculus.GRCm39.${release}
proseqfolder=/home/FCAM/rmukherjee/Spring2022Rotation/proseqfiles
annotation_bed_pTA=/home/FCAM/rmukherjee/Summer2022Rotation/R_analysis/primTrnscript/mastercoords.bed # from primary Transcript annotation
outputFolder=/home/FCAM/rmukherjee/Summer2022Rotation/proseqFiles
printf "entering proseqfolder: %s\n" "${proseqfolder}"
cd ${proseqfolder}

for filename in *_not_scaled_PE1.bed
do
    name=$(echo $filename | awk -F"_not_scaled_PE1.bed" '{print $1}')
    echo -e  "\t${name}" > ${outputFolder}/${name}_gene_counts_pTA.txt
    coverageBed -counts -s -a ${annotation_bed_pTA} -b $filename | \
        awk '{OFS="\t";} {print $4,$7}' >> ${outputFolder}/${name}_gene_counts_pTA.txt
done

printf "entering outputFolder: %s\n" "${outputFolder}"
cd ${outputFolder}
paste -d'\t' *_gene_counts_pTA.txt > Macrophages_PRO_gene_counts_pTA.txt


