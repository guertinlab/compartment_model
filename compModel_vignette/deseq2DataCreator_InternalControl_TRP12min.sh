#!/bin/bash
set -e
module load bedtools # specific for Xanadu cluster
proseqfolder=/home/FCAM/rmukherjee/compartmentModel_project_datasets/Jonkers_et_al_2014/fastqfiles
annotation_bed_TRP_sorted=/home/FCAM/rmukherjee/compartmentModel_project_datasets/Jonkers_et_al_2014/primaryTranscriptAnnotation/GenesWithRegionsAffected_TRP_12min_sorted.bed
annotationbed_TRP_internalControls=/home/FCAM/rmukherjee/compartmentModel_project_datasets/Jonkers_et_al_2014/primaryTranscriptAnnotation/InternalControlGenes_TRP_12min_sorted.bed
outputFolder=/home/FCAM/rmukherjee/compartmentModel_project_datasets/Jonkers_et_al_2014/deseq2

printf "entering proseqfolder: %s\n" "${proseqfolder}"
cd ${proseqfolder}
filenames=$(ls V65*rep?_not_scaled.bed | grep -v "FP")
for filename in ${filenames}
do
 name=$(echo $filename | awk -F"_not_scaled.bed" '{print $1}')
 echo -e  "\t${name}" > ${outputFolder}/${name}_gene_counts_Affected_Trp12min.txt
 mapBed -null "0" -s -a ${annotation_bed_TRP_sorted} -b $filename |\
	 awk '{OFS="\t";} {print $4,$7}' >> ${outputFolder}/${name}_gene_counts_Affected_Trp12min.txt
done

#internal spike ins
for filename in ${filenames}
do
 name=$(echo $filename | awk -F"_not_scaled.bed" '{print $1}')
 echo -e  "\t${name}" > ${outputFolder}/${name}_gene_counts_InternalControl_Trp12min.txt
 mapBed -null "0" -s -a ${annotationbed_TRP_internalControls} -b $filename |\
	 awk '{OFS="\t";} {print $4,$7}' >> ${outputFolder}/${name}_gene_counts_InternalControl_Trp12min.txt
done

printf "entering outputFolder: %s\n" "${outputFolder}"
cd ${outputFolder}
paste -d'\t' *_gene_counts_Affected_Trp12min.txt > PRO_TRP_V65_Affected_GeneCounts.txt
paste -d'\t' *_gene_counts_InternalControl_Trp12min.txt >  PRO_TRP_V65_InternalControl_GeneCounts.txt


