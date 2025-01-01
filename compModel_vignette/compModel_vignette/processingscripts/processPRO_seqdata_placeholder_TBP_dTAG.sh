#!/bin/bash
#SBATCH --job-name=processPRO_seqdata_placeholder
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=60G
#SBATCH --mail-user=rmukherjee@uchc.edu
#SBATCH -o processPRO_seqdata_placeholder_%j.out
#SBATCH -e processPRO_seqdata_placeholder_%j.err

set -e
# add path to the R files
export PATH=$PATH:/home/FCAM/rmukherjee/Spring2022Rotation/R_files
module load cutadapt;
module load seqtk;
module load bowtie2;
module load bedtools;
module load samtools;
module load bamtools;
module load genometools;

printf "running on %s\n" $(hostname)

filesDir=/home/FCAM/rmukherjee/compartmentModel_project_datasets/price2022_TBP
UMI_length=4
humanGenomeIndex=/home/FCAM/rmukherjee/humanGenome/hg38
#mothGenomeIndex=/home/FCAM/rmukherjee/mothGenome/moth_JQCY02
mothGenomeIndex=/home/FCAM/rmukherjee/mothGenome/moth_GCF_011064685.1_ZJU_Sfru_1_bowtie/moth_GCF_011064685.1_ZJU_Sfru_1
humanrDNAIndex=/home/FCAM/rmukherjee/humanGenome/human_rDNA
name=placeholder
cores=8
read_size=51 # as 51 was the smallest read found in fastq files
parentHumanGenomeDir=/home/FCAM/rmukherjee/humanGenome
tallymer=${parentHumanGenomeDir}/hg38.tal_${read_size}.gtTxt.gz
table=${parentHumanGenomeDir}/hg38_${read_size}.4.2.2.tbl
annotation_prefix=/home/FCAM/rmukherjee/humanGenome/annotations/Homo_sapiens.GRCh38.104

cd ${filesDir}
gunzip ${name}_PE1.fastq.gz
gunzip ${name}_PE2.fastq.gz
cutadapt --cores=${cores} -m $((UMI_length+2)) -O 1 -a TGGAATTCTCGGGTGCCAAGG ${name}_PE1.fastq \
	            -o ${name}_PE1_noadap.fastq --too-short-output ${name}_PE1_short.fastq > ${name}_PE1_cutadapt.txt
cutadapt --cores=${cores} -m $((UMI_length+10)) -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC ${name}_PE2.fastq \
	    -o ${name}_PE2_noadap.fastq --too-short-output ${name}_PE2_short.fastq > ${name}_PE2_cutadapt.txt

PE1_total=$(wc -l ${name}_PE1.fastq | awk '{print $1/4}')
PE1_w_Adapter=$(wc -l ${name}_PE1_short.fastq | awk '{print $1/4}')
AAligation=$(echo "scale=2 ; $PE1_w_Adapter / $PE1_total" | bc)
echo -e  "value\texperiment\tthreshold\tmetric" > ${name}_QC_metrics.txt
echo -e "$AAligation\t$name\t0.80\tAdapter/Adapter" >> ${name}_QC_metrics.txt

seqtk seq -L $((UMI_length+10)) ${name}_PE1_noadap.fastq > ${name}_PE1_noadap_trimmed.fastq
#remove PCR duplicates
fqdedup -i ${name}_PE1_noadap_trimmed.fastq -o ${name}_PE1_dedup.fastq

PE1_noAdapter=$(wc -l ${name}_PE1_dedup.fastq | awk '{print $1/4}')

#pair FASTQ files
fastq_pair -t $PE1_noAdapter ${name}_PE1_dedup.fastq ${name}_PE2_noadap.fastq
flash -q --compress-prog=gzip --suffix=gz ${name}_PE1_dedup.fastq.paired.fq \
	    ${name}_PE2_noadap.fastq.paired.fq -o ${name}
insert_size.R ${name}.hist ${UMI_length}

seqtk trimfq -b ${UMI_length} -e ${UMI_length} ${name}_PE1_dedup.fastq | seqtk seq -r - > ${name}_PE1_processed_UMI_${UMI_length}.fastq
seqtk trimfq -b ${UMI_length} -e ${UMI_length} ${name}_PE2_noadap.fastq | seqtk seq -r - > ${name}_PE2_processed_UMI_${UMI_length}.fastq

## REMOVE reads aligning to moth genome
bowtie2 -p $((cores-2)) -x ${mothGenomeIndex} -U ${name}_PE1_processed_UMI_${UMI_length}.fastq 2> ${name}_bowtie2_moth.log | \
	samtools sort -n - | samtools fastq -f 0x4 - > ${name}_PE1.nonMoth.fastq

## REMOVE reads aligning to human ribosome
bowtie2 -p $((cores-2)) -x ${humanrDNAIndex} -U ${name}_PE1.nonMoth.fastq 2> ${name}_bowtie2_rDNA.log | \
	samtools sort -n - | samtools fastq -f 0x4 - > ${name}_PE1.nonrDNA.fastq
reads=$(wc -l ${name}_PE1.nonrDNA.fastq | awk '{print $1/4}')
# find the pairs for these reads
fastq_pair -t $reads ${name}_PE1.nonrDNA.fastq ${name}_PE2_processed.fastq

## Now, ALIGN to human genome
bowtie2 -p $((cores-2)) --maxins 1000 -x $genome_index --rf -1 ${name}_PE1.nonrDNA.fastq.paired.fq \
	    -2 ${name}_PE2_processed.fastq.paired.fq 2>${name}_bowtie2.log | samtools view -b - | \
	        samtools sort - -o ${name}.nonrDNA.bam

#calculate the total number of rDNA-aligned reads
PE1_prior_rDNA=$(wc -l ${name}_PE1.nonMoth.fastq | awk '{print $1/4}')
PE1_post_rDNA=$(wc -l ${name}_PE1.nonrDNA.fastq | awk '{print $1/4}')
total_rDNA=$(echo "$(($PE1_prior_rDNA-$PE1_post_rDNA))") 

#calculate the total that concordantly align to hg38 and/or rDNA
concordant_pe1=$(samtools view -c -f 0x42 ${name}.nonrDNA.bam)
total=$(echo "$(($concordant_pe1+$total_rDNA))")

#rDNA alignment rate
rDNA_alignment=$(echo "scale=2 ; $total_rDNA / $total" | bc)

echo -e "$rDNA_alignment\t$name\t0.20\trDNA Alignment Rate" >> ${name}_QC_metrics.txt

# MAPPABIlity RATE
map_pe1=$(samtools view -c -f 0x42 ${name}.nonrDNA.bam)
pre_alignment=$(wc -l ${name}_PE1.nonrDNA.fastq.paired.fq | awk '{print $1/4}')
alignment_rate=$(echo "scale=2 ; $map_pe1 / $pre_alignment" | bc)
echo -e "$alignment_rate\t$name\t0.80\tAlignment Rate" >> ${name}_QC_metrics.txt

#Run-on efficiency
#convert to bigWig and BED6
seqOutBias scale $table ${name}.nonrDNA.bam --no-scale --stranded --bed-stranded-positive \
	    --bw=${name}_new.bigWig --bed=${name}_new.bed --out-split-pairends --only-paired \
	        --tail-edge --read-size=$read_size --tallymer=$tallymer

# Complexity and theoretical read depth
#fqComplexity -i ${name}_PE1_noadap_trimmed.fastq
#PE1_total=$(wc -l ${name}_PE1.fastq | awk '{print $1/4}')
#PE1_noadap_trimmed=$(wc -l ${name}_PE1_noadap_trimmed.fastq | awk '{print $1/4}')
#
#factorX=$(echo "scale=2 ; $PE1_noadap_trimmed / $PE1_total" | bc)
#
#echo fraction of reads that are not adapter/adapter ligation products or below 10 base inserts
#echo $factorX
#
##calculate PE1 deduplicated reads
#PE1_dedup=$(wc -l ${name}_PE1_dedup.fastq | awk '{print $1/4}')
#
##divide
#factorY=$(echo "scale=2 ; $concordant_pe1 / $PE1_dedup" | bc)
#
##re-run with factors
#fqComplexity -i ${name}_PE1_noadap_trimmed.fastq -x $factorX -y $factorY

##convert to bigWig and BED6
seqOutBias scale $table ${name}.nonrDNA.bam --no-scale --stranded --bed-stranded-positive \
	    --bw=${name}_new.bigWig --bed=${name}_new.bed --out-split-pairends --only-paired \
	        --tail-edge --read-size=$read_size --tallymer=$tallymer
#Remove chromosomes not in the gene annotation file and sort for use in mapBed
grep -v "random" ${name}_not_scaled_PE1.bed | grep -v "chrUn" | grep -v "chrEBV" | sort -k1,1 -k2,2n > ${name}_tmp.txt
mv ${name}_tmp.txt ${name}_not_scaled_PE1.bed 
#
##count reads in pause region
#mapBed -null "0" -s -a $annotation_prefix.pause.bed -b ${name}_not_scaled_PE1.bed | \
#	    awk '$7>0' | sort -k5,5 -k7,7nr | sort -k5,5 -u > ${name}_pause.bed
#
##discard anything with chr and strand inconsistencies
#join -1 5 -2 5 ${name}_pause.bed $annotation_prefix.bed | \
#	    awk '{OFS="\t";} $2==$8 && $6==$12 {print $2, $3, $4, $1, $6, $7, $9, $10}' | \
#	        awk '{OFS="\t";} $5 == "+" {print $1,$2+480,$8,$4,$6,$5} $5 == "-" {print $1,$7,$2 - 380,$4,$6,$5}' | \
#		    awk  '{OFS="\t";} $3>$2 {print $1,$2,$3,$4,$5,$6}' > ${name}_pause_counts_body_coordinates.bed
#
#sort -k1,1 -k2,2n ${name}_pause_counts_body_coordinates.bed > ${name}_pause_counts_tmp.bed
#mv ${name}_pause_counts_tmp.bed ${name}_pause_counts_body_coordinates.bed
#
##column ten is Pause index
#mapBed -null "0" -s -a ${name}_pause_counts_body_coordinates.bed \
#	    -b ${name}_not_scaled_PE1.bed | awk '$7>0' | \
#	        awk '{OFS="\t";} {print $1,$2,$3,$4,$5,$6,$7,$5/100,$7/($3 - $2)}' | \
#		    awk '{OFS="\t";} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$8/$9}' > ${name}_pause_body.bed
#
#pause_index.R ${name}_pause_body.bed
#
#mapBed -null "0" -s -a $annotation_prefix.introns.bed \
#	    -b ${name}_not_scaled_PE1.bed | awk '$7>0' | \
#	        awk '{OFS="\t";} {print $1,$2,$3,$5,$5,$6,$7,($3 - $2)}' > ${name}_intron_counts.bed
#
#mapBed -null "0" -s -a $annotation_prefix.no.first.exons.named.bed \
#	    -b ${name}_not_scaled_PE1.bed | awk '$7>0' | \
#	        awk '{OFS="\t";} {print $1,$2,$3,$4,$4,$6,$7,($3 - $2)}' > ${name}_exon_counts.bed
#
#exon_intron_ratio.R ${name}_exon_counts.bed ${name}_intron_counts.bed
#
rm ${name}_PE1_short.fastq
rm ${name}_PE2_short.fastq
rm ${name}_PE1_noadap.fastq
rm ${name}_PE2_noadap.fastq
rm ${name}_PE1_noadap_trimmed.fastq
rm ${name}_PE1_dedup.fastq
rm ${name}_PE1_processed_UMI_${UMI_length}.fastq
rm ${name}_PE2_processed_UMI_${UMI_length}.fastq
rm ${name}_PE1_dedup.fastq.paired.fq   
rm ${name}_PE2_noadap.fastq.paired.fq
rm ${name}_PE1_dedup.fastq.single.fq
rm ${name}_PE2_noadap.fastq.single.fq
rm ${name}_PE1.non.rDNA.fastq.paired.fq
rm ${name}_PE1.non.rDNA.fastq.single.fq
#rm ${name}_PE2_processed_UMI_${UMI_length}.fastq.paired.fq
rm ${name}_PE2_processed_UMI_${UMI_length}.fastq.single.fq
rm ${name}.extendedFrags.fastq.gz
rm ${name}.notCombined_1.fastq.gz
rm ${name}.notCombined_2.fastq.gz
gzip ${name}_PE1.fastq
gzip ${name}_PE2.fastq

