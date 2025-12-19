#!/bin/bash
#SBATCH --job-name=processPRO_seqdata
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=50G
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
module load genometools
seqBiasDir=/home/FCAM/rmukherjee/seqOutBias/target/release
filesDir=/home/FCAM/rmukherjee/compartmentModel_project_datasets/Jonkers_et_al_2014/fastqfiles
UMI_length=0
mouseGenomeIndex=/home/FCAM/rmukherjee/Spring2022Rotation/mm39_bt/mm39
mouserDNAIndex=/home/FCAM/rmukherjee/Spring2022Rotation/mouse_rDNA_bt2/mouse_rDNA
name=placeholder
annotation_prefix=/home/FCAM/rmukherjee/Spring2022Rotation/genome/annotations/Mus_musculus.GRCm39.105
aThalianaGenomeIndex=/home/FCAM/rmukherjee/arabidopsisThalianaGenome/AThaliana_GCF1735
cores=8
read_size=50
parentMouseGenomeDir=/home/FCAM/rmukherjee/Spring2022Rotation/genome
tallymer=${parentMouseGenomeDir}/mm39.tal_${read_size}.gtTxt.gz
table=${parentMouseGenomeDir}/mm39_${read_size}.4.2.2.tbl
aThalianaGenomeDir=/home/FCAM/rmukherjee/arabidopsisThalianaGenome
aThalianaTallymer=${aThalianaGenomeDir}/GCF_000001735.4.tal_${read_size}.gtTxt.gz
aThalianaTable=${aThalianaGenomeDir}/GCF_000001735.4_${read_size}.4.2.2.tbl

cd ${filesDir}
gunzip ${name}.fastq.gz

cutadapt --cores=${cores} -m $((UMI_length+2)) -O 1 -a TGGAATTCTCGGGTGCCAAGG ${name}.fastq \
	                            -o ${name}_noadap.fastq --too-short-output ${name}_short.fastq > ${name}_cutadapt.txt
## initialize the metrics file
#echo -e  "value\texperiment\tthreshold\tmetric" > ${name}_QC_metrics.txt
#PE1_total=$(wc -l ${name}.fastq | awk '{print $1/4}')
#PE1_w_Adapter=$(wc -l ${name}_short.fastq | awk '{print $1/4}')
#AAligation=$(echo "scale=2 ; $PE1_w_Adapter / $PE1_total" | bc)
#echo -e "$AAligation\t$name\t0.80\tAdapter/Adapter" >> ${name}_QC_metrics.txt
seqtk seq -L $((UMI_length+10)) ${name}_noadap.fastq > ${name}_noadap_trimmed.fastq
#
########################
####### since UMI length is zero
##### - no fqdedup as UMI is not present
##### - no reverse complementing for GRO-seq as `opposite adapters on the 3 and 5 ends for GRO`.
mv ${name}_noadap_trimmed.fastq ${name}_processed_UMI_${UMI_length}.fastq
#########################

## bowtie log showed highest alignment with both b and e removed
#seqtk trimfq -b ${UMI_length} -e ${UMI_length} ${name}_dedup.fastq | seqtk seq -r - > ${name}_processed_UMI_${UMI_length}.fastq

##bowtie2 -p $((cores)) -x ${aThalianaGenomeIndex} -U ${name}_processed_UMI_${UMI_length}.fastq 2>${name}_bowtie2_${UMI_length}_aThaliana.log | samtools view -b - | \
#	                       samtools sort - -o ${name}_${UMI_length}_aThaliana.bam
##parse bowtie log for alignment
#filename=${name}_bowtie2_${UMI_length}_aThaliana.log # logfile
#totalReads=$(wc -l ${name}_processed_UMI_${UMI_length}.fastq | awk '{ print $1/4 }')
#alignmentRate=$(grep "overall alignment rate" $filename | awk -F"%" '{ print $1 }')
#alignmentFraction=$(echo "scale=2 ; ${alignmentRate} / 100" | bc)
#spikeInReads=$(echo "scale=0; ${alignmentFraction} * ${totalReads}" | bc)
#mousealignedPercent=$(grep "overall alignment rate" ${name}_bowtie2_${UMI_length}_human.log | awk -F"%" '{ print $1 }')
#mouseReads=$(echo "scale=0; ${mousealignedPercent} * ${totalReads} * 0.01" | bc)
#echo -e "${alignmentFraction}\t$name\t0.30\tThaliana Alignment Rate" > ${name}_metrics_spikeIn.txt
#echo -e "spikeInReads\tmouseReads\ttotalReads" > ${name}_metrics_spikeIn.txt
#echo -e "${name}\t${spikeInReads}\t${mouseReads}\t${totalReads}" >> ${name}_metrics_spikeIn.txt

#
#bowtie2 -p $((cores)) -x ${mouseGenomeIndex} -U ${name}_processed_UMI_${UMI_length}.fastq 2>${name}_bowtie2_${UMI_length}_mouse.log | samtools view -b - | \
#	                       samtools sort - -o ${name}_${UMI_length}_mouse.bam
##parse bowtie log for alignment
#filename=${name}_bowtie2_${UMI_length}_mouse.log # logfile
#alignmentRate=$(grep "overall alignment rate" $filename | awk -F"%" '{ print $1 }')
#alignmentFraction=$(echo "scale=2 ; ${alignmentRate} / 100" | bc)
#echo -e "${alignmentFraction}\t$name\t0.80\tAlignment Rate" >> ${name}_QC_metrics.txt

# thaliana reads - collect all reads not aligning to mouse genome
bowtie2 -p $((cores)) --maxins 1000 -x ${mouseGenomeIndex} -U ${name}_processed_UMI_${UMI_length}.fastq  2>${name}_bowtie2_non_mouseDNA.log |  samtools sort -n - | samtools fastq -f 0x4 - > ${name}.non.mouse.fastq

# align to ribosomal reads and remove them
bowtie2 -p $((cores)) -x ${mouserDNAIndex} -U ${name}_processed_UMI_${UMI_length}.fastq 2>${name}_bowtie2_rDNA.log | \
	                    samtools sort -n - | samtools fastq -f 0x4 - > ${name}.non.rDNA.fastq
reads=$(wc -l ${name}.non.rDNA.fastq | awk '{print $1/4}')

filename=${name}_bowtie2_rDNA.log # logfile
rDNA_alignment=$(grep "overall alignment rate" $filename | awk -F"%" '{ print $1 }')
rDNA_alignmentFraction=$(echo "scale=2 ; ${rDNA_alignment} / 100" | bc)
echo -e "$rDNA_alignmentFraction\t$name\t0.10\trDNA Alignment Rate" >> ${name}_QC_metrics.txt

bowtie2 -p $((cores)) --maxins 1000 -x ${mouseGenomeIndex} -U ${name}.non.rDNA.fastq  2>${name}_bowtie2_non_rDNA.log | samtools view -b - | samtools sort - -o ${name}.nonrDNA.bam
# parse bowtie log to find number of non-rDNA reads mapping to mouse genome.
filename=${name}_bowtie2_non_rDNA.log # logfile
non_rDNA_reads=$(wc -l ${name}.non.rDNA.fastq | awk '{ print $1/4 }')
non_rDNA_alignment=$(grep "overall alignment rate" $filename | awk -F"%" '{ print $1 }')
nonrDNAalignmentFraction=$(echo "scale=2 ; ${non_rDNA_alignment} / 100" | bc)
non_rDNA_ReadCount=$(echo "scale=0 ; ( ${non_rDNA_alignment} * ${non_rDNA_reads} ) / 100" | bc)
echo -e "${nonrDNAalignmentFraction}\t$name\t0.80\tAlignment Rate" >> ${name}_QC_metrics.txt
#
### complexity and theoretical read depth
#fqComplexity -i ${name}_noadap_trimmed.fastq
#PE1_total=$(wc -l ${name}.fastq | awk '{print $1/4}')
#PE1_noadap_trimmed=$(wc -l ${name}_noadap_trimmed.fastq | awk '{print $1/4}')
#factorX=$(echo "scale=2 ; $PE1_noadap_trimmed / $PE1_total" | bc)
#echo fraction of reads that are not adapter/adapter ligation products or below 10 base inserts
#echo $factorX
###calculate PE1 deduplicated reads
#PE1_dedup=$(wc -l ${name}_dedup.fastq | awk '{print $1/4}')
###divide
#factorY=$(echo "scale=2 ; $non_rDNA_ReadCount / $PE1_dedup" | bc)
##re-run with factors
#fqComplexity -i ${name}_noadap_trimmed.fastq -x $factorX -y $factorY
#
bowtie2 -p $((cores)) -x ${aThalianaGenomeIndex} -U ${name}.non.mouse.fastq  2>${name}_bowtie2_${UMI_length}_aThaliana.log | samtools view -b - | \
	                               samtools sort - -o ${name}_${UMI_length}_aThaliana.bam
##parse bowtie log for alignment
filename=${name}_bowtie2_${UMI_length}_aThaliana.log # logfile
totalReads=$(wc -l ${name}_processed_UMI_${UMI_length}.fastq | awk '{ print $1/4 }')
alignmentRate=$(grep "overall alignment rate" $filename | awk -F"%" '{ print $1 }')
alignmentFraction=$(echo "scale=2 ; ${alignmentRate} / 100" | bc)
spikeInReads=$(echo "scale=0; ${alignmentFraction} * ${totalReads}" | bc)
mousealignedPercent=$(grep "overall alignment rate" ${name}_bowtie2_${UMI_length}_human.log | awk -F"%" '{ print $1 }')
mouseReads=$(echo "scale=0; ${mousealignedPercent} * ${totalReads} * 0.01" | bc)
#echo -e "${alignmentFraction}\t$name\t0.30\tThaliana Alignment Rate" > ${name}_metrics_spikeIn.txt
echo -e "exprName\tspikeInReads\tmouseReads\ttotalReads" > ${name}_metrics_spikeIn.txt
echo -e "${name}\t${spikeInReads}\t${mouseReads}\t${totalReads}" >> ${name}_metrics_spikeIn.txt

##TBD: convert to bigWig and BED6 using seqOutBias
${seqBiasDir}/seqOutBias scale $table ${name}.nonrDNA.bam --no-scale --stranded --bed-stranded-positive \
	                --bw=${name}.bigWig --bed=${name}.bed \
			--tail-edge --read-size=${read_size} --tallymer=$tallymer

grep -v "random" ${name}_not_scaled.bed | grep -v "chrUn" | grep -v "chrEBV" | sort -k1,1 -k2,2n > ${name}_tmp.txt
mv ${name}_tmp.txt ${name}_not_scaled.bed
# run seqoutBias on aThaliana bam
${seqBiasDir}/seqOutBias scale ${aThalianaTable} ${name}_${UMI_length}_aThaliana.bam --no-scale --stranded --bed-stranded-positive \
	                --bw=${name}_aThaliana.bigWig --bed=${name}_aThaliana.bed \
			 --read-size=${read_size} --tallymer=${aThalianaTallymer}

sort -k1,1 -k2,2n ${name}_aThaliana_not_scaled.bed > ${name}_pause_counts_tmp.bed
mv ${name}_pause_counts_tmp.bed ${name}_aThaliana_not_scaled.bed

##count reads in pause region
mapBed -null "0" -s -a $annotation_prefix.pause.bed -b ${name}_not_scaled.bed | \
	                            awk '$7>0' | sort -k5,5 -k7,7nr | sort -k5,5 -u > ${name}_pause.bed
#discard anything with chr and strand inconsistencies
join -1 5 -2 5 ${name}_pause.bed $annotation_prefix.bed | \
	                            awk '{OFS="\t";} $2==$8 && $6==$12 {print $2, $3, $4, $1, $6, $7, $9, $10}' | \
				                                awk '{OFS="\t";} $5 == "+" {print $1,$2+480,$8,$4,$6,$5} $5 == "-" {print $1,$7,$2 - 380,$4,$6,$5}' | \
								                            awk  '{OFS="\t";} $3>$2 {print $1,$2,$3,$4,$5,$6}' > ${name}_pause_counts_body_coordinates.bed

sort -k1,1 -k2,2n ${name}_pause_counts_body_coordinates.bed > ${name}_pause_counts_tmp.bed
mv ${name}_pause_counts_tmp.bed ${name}_pause_counts_body_coordinates.bed

#column ten is Pause index
mapBed -null "0" -s -a ${name}_pause_counts_body_coordinates.bed \
	        -b ${name}_not_scaled.bed | awk '$7>0' | \
		                            awk '{OFS="\t";} {print $1,$2,$3,$4,$5,$6,$7,$5/100,$7/($3 - $2)}' | \
					                                awk '{OFS="\t";} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$8/$9}' > ${name}_pause_body.bed
pause_index.R ${name}_pause_body.bed
#
mapBed -null "0" -s -a $annotation_prefix.introns.bed \
	                    -b ${name}_not_scaled.bed | awk '$7>0' | \
			awk '{OFS="\t";} {print $1,$2,$3,$5,$5,$6,$7,($3 - $2)}' > ${name}_intron_counts.bed

mapBed -null "0" -s -a $annotation_prefix.no.first.exons.named.bed \
	        -b ${name}_not_scaled.bed | awk '$7>0' | \
		        awk '{OFS="\t";} {print $1,$2,$3,$4,$4,$6,$7,($3 - $2)}' > ${name}_exon_counts.bed

exon_intron_ratio.R ${name}_exon_counts.bed ${name}_intron_counts.bed

gzip ${name}.fastq

#cat *_QC_metrics.txt | awk '!x[$0]++' > project_QC_metrics.txt 
#plot_all_metrics.R project_QC_metrics.txt Estrogen_treatment_PRO



