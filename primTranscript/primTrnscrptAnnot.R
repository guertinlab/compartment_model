# for bigWig library
# install.packages(c("devtools", "NMF", "pracma", "RColorBrewer"))
#devtools::install_github("andrelmartins/bigWig", subdir="bigWig")
#devtools::install_github("WarrenDavidAnderson/genomicsRpackage/primaryTranscriptAnnotation")
library(bigWig)
library(NMF)
library(dplyr)
library(pracma)
library(RColorBrewer)
#library(primaryTranscriptAnnotation)
# source primaryTranscriptAnnotation package files here
source("/home/FCAM/rmukherjee/Summer2022Rotation/genomicsRpackage-master/primaryTranscriptAnnotation/R/gene_ann.R")
source("/home/FCAM/rmukherjee/Summer2022Rotation/genomicsRpackage-master/primaryTranscriptAnnotation/R/map_tu.R")
# load transcript, firstExon files
#annotDir="/home/FCAM/rmukherjee/Summer2022Rotation/annotations"
annotDir="../data"
# import data for first exons, annotate, and remove duplicate transcripts
fname = sprintf("%s/gencode.mm39.firstExon.bed", annotDir)
dat0 = read.table(fname,header=F,stringsAsFactors=F)
names(dat0) = c('chr', 'start', 'end', 'gene', 'xy', 'strand')
dat0 = unique(dat0)
gencode.firstExon = dat0
# import data for all transcripts, annotate, and remove duplicate transcripts
fname = sprintf("%s/gencode.mm39.transcript.bed", annotDir)
dat0 = read.table(fname,header=F,stringsAsFactors=F)
names(dat0) = c('chr', 'start', 'end', 'gene', 'xy', 'strand')
dat0 = unique(dat0)
gencode.transcript = dat0

# location of bigWig files
#parentDir="/home/FCAM/rmukherjee/Spring2022Rotation/proseqfiles/merged_sorted_bigWig_files"
parentDir="../data"
plus.file = sprintf("%s/Macrophages_merged_sorted_plus_PE1.bigWig", parentDir)
minus.file = sprintf("%s/Macrophages_merged_sorted_minus_PE1.bigWig", parentDir) 
bw.plus = load.bigWig(plus.file)
bw.minus = load.bigWig(minus.file)

# (using gencode transcript bed) get intervals for furthest TSS and TTS +/- interval
# stitch together the lowest start and highest end of each gene
largest.interval.bed = get.largest.interval(bed=gencode.transcript)
save(largest.interval.bed, file="gencode.largest.interval.Rdata")
# (use gencode transcript bed) get read counts and densities for each annotated gene
transcript.reads = read.count.transcript(bed=gencode.transcript,
bw.plus=bw.plus, bw.minus=bw.minus)
save(transcript.reads, file="gencode.transcript.reads.Rdata")

# decide on cutoffs for discarding reads - evaluate count and density distributions
pdf("read_count_density.pdf")
par(mfrow=c(1,2))
hist(log(transcript.reads$density), breaks=200,
col="black",xlab="log read density",main="")
hist(log(transcript.reads$counts), breaks=200,
col="black",xlab="log read count",main="")
dev.off()

# create cutoffs and discard gene with low expression
den.cut = -3.5 # density < exp(-3.5)
cnt.cut = 3    # count < exp (3)
# plot the graph
pdf("read_count_density_thresholdLine.pdf")
par(mfrow=c(1,2))
hist(log(transcript.reads$density), breaks=200,
col="black",xlab="log read density",main="")
abline(v=den.cut, col="red")
hist(log(transcript.reads$counts), breaks=200,
col="black",xlab="log read count",main="")
abline(v=cnt.cut, col="red")
dev.off()
# remove "unexpressed" genes
# get indices based on density and cnt threshold
ind.cut.den = which(log(transcript.reads$density) < den.cut)
ind.cut.cnt = which(log(transcript.reads$counts) < cnt.cut)
ind.cut = union(ind.cut.den, ind.cut.cnt) # take union of indices to be removed 
unexp = names(transcript.reads$counts)[ind.cut] # get all unexpressed genes
largest.interval.expr.bed = largest.interval.bed[!(largest.interval.bed$gene %in% unexp),] # get interval bed of all genes which are expressed in sufficient amount
save(largest.interval.expr.bed, file="largest_interval_expr_bed.Rdata")

# select the TSS (defined as having the highest count) for each gene and incorporate these TSSs
# into the largest interval coordinates
bp.range = c(20,120)
cnt.thresh = 5 #these number of counts are needed for TSS to be evaluated
bed.out = largest.interval.expr.bed
# select exon 1 having genes in largest.interval. Probably to enforce same genes?
bed.in = gencode.firstExon[gencode.firstExon$gene %in% bed.out$gene,]
# bp.range is the interval within which the greatest read density is searched for
# the following regions on gene are searched through bigwig.query 
# for - strand, the start and end are set as "end-120", "end-20"
# for + strand, the start and end are set as "start+20" and "start+120" 
# after the max coverage is found,
# the + strand TSS is set to "start"
# whereas - strand TSS is set to "end"
# finally, the gene names from TSS calculated are used to set the same genes in largest 
# interval bed, where genes in plus(+) strand get their "start" set as TSS whereas genes in minus(-) strand get their "end" set as TSS 
# get.TSS also returns an "issues" object which contain all genes whose ends are lower than starts and whose start are less than zero
# the "bed" object returned contains the interval data patched with TSS as start.
gene.TSS = get.TSS(bed.in=bed.in, bed.out=bed.out,
bw.plus=bw.plus, bw.minus=bw.minus,
bp.range=bp.range, cnt.thresh=cnt.thresh)
TSS.gene = gene.TSS$bed
# pre-filter to remove poorly annotated genes
for (substr in c("Gm", "Rik", "AC", "AW")) {
   grep.result = grep(substr, TSS.gene$gene)
   if (length(grep.result) > 0) { 
      TSS.gene = TSS.gene[-grep.result,]
   }
}
# 11681 genes found in total
save(TSS.gene, file="TSS_gene_filtered.Rdata")

# remove overlaps by selecting the gene with the highest density
# (for both plus and minus strand) gene.overlaps searches for genes within the same chromosome which has a start inside other genes. A dataframe is created with gene list and which genes start within them ("has.start.inside") or where they start within other genes ("is.start.inside").
# also, number of cases of such occurence of inner and outer genes are documented for each of the genes. (integer overlap.)
overlap.data = gene.overlaps( bed = TSS.gene )
# the cases from gene.overlaps is used again within remove.overlaps function
# only those genes are kept which have the maximum density within the transcript bed
# the Tss.gene object with rows filtered out are returned
filtered.id.overlaps = remove.overlaps(bed=TSS.gene,
overlaps=overlap.data$cases,
transcripts=gencode.transcript,
bw.plus=bw.plus, bw.minus=bw.minus, by="den")

# now start finding TTS coordinates
# get intervals for TTS evaluation
# the output from the function is an input for get.TTS()
# for plus strand,
# add.to.end is number of bases added to continue searching after gene end annotation
# start of TTS search is set by subtracting fraction.end*(gene_length_from_annotation) from (end+add.to.end)
# for minus strand,
# add.to.end is subtracted from start
# the end of TTS search is set by adding fraction.end*(fraction.end*gene_length) from (start - add.to.end) 
# (for both plus and minus strand) dist.to.start is the number of bases left from gene start of the next gene in annotati
add.to.end = 100000
fraction.end = 0.2
dist.from.start = 50
bed.for.tts.eval = get.end.intervals(bed=filtered.id.overlaps,
add.to.end=add.to.end,
fraction.end=fraction.end,
dist.from.start=dist.from.start)

# identify gene ends
# @param knot.div (in get.TTS) the number of knots for the spline fit is defined as the number of bins covering the gene end divided by this number;
#' increasing this parameter results and a smoother curve
add.to.end = max(bed.for.tts.eval$xy)
knot.div = 40
pk.thresh = 0.05
bp.bin = 50
knot.thresh = 5
cnt.thresh = 5
tau.dist = 50000
frac.max = 1
frac.min = 0.3
# get.TTS tries to figure how to consider the area with the last transcription signal
# as count data is disordered, we need to find proper peaks. For that purpose a spline is fitted with number of knots (the number of polynomial functions.)
# then we search for the maximum value reached by spline. From that peak of spline,
# we move to the last peak of the spline which at a threshold decided by equations in the code.
# we call the x coodirnate of that peak as TTS.
# (all the above process in reverse direction for minus strand.) 
inferred.coords = get.TTS(bed=bed.for.tts.eval, tss=filtered.id.overlaps,
bw.plus=bw.plus, bw.minus=bw.minus,
bp.bin=bp.bin, add.to.end=add.to.end,
pk.thresh=pk.thresh, knot.thresh=knot.thresh,
cnt.thresh=cnt.thresh, tau.dist=tau.dist,
frac.max=frac.max, frac.min=frac.min, knot.div=knot.div)
save(inferred.coords, file="inferred_coords.Rdata")
coords = inferred.coords$bed

# (from section 8, page 27 of primaryTranscriptAnnotation vignette)
# filter for chrom sizes - if any ends are greater than the species chromosome size,
# reduce the end to the size limit
# any negative starts will be set to 0
chr.sizes.file="/home/FCAM/rmukherjee/Spring2022Rotation/genome/mm39.chrom.sizes"
chrom.sizes = read.table(chr.sizes.file,stringsAsFactors=F,header=F)
names(chrom.sizes) = c("chr","size")

coords.filt = chrom.size.filter(coords=coords, chrom.sizes=chrom.sizes)
coords.filt$log
head(coords.filt$bed)
# process output, remove mitochondrial coordinates, and write out data
out = coords.filt$bed %>% select(chr,start,end,gene)
out = cbind(out, c(500), coords.filt$bed$strand) %>% as.data.frame(stringsAsFactors=F)
# create a gene list for use while merging
pta.gene.list = out$gene # gene lists from primary Transcript annotation workflow 
# save the previous data as inferred coords
out = out[-which(out$chr == "chrM"),]
write.table(out,"inferredcoords.bed",quote=F,sep="\t",col.names=F,row.names=F)

# now, begin to merge together a full annotation list.
removed.gene.interval = largest.interval.bed %>% filter(!(gene %in% pta.gene.list | chr == "chrM"))
names(removed.gene.interval) = names(out) # set the colnames equal
removed.gene.interval[c("c(500)")] = 100 # set score to something low as these were discarded
master.annotation.bed = rbind(out, removed.gene.interval)
write.table(master.annotation.bed, "mastercoords.bed", quote=F,sep="\t",col.names=F,row.names=F)

