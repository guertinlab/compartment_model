library(bigWig)
library(NMF)
library(dplyr)
library(pracma)
library(RColorBrewer)
library(dplyr)

# source primaryTranscriptAnnotation package files here
source("/home/FCAM/rmukherjee/Summer2022Rotation/genomicsRpackage-master/primaryTranscriptAnnotation/R/gene_ann.R")
source("/home/FCAM/rmukherjee/Summer2022Rotation/genomicsRpackage-master/primaryTranscriptAnnotation/R/map_tu.R")

# modified from get.TSS()
get.counts.around.potential.TSS <- function(bed.in=NULL, bed.out=NULL, bw.plus=NULL, bw.minus=NULL,
					                       bp.range=NULL, cnt.thresh=NULL, low.read=FALSE, by="cnt"){
     final.df = c()
	  # error handling. note that input and output genes must be matched
	  if( length(unique(bed.in$gene)) != length(unique(bed.out$gene)) ){stop("input genes do not match output genes")}
  if( setequal(unique(bed.in$gene),unique(bed.out$gene)) == FALSE ){stop("input genes do not match output genes")}
    # loop through each unique gene
    # map reads a range for each annotated TSS
    bed0 = bed.out
    bed = bed.in
      genes = unique(bed$gene)
      tss = rep(0,length(genes))
        for(ii in 1:length(genes)){
	ind.gene = which(bed$gene == genes[ii])
    strand.gene = unique(bed$strand[ind.gene])[1]
    chr.gene = unique(bed$chr[ind.gene])[1]
        # get all sets of intervals around each annotated TSS
        # sort upstream to downstream
        if(strand.gene == "-"){
		      starts = unique(bed$end[ind.gene])
          ranges = data.frame(chr=chr.gene, start=starts-bp.range[2],
			                        end=starts-bp.range[1], stringsAsFactors=FALSE)
	        ranges = ranges[order(ranges$end,decreasing=TRUE),]
	      }
        if(strand.gene == "+"){
		      starts = unique(bed$start[ind.gene])
	      ranges = data.frame(chr=chr.gene, start=starts+bp.range[1],
				                    end=starts+bp.range[2], stringsAsFactors=FALSE)
	            ranges = ranges[order(ranges$start),]
	          }

	    # get counts and densities for each interval of interest
	    if(strand.gene == "-"){ counts = abs(bed.region.bpQuery.bigWig(bw.minus, ranges)) }
	    if(strand.gene == "+"){ counts = abs(bed.region.bpQuery.bigWig(bw.plus, ranges)) }
	        dens = counts / (ranges$end - ranges$start)

	    # select the best TSS, select the most upstream for ties
	    if(by == "cnt"){ind.max = which(counts == max(counts))[1]}
	    if(by == "den"){ind.max = which(dens == max(dens))[1]}
           final.df = rbind(final.df, c(genes[ii], counts[ind.max], ranges[ind.max, ], strand.gene))
	}
      final.df = as.data.frame(final.df)
      return(final.df)
}
# load transcript, firstExon files
annotDir="/home/FCAM/rmukherjee/Summer2022Rotation/annotations"
outputDir = "/home/FCAM/rmukherjee/compartmentModel_project_datasets/Jonkers_et_al_2014/primaryTranscriptAnnotation"

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
# chromosome sizes
chrom.sizes = read.table("/home/FCAM/rmukherjee/Spring2022Rotation/genome/mm39.chrom.sizes",stringsAsFactors=F,header=F)
names(chrom.sizes) = c("chr","size")

bigwigDir = "/home/FCAM/rmukherjee/compartmentModel_project_datasets/Jonkers_et_al_2014/fastqfiles"
#plus.file = sprintf("%s/TRP_FP_control_merged_plus.bigWig", bigwigDir)
#minus.file = sprintf("%s/TRP_FP_control_merged_minus.bigWig", bigwigDir)
#plus.file = sprintf("%s/V65_CNTRL_All_plus_normalized.bigWig", bigwigDir)
#minus.file = sprintf("%s/V65_CNTRL_All_minus_normalized.bigWig", bigwigDir)
plus.file = sprintf("%s/V65_cntrl_plus_normalized.bigWig", bigwigDir)
minus.file = sprintf("%s/V65_cntrl_minus_normalized.bigWig", bigwigDir)
bw.plus = load.bigWig(plus.file)
bw.minus = load.bigWig(minus.file)
cntrl.rep1.plus = load.bigWig(sprintf("%s/V65_cntrl_rep1_plus.bigWig", bigwigDir))
cntrl.rep1.minus = load.bigWig(sprintf("%s/V65_cntrl_rep1_minus.bigWig", bigwigDir))
cntrl.rep2.plus = load.bigWig(sprintf("%s/V65_cntrl_rep2_plus.bigWig", bigwigDir))
cntrl.rep2.minus = load.bigWig(sprintf("%s/V65_cntrl_rep2_minus.bigWig", bigwigDir))

# (using gencode transcript bed) get intervals for furthest TSS and TTS +/- interval
# stitch together the lowest start and highest end of each gene
largest.interval.bed = get.largest.interval(bed=gencode.transcript)
saveRDS(largest.interval.bed, file=sprintf("%s/gencode.largest.interval.rds", outputDir))
# (use gencode transcript bed) get read counts and densities for each annotated gene
transcript.reads = read.count.transcript(bed=gencode.transcript,
					 bw.plus=bw.plus, bw.minus=bw.minus)
saveRDS(transcript.reads, file = sprintf("%s/gencode.transcript.reads.rds", outputDir))
# decide on cutoffs for discarding reads - evaluate count and density distributions
pdf(sprintf("%s/read_count_density.pdf", outputDir))
par(mfrow=c(1,2))
hist(log(transcript.reads$density), breaks=200,
          col="black",xlab="log read density",main="")
hist(log(transcript.reads$counts), breaks=200,
          col="black",xlab="log read count",main="")
dev.off()
# create cutoffs and discard gene with low expression
den.cut = -7.2 # density < exp(-7.5)
cnt.cut = 0.3    # count < exp (0)
# plot the graph
pdf(sprintf("%s/read_count_density_thresholdLine.pdf", outputDir))
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
saveRDS(largest.interval.expr.bed, file=sprintf("%s/largest_interval_expr_bed.rds", outputDir))

# select the TSS (defined as having the highest count) for each gene and incorporate these TSSs
# into the largest interval coordinates
bp.range = c(20,160) # estimated using plotMetaProfileOfBigwigs.R
cnt.thresh = 0.9 # chose this because the Stat4 correct annotation was at 0.1805. Also doesn't make sense to use this threshold as most reads are really low 
#original count threshold was five. Took few cases with read-depth normalized values and computed the threshold
## normalization factor
# control_rep1 = 0.2495727
# control_rep2 = 0.3610368
# the nromalized count calculated using the following cases
#  control_rep1 control_rep2 normalizedRead
#1            5            0      0.6239319
#2            0            5      0.9025920
#3            5            5      1.5265239
#4            2            3      0.7911280
#5            3            2      0.7353959
#6            1            4      0.8468600
#7            4            1      0.6796639
bed.out = largest.interval.expr.bed
# select exon 1 having genes in largest.interval. Probably to enforce same genes?
bed.in = gencode.firstExon[gencode.firstExon$gene %in% bed.out$gene,]
saveRDS(bed.in, sprintf("%s/V65.Jonkers.bed.in.rds", outputDir))
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

##################### EXERCISE FOR estimating the counts around TSS 
#counts.df = get.counts.around.potential.TSS(bed.in=bed.in, bed.out=bed.out,
#					                       bw.plus=bw.plus, bw.minus=bw.minus,
#							                          bp.range=bp.range, cnt.thresh=cnt.thresh)
#counts.df = as.data.frame(counts.df)
#counts.df[,2] =as.numeric(counts.df[,2])
#saveRDS(counts.df, file=sprintf("%s/max.counts.around.potential.tss.rds", outputDir))
## control rep1
#counts.df.cntrl.rep1 = get.counts.around.potential.TSS(bed.in=bed.in, bed.out=bed.out,
#					                       bw.plus=cntrl.rep1.plus, bw.minus=cntrl.rep1.minus,
#							                          bp.range=bp.range, cnt.thresh=cnt.thresh)
#counts.df.cntrl.rep1 = as.data.frame(counts.df.cntrl.rep1)
#colnames(counts.df.cntrl.rep1) = c("gene", "max.count", "chr", "start", "end", "strand")
#counts.df.cntrl.rep1[,2] =as.numeric(counts.df.cntrl.rep1[,2])
#saveRDS(counts.df.cntrl.rep1, file=sprintf("%s/cntrl.rep1.max.counts.potential.tss.rds", outputDir))
## control rep2
#counts.df.cntrl.rep2 = get.counts.around.potential.TSS(bed.in=bed.in, bed.out=bed.out,
#					                       bw.plus=cntrl.rep2.plus, bw.minus=cntrl.rep2.minus,
#							                          bp.range=bp.range, cnt.thresh=cnt.thresh)
#counts.df.cntrl.rep2 = as.data.frame(counts.df.cntrl.rep2)
#counts.df.cntrl.rep2[,2] =as.numeric(counts.df.cntrl.rep2[,2])
#colnames(counts.df.cntrl.rep2) = c("gene", "max.count", "chr", "start", "end", "strand")
#saveRDS(counts.df.cntrl.rep2, file=sprintf("%s/cntrl.rep2.max.counts.potential.tss.rds", outputDir))
#############################

gene.TSS = get.TSS(bed.in=bed.in, bed.out=bed.out, 
		   bw.plus=bw.plus, bw.minus=bw.minus,
		   bp.range=bp.range, cnt.thresh=cnt.thresh)
TSS.gene = gene.TSS$bed
saveRDS(TSS.gene, file=sprintf("%s/TSS_gene_unfiltered.cnt.thresh.%0.2f.rds", outputDir, cnt.thresh))
# pre-filter to remove poorly annotated genes
for (substr in c("Gm", "Rik", "AC", "AW")) {
	           grep.result = grep(substr, TSS.gene$gene)
   if (length(grep.result) > 0) {
	                    TSS.gene = TSS.gene[-grep.result,]
         }
}
saveRDS(TSS.gene, file=sprintf("%s/TSS_gene_filtered.cnt.thresh.%0.2f.rds", outputDir, cnt.thresh))
# remove overlaps by selecting the gene with the highest density
# (for both plus and minus strand) gene.overlaps searches for genes within the same chromosome which has a start inside other genes. A dataframe is created with gene list and which genes start within them ("has.start.inside") or where they start within other genes ("is.start.inside").
# also, number of cases of such occurence of inner and outer genes are documented for each of the genes. (integer overlap.)
overlap.data = gene.overlaps( bed = TSS.gene )
# the cases from gene.overlaps is used again within remove.overlaps function
# only those genes are kept which have the maximum density within the transcript bed
# the Tss.gene object with rows filtered out are returned
filtered.id.overlaps = remove.overlaps(bed=TSS.gene,
				       overlaps=overlap.data$cases,	                                  
				       transcripts=gencode.transcript,bw.plus=bw.plus, bw.minus=bw.minus, by="den")
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
				     add.to.end=add.to.end, fraction.end=fraction.end,						                                            dist.from.start=dist.from.start)
saveRDS(bed.for.tts.eval, file=sprintf("%s/bed.for.tts.eval.cnt.thresh.%0.2f.rds", outputDir, cnt.thresh))
# identify gene ends
# @param knot.div (in get.TTS) the number of knots for the spline fit is defined as the number of bins covering the gene end divided by this number;
#' increasing this parameter results and a smoother curve
add.to.end = max(bed.for.tts.eval$xy)
knot.div = 40
pk.thresh = 0.05
bp.bin = 50
knot.thresh = 5
Cnt.thresh = 4 # same threshold as for get.TSS
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
			     pk.thresh=pk.thresh, knot.thresh=knot.thresh,cnt.thresh=Cnt.thresh, tau.dist=tau.dist,
			     frac.max=frac.max, frac.min=frac.min, knot.div=knot.div)
saveRDS(inferred.coords, file=sprintf("%s/inferred_coords.cnt.thresh.%0.2f.rds", outputDir, cnt.thresh))
coords = inferred.coords$bed
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
write.table(out,sprintf("%s/inferredcoords.cnt.thresh.%0.2f.bed", outputDir, cnt.thresh),quote=F,sep="\t",col.names=F,row.names=F)
# now, begin to merge together a full annotation list.
removed.gene.interval = largest.interval.bed %>% filter(!(gene %in% pta.gene.list | chr == "chrM"))
names(removed.gene.interval) = names(out) # set the colnames equal
removed.gene.interval[c("c(500)")] = 100 # set score to something low as these were discarded
master.annotation.bed = rbind(out, removed.gene.interval)
write.table(master.annotation.bed, sprintf("%s/mastercoords.cnt.thresh.%0.2f.bed", outputDir, cnt.thresh), quote=F,sep="\t",col.names=F,row.names=F)





