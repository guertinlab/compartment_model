library(MatchIt)
library(bigWig)
library(DESeq2)

source("/home/FCAM/rmukherjee/Spring2022Rotation/R_files/ZNF143_functions.R")
setwd("/home/FCAM/rmukherjee/compartmentModel_project_datasets/Jonkers_et_al_2014/deseq2")

outputDir = "/home/FCAM/rmukherjee/compartmentModel_project_datasets/Jonkers_et_al_2014/deseq2"
pTA_dir = "/home/FCAM/rmukherjee/compartmentModel_project_datasets/Jonkers_et_al_2014/primaryTranscriptAnnotation"
fastqfilesdir = "/home/FCAM/rmukherjee/compartmentModel_project_datasets/Jonkers_et_al_2014/fastqfiles"
genecountsDir = "/home/FCAM/rmukherjee/compartmentModel_project_datasets/Jonkers_et_al_2014/deseq2"
spikeInCountsFile = sprintf("%s/PRO_TRP_V65_InternalControl_GeneCounts.txt", genecountsDir)
mouseGeneCountFile = sprintf("%s/PRO_TRP_V65_Affected_GeneCounts.txt", genecountsDir) # made from primary transcript annotation
counts.obj.internal = read.table(spikeInCountsFile, sep = '\t', header = T)
row.names.internal = counts.obj.internal[,1]
counts.obj.internal = counts.obj.internal[,seq(2,to=ncol(counts.obj.internal),by=2)]
rownames(counts.obj.internal) = row.names.internal

## only consider TRP and control files - regex V65_[cT]
## only consider genes with length upto 25kb 
inferred.coords=read.table(sprintf('%s/GenesWithRegionsAffected_TRP_12min_sorted.bed', pTA_dir), sep="\t", header =FALSE)
inferred.coords[,"V4"] = gsub("_affected", "", inferred.coords[,"V4"])
counts.df = abs(get.raw.counts.interval(inferred.coords, fastqfilesdir, file.prefix = "V65_[cT]", 
		file.plus.suffix = "*rep?_plus.bigWig", file.minus.suffix = "_minus.bigWig"))
saveRDS(counts.df, "raw.counts.TRP12min_affected.rds")
merged.counts.obj = rbind(counts.obj.internal, counts.df)
saveRDS(merged.counts.obj, file = "merged.internal.plus.raw.counts.rds")
cond.list = c("control", "control", "TRP_12min", "TRP_12min")
cond.levels=c("TRP_12min", "control")
sample.conditions = factor(cond.list, levels = cond.levels)
sample.conditions = relevel(sample.conditions, ref="control")

# #PCA for experiments
dds = run.deseq.list.dds(merged.counts.obj, sample.conditions, cond.levels, ControlGenes.Indices=1:nrow(counts.obj.internal))
rld <- rlog(dds)
plotPCA(rld) # basic plot

pca.plot = plotPCA(rld, intgroup="sample.conditions", returnData=TRUE)
pca.plot$sample.conditions = rownames(pca.plot)
plotPCAlattice(pca.plot, file = 'PCA_TRP_12min_PRO.pdf')

rep = factor(sapply(strsplit(colnames(counts.df), 'rep'), '[', 2))
# formaula: ~ rep + sample.conditions
deseq.df = DESeqDataSetFromMatrix(counts.df,colData = cbind.data.frame(sample.conditions, rep), ~ rep + sample.conditions)
deseq.df$sizeFactor = estimateSizeFactorsForMatrix(merged.counts.obj, controlGenes = 1:nrow(counts.obj.internal))
saveRDS(deseq.df$sizeFactor, file = sprintf("%s/size.factors.rawCounts.TRPInternalGenes.rds", outputDir))
deseq.df = DESeq(deseq.df)
saveRDS(deseq.df, "deseq.df.TRPinternal.RAW.counts.rds")

# formula: ~ sample.conditions
deseq.df = DESeqDataSetFromMatrix(counts.df,colData = DataFrame(sample.conditions), ~ sample.conditions)
deseq.df$sizeFactor = estimateSizeFactorsForMatrix(merged.counts.obj, controlGenes = 1:nrow(counts.obj.internal))
saveRDS(deseq.df$sizeFactor, file = sprintf("%s/size.factors.modifiedFormula.rawCounts.TRPInternalGenes.rds", outputDir))
deseq.df = DESeq(deseq.df)
saveRDS(deseq.df, "deseq.df.modifiedFormula.TRPinternal.RAW.counts.rds")

# Create tables of activated/repressed genes and their matching genes 
for (i in levels(sample.conditions)[-1]) {
	res.deseq = results(deseq.df, contrast = c("sample.conditions", i, "control"))
  sum(res.deseq$padj < 0.1 & !is.na(res.deseq$padj))
    activated = res.deseq[res.deseq$padj < 0.1 & !is.na(res.deseq$padj) & res.deseq$log2FoldChange > 0,]
    activated.strand = merge(cbind(tighten_summit_window(activated), "Activated"), inferred.coords, by.x = "gene", by.y = "V4")[,c(2, 3, 4, 1, 5, 10)]    
      write.table(activated.strand, file = paste0(i, '_activated_genes.bed'), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
      unchanged = res.deseq[!is.na(res.deseq$padj) & res.deseq$padj > 0.1 & abs(res.deseq$log2FoldChange) < 0.01,]
        #unchanged = unchanged[sample(x, size, replace = FALSE, prob = NULL),]
        unchanged$treatment = 0
        activated$treatment = 1
	  df.deseq.effects.lattice = rbind(unchanged, activated)
	  out = matchit(treatment ~ baseMean, data = df.deseq.effects.lattice, method = "optimal", ratio = 1)
	  summary(out, un=F) # don't show stats from pre-matching data
	  pdf("matchIt.optimal.activated.vs.unchanged.TRP12min.pdf")
	  plot(out, type="density", interactive = FALSE)
	  dev.off()  
	  
	  # chose unchanged genes which are matched to activated
	  unchanged = df.deseq.effects.lattice[rownames(df.deseq.effects.lattice) %in% out$match.matrix,]
	  unchanged.strand = merge(cbind(tighten_summit_window(unchanged), "Matched to Activated"), inferred.coords, by.x = "gene", by.y = "V4")[,c(2, 3, 4, 1, 5, 10)] 
	      write.table(unchanged.strand, file = paste0(i, '_activated_matched_unchanged_genes.bed'), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
	      unchanged$treatment = "Matched to Activated"
        activated$treatment = "Activated"
        df.x = rbind(activated, unchanged)
	# do the similar procedure using repressed genes
	  repressed = res.deseq[res.deseq$padj < 0.1 & !is.na(res.deseq$padj) & res.deseq$log2FoldChange < 0,]
	  repressed.strand = merge(cbind(tighten_summit_window(repressed), "Repressed"), inferred.coords, by.x = "gene", by.y = "V4")[,c(2, 3, 4, 1, 5, 10)]    
	    write.table(repressed.strand, file = paste0(i, '_repressed_genes.bed'), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
	    unchanged = res.deseq[!is.na(res.deseq$padj) & res.deseq$padj > 0.1 & abs(res.deseq$log2FoldChange) < 0.01,]

      # unchanged = res.deseq[rownames(res.deseq) %notin% rownames(repressed) & !is.na(res.deseq$padj) & res.deseq$log2FoldChange > 0,]
       unchanged$treatment = 0
      repressed$treatment = 1
        df.deseq.effects.lattice = rbind(unchanged, repressed)
        out = matchit(treatment ~ baseMean, data = df.deseq.effects.lattice, method = "optimal", ratio = 1)
	summary(out, un=F)
        ## Warning message: Fewer control units than treated units; not all treated units will get a match.
	pdf("matchIt.optimal.repressed.vs.unchanged.TRP12min.pdf")
	  plot(out, type="density", interactive = FALSE)
	  dev.off()
	  unchanged = res.deseq[!is.na(res.deseq$padj) & res.deseq$padj > 0.1 & abs(res.deseq$log2FoldChange) < 0.01,]  
	  unchanged = df.deseq.effects.lattice[rownames(df.deseq.effects.lattice) %in% out$match.matrix,]
	  unchanged.strand = merge(cbind(tighten_summit_window(unchanged), "Matched to Repressed"), inferred.coords, by.x = "gene", by.y = "V4")[,c(2, 3, 4, 1, 5, 10)] 

	    write.table(unchanged.strand, file = paste0(i, '_repressed_matched_unchanged_genes.bed'), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
	    unchanged$treatment = "Matched to Repressed"
	      repressed$treatment = "Repressed"
	      df.x = rbind(df.x, unchanged)
	        df.x = rbind(df.x, repressed)
}

## MA Plot for activated, repressed, unchanged

temp = df.x
unique(temp$treatment)

dim(res.deseq[is.na(res.deseq$padj),])
otherGenes = res.deseq[!(rownames(res.deseq) %in% rownames(temp)) & !(is.na(res.deseq$padj)),]
otherGenes$treatment = "Other"

ma.df = rbind(temp,otherGenes)
dim(ma.df)
# `matched to activated` and `matched to repressed` are labels of unchanged genes
levelsMA = c("Other","Matched to Activated", "Matched to Repressed","Activated","Repressed")
ma.df$treatment = factor(ma.df$treatment, levels = levelsMA)
saveRDS(ma.df, file = "TRP_treat_matched_DESeq_genes_categorized.rds")
ma.plot.lattice(ma.df, filename='TRP_treat_matched',title.main = "Diff Exp (TRP 12.5min)", col=c("grey90","grey30","grey30", "#ce228e" ,"#2290cf"))
deseqGenes = tighten_summit_window(ma.df)
dim(deseqGenes)

deseqGenesInfo = cbind(deseqGenes, ma.df)
dim(deseqGenesInfo)
# save cateogrized data
saveRDS(deseqGenesInfo, "deseqGenesInfo.TRP12min.matched.rds")
write.table(deseqGenesInfo,file="TRP_treat_matched_DESeq_genes_categorized.bed",quote = FALSE, sep ='\t', row.names = FALSE, col.names = TRUE)

