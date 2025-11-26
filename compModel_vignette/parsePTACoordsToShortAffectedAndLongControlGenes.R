library(dplyr)

# wrong # df.pause.body <- readRDS("pauseWindowsAndBWplots/TRP_FP_vscntrl.pTA.pause.body.rds")
df.pause.body = read.table("primaryTranscriptAnnotation/mastercoords.bed", header=F)
colnames(df.pause.body) = c("chr", "start", "end", "gene", "score", "strand")

# break long genes into long genes with no signal
controlgenes.len.gt.50kb =  subset(df.pause.body, end - start >= 49999)

control.genes.df = c()
genes.for.affect = c()
colnames.of.interest = c("chr","start", "end", "gene", "score", "strand")
for (i in c(1:nrow(controlgenes.len.gt.50kb))) {
	this.row = controlgenes.len.gt.50kb[i,]
	if ( this.row$strand == "+" ) {
		copy.row = this.row
		this.row$start = this.row$start + 49999
		this.row$gene = sprintf("%s_cntrlGene", this.row$gene)
		control.genes.df = rbind(control.genes.df, this.row[,colnames.of.interest])
		copy.row$end = copy.row$start + 24999
		copy.row$gene = sprintf("%s_affected", copy.row$gene)
		genes.for.affect = rbind(genes.for.affect, copy.row[,colnames.of.interest])
	} else {
		copy.row = this.row
		this.row$end = this.row$end - 49999
		this.row$gene = sprintf("%s_cntrlGene", this.row$gene)
		control.genes.df = rbind(control.genes.df, this.row[,colnames.of.interest])
		copy.row$start = copy.row$end - 24999
		copy.row$gene = sprintf("%s_affected", copy.row$gene)
		genes.for.affect = rbind(genes.for.affect, copy.row[,colnames.of.interest])
       }
}

# shorter genes
controlgenes.len.short.50kb =  subset(df.pause.body, end - start < 49999)
for (i in c(1:nrow(controlgenes.len.short.50kb))) {
	this.row = controlgenes.len.short.50kb[i,]
	if ( this.row$strand == "+" ) {
	copy.row = this.row
	copy.row$end = min(copy.row$end, copy.row$start + 24999)
	copy.row$gene = sprintf("%s_affected", copy.row$gene)
	genes.for.affect = rbind(genes.for.affect, copy.row[,colnames.of.interest])
	} else {
	copy.row = this.row
	copy.row$start = max(copy.row$start, copy.row$end - 24999)
	copy.row$gene = sprintf("%s_affected", copy.row$gene)
	genes.for.affect = rbind(genes.for.affect, copy.row[,colnames.of.interest])
	}
}
write.table(as.data.frame(control.genes.df), file="primaryTranscriptAnnotation/InternalControlGenes_TRP_12min.bed", quote=F, sep="\t", row.names=F, col.names=F)
write.table(as.data.frame(genes.for.affect), file="primaryTranscriptAnnotation/GenesWithRegionsAffected_TRP_12min.bed", quote=F, sep="\t", row.names=F, col.names=F)


