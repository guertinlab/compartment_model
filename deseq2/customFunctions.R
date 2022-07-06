# collection of functions to ease flow and reduce redundancy of code
library(fgsea)
library(apeglm)
library(DESeq2)
library(lattice)
library(dplyr)
source('/home/FCAM/rmukherjee/Spring2022Rotation/R_files/ZNF143_functions.R')

# using fgsea and allowing the use of lfc.shrink 
run.gesea.analysis <- function(dds, baseFilename="test", use.lfcShrink = T, gencode.key, fdr = 0.05, treat = "treatment", log2fold = 0, ma.plot.title, gene.ontology.file = "msigdb.v7.0.symbols.gmt", gsea.analysis.name = "test_vs_control") {
   #dds = DESeq(deseq.df) # dds is already passed to the function
   #save(dds, file = sprintf("dds.%s.Rdata", baseFilename))
   ddsObj <- dds
   if (use.lfcShrink) {
    ddsObj = lfcShrink(dds = dds, coef = 2)
    baseFilename = sprintf("%s.lfc", baseFilename)
    filename = sprintf("dds.%s.Rdata", baseFilename)
    print(sprintf("saving %s", filename))
   save(ddsObj, file = filename)   
  ddsObj.lattice = categorize.deseq.df(ddsObj, fdr = fdr, log2fold = log2fold, treat = treat)
  } else {
  ddsObj.lattice = categorize.deseq.df(results(ddsObj), fdr = fdr, log2fold = log2fold, treat = treat)
}
#ma.plot.lattice(ddsObj.lattice, filename = sprintf("MM_macrophage_HumanT_%s", baseFilename), title.main = ma.plot.title)
merged.counts.pre = merge(as.data.frame(ddsObj.lattice), gencode.key, by = "row.names", all.x = F)
merged.counts = merged.counts.pre[, c(2:(ncol(merged.counts.pre) - 2))]

rownames(merged.counts) = make.names(merged.counts.pre$gene, unique = TRUE)

 filename =  sprintf("merged.counts.%s.Rdata", baseFilename) 
  print(sprintf("saving %s", filename))
  save(merged.counts, file = filename)
  # ordering the genes
  x = merged.counts$log2FoldChange
  x[x < 0] <- -1
  x[x >= 0] <- 1
  gene.list = x * (2^abs(merged.counts$log2FoldChange)) * -log(merged.counts$pvalue,
    base = 10) 
  names(gene.list) = row.names(merged.counts)
  # rank order
  gene.list = sort(gene.list, decreasing = TRUE)
  # remove NA entries
  gene.list = gene.list[!is.na(gene.list)]
  # remove duplicated genes
  gene.list = gene.list[!duplicated(names(gene.list))]
  filename = sprintf("gene.list.%s.Rdata", baseFilename)
  print(sprintf("saving %s", filename))
  save(gene.list, file = filename)
   
  ## gene enrichment analysis
  myGO = gmtPathways(gene.ontology.file)
  fgRes = as.data.frame(fgsea(pathways = myGO, stats = gene.list,
    minSize = 10, maxSize = 600))
  filename = sprintf("fgsea.list.%s.Rdata", baseFilename)
  save(fgRes, file=filename)
  print(sprintf("saved %s", filename))
  fgRes[fgRes$padj < 0.1, ][order(fgRes$NES[fgRes$padj < 0.1],
    decreasing = TRUE), ]
  filename = sprintf("fgsea.list.%s.signif.Rdata", baseFilename)
  save(fgRes, file=filename)
  print(sprintf("saved %s", filename))

   filename =paste("GSEA_", baseFilename, "_", gsea.analysis.name, ".pdf", sep = "") 
   print(sprintf("saving %s", filename))
 #  return(fgRes) # return for debugging
   
   pdf(filename, useDingbats = FALSE,
    width = 6.83, height = 3.83)
    print(plotEnrichment(myGO[[fgRes[fgRes$padj < 0.1, ][order(fgRes$NES[fgRes$padj <
    0.1], decreasing = TRUE), ][[1]][1]]], gene.list) + labs(title = gsub("_", " ", fgRes[fgRes$padj <
    0.1, ][order(fgRes$NES[fgRes$padj < 0.1], decreasing = TRUE),][[1]][1])))
dev.off()
   return(fgRes)
}

# implement differential_expression_v2.R as a function to prevent repeated code
# merged.counts.file contains gene counts across all conditions. See section for DESeq2 in Guertin Lab's Nascent_RNA methods
# condition.list is the name of factors corresponding to columns in the merged.counts file. Say, 45min, 45min, 90min, 90min
# condition.levels is the number of levels in the condition.list
# reference.level is the control level which is explicitly set by relevel
# fileidentifier is the substring to be used while saving R object
# treatment is a descriptive name of the treatment used by categorize.deseq.df() in ZNF143_functions.R to name active and repressed genes
# fdr is the value passed to categorize.deseq.df() function
# CAUTION: two levels (control and treatment) currently supported as lfcShrink object doesn't allow an easy way to toggle between more than two conditions
# merged.counts.obj has been prepared in the following way
# meged.counts.obj = read.table(merged.counts.file, sep = '\t', header = TRUE) # built in DESeq2 section from Guertin lab's Nascent RNA methods
# rownames(merged.counts.obj) = merged.counts.obj[, 1]
# merged.counts.obj = merged.counts.obj[,seq(2,to=ncol(merged.counts.obj),by=2)] # append the columns with gene counts (side by side) in the final table
run.deseq2.lfcShrink = function(merged.counts.obj, condition.list, condition.levels, reference.level, fileidentifier="test", treatment="treatment", fdr=fdr) {
  sample.conditions = factor(condition.list, levels=condition.levels)
  sample.conditions = relevel(sample.conditions, ref=reference.level)
  deseq.df = DESeqDataSetFromMatrix(merged.counts.obj, DataFrame(sample.conditions), ~ sample.conditions)
  dds = DESeq(deseq.df)
  dds.filename = sprintf("dds.%s.gene.Rdata", fileidentifier)
  print(sprintf("saving dds object as %s", dds.filename))
  save(dds, file= dds.filename)
  #----------
  dds.lfcShrink = lfcShrink(dds = dds, coef = 2)
  dds.lfc.filename = sprintf("dds.%s.lfcshrink.Rdata", fileidentifier)
  print(sprintf("saving dds.lfcShrink object as %s", dds.lfc.filename))
  save(dds.lfcShrink, file=dds.lfc.filename)
  #---------
  dds.lfc.categorized = categorize.deseq.df(dds.lfcShrink, fdr = fdr, log2fold = 0, treat = treatment)
  dds.categorized.filename = sprintf("dds.%s.categorized.Rdata", fileidentifier)
  print(sprintf("saving dds.categorized object as %s", dds.categorized.filename))
  save(dds.lfc.categorized, file=dds.categorized.filename)
  #---------
  return(dds.lfc.categorized)  
}

# without lfcShrink
run.deseq2 = function(merged.counts.obj, condition.list, condition.levels, reference.level, fileidentifier="test") {
  sample.conditions = factor(condition.list, levels=condition.levels)
  sample.conditions = relevel(sample.conditions, ref=reference.level)
  deseq.df = DESeqDataSetFromMatrix(merged.counts.obj, DataFrame(sample.conditions), ~ sample.conditions)
  dds = DESeq(deseq.df)
  dds.filename = sprintf("dds.%s.gene.Rdata", fileidentifier)
  print(sprintf("saving dds object as %s", dds.filename))
  save(dds, file= dds.filename)
  return(dds)
  # categorize.deseq.df can be used from ZNF143_functions.R to add a response column with categories with proper factor level  
}

# load Rdata and send it to the caller, who can assign it to a variable name of personal choice
# source: https://stackoverflow.com/questions/5577221/how-can-i-load-an-object-into-a-variable-name-that-i-specify-from-an-r-data-file
loadRData <- function(fileName){
#loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}

