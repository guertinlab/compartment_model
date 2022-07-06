library(dplyr)
source("customFunctions.R")
#-----
#setwd("/home/FCAM/rmukherjee/Summer2022Rotation/R_analysis/deseq2")
merged.counts.file="/home/FCAM/rmukherjee/Summer2022Rotation/proseqFiles/Macrophages_PRO_gene_counts_pTA.txt"
merged.counts.obj = read.table(merged.counts.file, sep = '\t', header = TRUE) # verify column names and spacing by importing it in a spreadsheet
rownames(merged.counts.obj) = merged.counts.obj[, 1] # set row names to gene names
merged.counts.obj = merged.counts.obj[,seq(2,to=ncol(merged.counts.obj),by=2)]
save(merged.counts.obj, file="merged.gene.counts.pTA.Rdata")

# run DESeq2 for 45mn vs control and 90 min vs control comparisons
# ----------
cond.list = c(rep("45min", 4), rep("90min", 4), rep("control", 4))
cond.levels=c("45min", "90min", "control")
deseq.dds = run.deseq2(merged.counts.obj, cond.list, cond.levels, "control", fileidentifier="allcond.no.lfc")
res.45minvscontrol = results(deseq.dds, contrast=c("sample.conditions", "45min", "control"))
#deseq.dds = loadRData("dds.allcond.no.lfc.gene.Rdata") # func in customFunctions.R
res.90minvscontrol = results(deseq.dds, contrast=c("sample.conditions", "90min", "control"))
save(res.45minvscontrol, file="res_deseq2_45minvscontrol.Rdata")
save(res.90minvscontrol, file="res_deseq2_90minvscontrol.Rdata")
#--------
res.45min.categories = categorize.deseq.df(res.45minvscontrol, fdr = 0.05, log2fold = 0, treat = "apoptotic_45min")
res.45min.categories.df = as.data.frame(res.45min.categories)
saveRDS(res.45min.categories.df, "res.45minvscontrol.categorized.full.rds")
res.45min.active = as.data.frame(res.45min.categories) %>% filter(response=="apoptotic_45min Activated")

res.90min.categories = categorize.deseq.df(res.90minvscontrol, fdr = 0.05, log2fold = 0, treat = "apoptotic_90min")
res.90min.categories.df = as.data.frame(res.90min.categories)
saveRDS(res.90min.categories.df, "res.90minvscontrol.categorized.full.rds")

res.90min.active = as.data.frame(res.90min.categories) %>% filter(response=="apoptotic_90min Activated")
save(res.45min.active, file="active_deseq2_45minvscontrol.Rdata")
save(res.90min.active, file="active_deseq2_90minvscontrol.Rdata")
#--------
#the comparison between 90min vs 45min case
res.90minvs45min = results(deseq.dds, contrast=c("sample.conditions", "90min", "45min"))
res.90minvs45min.categories = categorize.deseq.df(res.90minvs45min, fdr = 0.05, log2fold = 0, treat = "90minvs45min")
res.90minvs45min.category.df = as.data.frame(res.90minvs45min.categories)
saveRDS(res.90minvs45min.category.df, "res.90minvs45min.categorized.full.rds") # can be read using readRDS() into any object
res.90minvs45min.active = as.data.frame(res.90minvs45min.categories) %>% filter(response=="90minvs45min Activated")
save(res.90minvs45min.active, file="active_deseq2_90minvs45min.Rdata")
#--------
#get gene list for each case as txt file
writeLines(rownames(res.45min.active), "activegene_45minvscontrol.txt", sep="\n")
writeLines(rownames(res.90min.active), "activegene_90minvscontrol.txt", sep="\n")
writeLines(rownames(res.90minvs45min.active), "activegene_90minvs45min.txt", sep="\n")

