
### For TRP 12.5min treatment
#genecountsDir = "/Users/rudradeepmukherjee/Documents/UConn Health/comparment.model.project/datasets, genome/Jonkers_2014_Trp_FP_paper"
genecountsDir = "/home/FCAM/rmukherjee/compartmentModel_project_datasets/Jonkers_et_al_2014"
outputDir = sprintf("%s/deseq2", genecountsDir)
size.factors = readRDS(sprintf("%s/size.factors.rawCounts.TRPInternalGenes.rds", outputDir))
#size.factors = readRDS(sprintf("%s/size.factors.modifiedFormula.rawCounts.TRPInternalGenes.rds", outputDir))
size.factors.transpose = t(t(size.factors))
write.table(size.factors.transpose, sprintf("%s/sizefactors_TRP_12min_InternalGenes.txt",outputDir), sep="\t", quote=F, col.names=F)

### For FP 25min treatment
#genecountsDir = "/Users/rudradeepmukherjee/Documents/UConn Health/comparment.model.project/datasets, genome/Jonkers_2014_Trp_FP_paper"
genecountsDir = "/home/FCAM/rmukherjee/compartmentModel_project_datasets/Jonkers_et_al_2014"
outputDir = sprintf("%s/deseq2", genecountsDir)
#size.factors = readRDS(sprintf("%s/size.factors.FPInternalGenes.rds", outputDir))
size.factors = readRDS(sprintf("%s/size.factors.modifiedFormula.rawCounts.FPInternalGenes.rds", outputDir))
size.factors.transpose = t(t(size.factors))
write.table(size.factors.transpose, sprintf("%s/sizefactors_FP_25min_InternalGenes.txt",outputDir), sep="\t", quote=F, col.names=F)



