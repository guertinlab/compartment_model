# calls functions in pause_workflow_v2.R to plot composites from conditions.
# this code was written for a project where there were three conditions - control, 45min and 90min.
library(dplyr)
library(bigWig)

# source the functions
source("pause_workflow_funcs_v2.R")

## NOTE: CHANGE these PATHS before running on your own platform.
# set the working directory
outputDir = '/home/FCAM/rmukherjee/compartmentModel_project_datasets/Jonkers_et_al_2014/pauseWindowsAndBWplots'
# load the merged bigWigs from all conditions x replicates
merged.bw.dir = "/home/FCAM/rmukherjee/compartmentModel_project_datasets/Jonkers_et_al_2014/fastqfiles"
#bw.plus <- load.bigWig(sprintf('%s/TRP_FP_control_merged_plus.bigWig', merged.bw.dir))
#bw.minus <- load.bigWig(sprintf('%s/TRP_FP_control_merged_minus.bigWig', merged.bw.dir))

normalized.bw.dir = "/home/FCAM/rmukherjee/compartmentModel_project_datasets/Jonkers_et_al_2014/fastqfiles"

bw.cntrl.plus = load.bigWig(sprintf('%s/V65_cntrl_plus_normalized_sizeFactor_TRP.bigWig', normalized.bw.dir))
bw.cntrl.minus = load.bigWig(sprintf('%s/V65_cntrl_minus_normalized_sizeFactor_TRP.bigWig', normalized.bw.dir))

bw.TRP.plus = load.bigWig(sprintf('%s/V65_TRP_treat_12min_plus_normalized_sizeFactor_TRP.bigWig', normalized.bw.dir))
bw.TRP.minus = load.bigWig(sprintf('%s/V65_TRP_treat_12min_minus_normalized_sizeFactor_TRP.bigWig', normalized.bw.dir))

deseq2.files.dir = "/home/FCAM/rmukherjee/compartmentModel_project_datasets/Jonkers_et_al_2014/deseq2"
categorized.TRPtreatvscontrol = readRDS(sprintf("%s/deseqGenesInfo.TRP12min.matched.rds",deseq2.files.dir))
#categorized.TRPtreatvscontrol$gene = categorized.TRPtreatvscontrol$gene %>% str_replace("_affected", "")
categorized.TRPtreatvscontrol  = categorized.TRPtreatvscontrol[!duplicated(categorized.TRPtreatvscontrol$gene),]
rownames(categorized.TRPtreatvscontrol) = categorized.TRPtreatvscontrol$gene

# find pause region for all genes listed in primaryTranscriptAnnotation bed file using the merged bigWigs 
# load the coordinates file
affectedGenes.newTTS.bed = "/home/FCAM/rmukherjee/compartmentModel_project_datasets/Jonkers_et_al_2014/primaryTranscriptAnnotation/GenesWithRegionsAffected_TRP_12min_sorted_renamed.bed"
pTA <- read.table(affectedGenes.newTTS.bed, sep="\t")
df.pause.body <- find.pause.regions(pTA,bw.cntrl.plus,bw.cntrl.minus, searchWindow = 300)
rownames(df.pause.body) = df.pause.body$gene
saveRDS(df.pause.body, file=sprintf("%s/Internal_TRP_vscntrl.pause.body.rds",outputDir))

merged.pbody.TRPtreat.vs.control = merge(df.pause.body, as.data.frame(categorized.TRPtreatvscontrol))
merged.pbody.TRPtreat.vs.control = merged.pbody.TRPtreat.vs.control[,c(colnames(df.pause.body), colnames(categorized.TRPtreatvscontrol))]
rownames(merged.pbody.TRPtreat.vs.control) = merged.pbody.TRPtreat.vs.control$gene
saveRDS(merged.pbody.TRPtreat.vs.control, file=sprintf("%s/Internal.merged.pbody.TRP.vs.control.rds",outputDir))

## find genes with high variability
result_df = find.density.across.windows(df.pause.body, bw.cntrl.plus, bw.cntrl.minus, start_bp=0, end_bp=15000,windowsize_bp=200, stepsize_bp=50)
saveRDS(result_df, sprintf("%s/TRP.body.density.rollingwindow.windowsize=200,stepsize=50,bp_end=15000.rds", outputDir))
summary.TRP.windows = group_by(result_df, gene) %>% summarise_at(c("body.avg"), c("mean", "median", "var", "min", "max"))
summary.TRP.windows$var_by_mean = summary.TRP.windows$var / summary.TRP.windows$mean
saveRDS(summary.TRP.windows, sprintf("%s/summary.TRP.windows.rds", outputDir))
#summary.TRP.windows = readRDS( sprintf("%s/summary.TRP.windows.rds", outputDir))
filteredGenes = subset(summary.TRP.windows, var_by_mean <= 0.1)$gene

setwd(outputDir)
# call plotting steps for all genes
run.plotting.steps(merged.pbody.TRPtreat.vs.control, "control_all", bw.cntrl.plus, bw.cntrl.minus,
		        "TRP_12.5min_all", bw.TRP.plus, bw.TRP.minus,
			"TRP_12.5min (all genes, Internal Control)", color.names=c(rgb(1,0,0,1/2), rgb(0,0,1,1/2))) # red, blue from rgb value. Indices set in alphabetical order of column names
density.TRP = readRDS("density.TRP_12.5min_all.v.control_all.rds")
repressed.TRP.baseline = subset(density.TRP, cond == "control_all" & treatment == "Repressed")
repressed.TRP.treatment = subset(density.TRP, cond == "TRP_12.5min_all" & treatment == "Repressed")

colNames = c("gene", "pause.sum", "body.avg")
write.table(repressed.TRP.baseline[,colNames], "InternalControl_repressed_TRP125min_baseline.txt", sep="\t", quote=F, row.names=F)
write.table(repressed.TRP.treatment[,colNames], "InternalControl_repressed_TRP125min_treatment.txt", sep="\t", quote=F, row.names=F)

# call the run.plotting.steps on pause body object containing genes of our interest.
repressed.TRPvscontrol = subset(merged.pbody.TRPtreat.vs.control, treatment == "Repressed")
saveRDS(repressed.TRPvscontrol, file=sprintf("%s/InternalControl.repressed.pbody.TRP.vs.control.rds", outputDir))
run.plotting.steps(repressed.TRPvscontrol, "control", bw.cntrl.plus, bw.cntrl.minus,
		        "TRP_12.5min (repressed, internal spikein)", bw.TRP.plus, bw.TRP.minus,
			"TRP_12.5min(rrepressed, internal control)", color.names=c(rgb(1,0,0,1/2), rgb(0,0,1,1/2))) # red, blue from rgb value. Indices set in alphabetical order of column names

## USE GENES with body variance / mean <= 0.1
repressed.TRP.baseline = subset(density.TRP, cond == "control_all" & treatment == "Repressed")
repressed.TRP.treatment = subset(density.TRP, cond == "TRP_12.5min_all" & treatment == "Repressed")
control.TRP.allgenes.varlessthanpoint1 = subset(density.TRP, cond == "control_all" & gene %in% filteredGenes)
repressed.TRP.baseline.varlessthanpoint1 = subset(density.TRP, cond == "control_all" & treatment == "Repressed" & gene %in% filteredGenes)
repressed.TRP.treatment.varlessthanpoint1 = subset(density.TRP, cond == "TRP_12.5min_all" & treatment == "Repressed"& gene %in% filteredGenes)

colNames = c("gene", "pause.sum", "body.avg")
write.table(subset(control.TRP.allgenes.varlessthanpoint1,pause.sum>0)[, colNames], "TRP125min_control_allgenes_varlessthanpoint1.txt", sep="\t", quote=F, row.names=F)
write.table(repressed.TRP.baseline.varlessthanpoint1[, colNames], "TRP125min_repressed_control_varlessthanpoint1.txt", sep="\t", quote=F, row.names=F)
write.table(repressed.TRP.treatment.varlessthanpoint1[, colNames], "TRP125min_repressed_treatment_varlessthanpoint1.txt", sep="\t", quote=F, row.names=F)

