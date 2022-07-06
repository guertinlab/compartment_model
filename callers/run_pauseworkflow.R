# calls functions in pause_workflow_v2.R to plot composites from conditions.
# this code was written for a project where there were three conditions - control, 45min and 90min.
library(dplyr)
library(bigWig)

# source the functions
source("../R/pause_workflow_funcs_v2.R")

# set the working directory
#setwd('/home/FCAM/rmukherjee/Summer2022Rotation/R_analysis/pauseWindowsAndBWplots')

# load the merged bigWigs from all conditions x replicates
#merged.bw.dir = "/home/FCAM/rmukherjee/Spring2022Rotation/proseqfiles/merged_sorted_bigWig_files/"
merged.bw.dir = "../data"
bw.plus <- load.bigWig(sprintf('%s/Macrophages_merged_sorted_plus_PE1.bigWig', merged.bw.dir))
bw.minus <- load.bigWig(sprintf('%s/Macrophages_merged_sorted_minus_PE1.bigWig', merged.bw.dir))

# get normalized bigwigs from different conditions
# generated from PRO_normalization in "modeling_PRO_composites" section on Github
#normalized.bw.dir = "/home/FCAM/rmukherjee/Spring2022Rotation/proseqfiles/"
normalized.bw.dir = "../data"
bw.cntrl.plus = load.bigWig(sprintf('%s/Macrophages_plus_PE1_normalized.bigWig', normalized.bw.dir))
bw.cntrl.minus = load.bigWig(sprintf('%s/Macrophages_minus_PE1_normalized.bigWig', normalized.bw.dir))

bw.cond1.plus = load.bigWig(sprintf('%s/Macrophages_apoptotic_45min_plus_PE1_normalized.bigWig', normalized.bw.dir))
bw.cond1.minus = load.bigWig(sprintf('%s/Macrophages_apoptotic_45min_minus_PE1_normalized.bigWig', normalized.bw.dir))

bw.cond2.plus =  load.bigWig(sprintf('%s/Macrophages_apoptotic_90min_plus_PE1_normalized.bigWig', normalized.bw.dir))
bw.cond2.minus = load.bigWig(sprintf('%s/Macrophages_apoptotic_90min_minus_PE1_normalized.bigWig', normalized.bw.dir))

# load deseq2 result objects to be used in subsequent analysis
# these objects were created in this flow:-- results obj (DESeq2) -> categorize.deseq.df() [ZNF143_functions.R] -> as.data.frame() 

#deseq2.files.dir = "/home/FCAM/rmukherjee/Summer2022Rotation/R_analysis/deseq2" # the files in this directory was created by callDESeq2.R
deseq2.files.dir = "../data"
# we have three comparisons - 45min vs control, 90min vs control and 90min vs 45min
categorized.45minvscontrol = readRDS(sprintf("%s/res.45minvscontrol.categorized.full.rds", deseq2.files.dir))

categorized.90minvscontrol = readRDS(sprintf("%s/res.90minvscontrol.categorized.full.rds", deseq2.files.dir))

categorized.90minvs45min = readRDS(sprintf("%s/res.90minvs45min.categorized.full.rds", deseq2.files.dir))

# find pause region for all genes listed in primaryTranscriptAnnotation bed file using the merged bigWigs 
# load the coordinates file
#master.pTA.coords.file= '/home/FCAM/rmukherjee/Summer2022Rotation/R_analysis/primTrnscript/mastercoords.bed'
master.pTA.coords.file = '../data/mastercoords.bed'
pTA <- read.table(master.pTA.coords.file)
df.pause.body <- find.pause.regions(pTA,bw.plus,bw.minus)
rownames(df.pause.body) = df.pause.body$gene
saveRDS(df.pause.body, file="../data/Macrophages.merged.pause.body.rds")

# merge the pause body object with deseq2 objects based on gene name
# ASSUMES that deseq2.df and df.pause.body have rownames set to gene names
# merge.pbody.deseq2 returns an object where first set of columns comes from df.pause.body and second set of columns comes from the deseq2 object
merged.pbody.45minvscontrol = merge.pbody.deseq2(df.pause.body, categorized.45minvscontrol)
saveRDS(merged.pbody.45minvscontrol, file="../data/merged.pausebody.45minvscontrol.rds")

merged.pbody.90minvscontrol = merge.pbody.deseq2(df.pause.body, categorized.90minvscontrol)
saveRDS(merged.pbody.90minvscontrol, file="../data/merged.pausebody.90minvscontrol.rds")

merged.pbody.90minvs45min = merge.pbody.deseq2(df.pause.body, categorized.90minvs45min)
saveRDS(merged.pbody.90minvs45min, file="../data/merged.pausebody.90minvs45min.rds")

# call the run.plotting.steps on pause body object containing genes of our interest.
activated.45minvscontrol = filter(merged.pbody.45minvscontrol, response == "apoptotic_45min Activated")
saveRDS(activated.45minvscontrol, file="../data/activated.45minvscontrol.rds")
run.plotting.steps(activated.45minvscontrol, "control", bw.cntrl.plus, bw.cntrl.minus,
                    "45min (activated)", bw.cond1.plus, bw.cond1.minus, "efferocytosis",
                   color.names=c(rgb(1,0,0,1/2), rgb(0,0,1,1/2))) # red, blue from rgb value. Indices set in alphabetical order of column names

activated.90minvscontrol = filter(merged.pbody.90minvscontrol, response == "apoptotic_90min Activated")
saveRDS(activated.90minvscontrol, file="../data/activated.90minvscontrol.rds")
run.plotting.steps(activated.90minvscontrol, "control", bw.cntrl.plus, bw.cntrl.minus,
                    "90min (activated)", bw.cond2.plus, bw.cond2.minus, "efferocytosis",
                                        color.names=c(rgb(1,0,0,1/2), rgb(0,0,1,1/2))) # R and B color. 

activated.90minvs45min = filter(merged.pbody.90minvs45min, response == "90minvs45min Activated")
saveRDS(activated.90minvs45min, file="../data/activated.90minvs45min.rds")
run.plotting.steps(activated.90minvs45min, "45min", bw.cond1.plus, bw.cond1.minus,
                    "90min (activated)", bw.cond2.plus, bw.cond2.minus, "efferocytosis",
                                        color.names=c(rgb(0,0,1,1/2), rgb(1,0,0,1/2))) # B and R color.

