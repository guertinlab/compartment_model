library(ggplot2)
library(dplyr)
# contains source colde of plotRates and plotRates_noKpre
#source("/home/rudramukh/Documents/COPASI_compartmentModel/plotRateFuncs.R")
source("plotRateFuncs.R")
#repeatScanParamDir = "/home/rudramukh/Documents/COPASI_compartmentModel/compPolModelpkg/paramfitdir_TBP_dTAGv1/repeatedParamScans" # Documents/COPASI_compartmentModel/compPolModelpkg/paramfitdir

## TRP 12.5 min - with kpre
# for all params between 1e-8 to 1.
combnResultsFile = "InternalControl_TRP12min_repressed/paramEstimationResults_combined.txt"
plotRates(combnResultsFile, "baseline", "TRP_12.5min", "TRP_12.5min (repressed)", "TRP_12min_repressed_allParamsConstrained", 12)


