#!/usr/bin/env Rscript

# modified from normalization_factor.R because of different file names 
# assumes the naming convention CELL_COND_TIME_ETC_rep<#>_<plus><minus>.bigWig
# assumes that you are running the script in the directory with your bigWigs 

if(!require(devtools)) { install.packages("devtools", repos='http://cran.us.r-project.org') }
if(!require(bigWig)) {devtools::install_github("andrelmartins/bigWig", subdir="bigWig")}

library(bigWig)

coverage = c()
file.name = c()

# commented out from the original code as it seems not to serve any purpose
#for (plus.bigWig in Sys.glob(file.path("*_rep1_plus_PE1.bigWig"))) { # name changed
#	file.prefix = strsplit(plus.bigWig, '_rep1')[[1]][1]
#	#print(file.prefix)
#    minus.bigWig = paste0(strsplit(plus.bigWig, 'plus')[[1]][1], "plus.bigWig") # should actually be minus.bigWig
#    #print(minus.bigWig)
#    }
# these bigWigs were creataed using seqOutBias
for (replicate.plus.bigWig in Sys.glob(file.path("*_plus.bigWig"))) { # name changed
   	#print(replicate.plus.bigWig)
   	replicate.prefix = strsplit(replicate.plus.bigWig, '_plus.bigWig')[[1]][1] # name changed
	replicate.minus.bigWig = paste0(strsplit(replicate.prefix, '_plus.bigWig')[[1]][1], "_minus.bigWig") # name changed
	#print(replicate.minus.bigWig)
	bigwig.plus = load.bigWig(replicate.plus.bigWig)
	bigwig.minus = load.bigWig(replicate.minus.bigWig)
	coverage.bigwig.plus = bigwig.plus$basesCovered*bigwig.plus$mean # coveredBases * mean (mysterious because I don't know the definitions)
	coverage.bigwig.minus = abs(bigwig.minus$basesCovered*bigwig.minus$mean) # abs used because of negstive values
	coverage = append(coverage, coverage.bigwig.plus)
	coverage = append(coverage, coverage.bigwig.minus)
	file.name = append(file.name, replicate.plus.bigWig)
	file.name = append(file.name, replicate.minus.bigWig)
	#print(coverage)
	unload.bigWig(bigwig.plus)
	unload.bigWig(bigwig.minus)
}	
# NOTE: CHANGE THE FOllowing prefix to `basename` in PRO_normalization (line 34)
#celltype.prefix = strsplit(file.prefix, '_')[[1]][1]
celltype.prefix = "V65"	
write.table(data.frame(file.name, coverage), file = paste0(celltype.prefix, '_normalization.txt'), sep = '\t', row.names = FALSE, col.names=FALSE, quote =FALSE)

