library(ggplot2)
library(scales)

# genefile = the file to read gene names from
# repeatScanParamDir = the directory where outputFiles of the form `repeatedParamEstimates_<gene>_onlyvalidfits.txt` exists
# exprName = name of experiment to be used for referencing column names.
# filePrefix = an identifier, used while saving plots and files
# colNames = columns to parsse from the repeatedParamScans
# uniqGeneNum = how many unique genes to use to create plots. Default is 10.
plotRates <- function(combinedResultsFile, baselineName, exprName, yLabel, fileidentifier=NULL, uniqGeneNum = 10, ALPHA=0.03) {
 #P_base = "Values.P_base."  #'Values[P_base]'
 #P_expr = "Values.P_expr." # 'Values[P_expr]'
 fontTheme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 16),
  axis.text.y = element_text(size = 16),
  axis.title.y = element_text(size = 18),
  legend.title = element_text(size=14),
  legend.text = element_text(size=13),
 legend.position="top")

  if (is.null(fileidentifier)) {
	 fileidentifier = sprintf("%s_vs_%s", exprName, baselineName)
	 print(sprintf("setting filePrefix to %s", fileidentifier))
 }
 kelong = "elongationRate" # 'Values[kelong]'
 kinit_base = sprintf("initiationRate_%s", baselineName) #'Values[kinit_base]'
 kinit_expr = sprintf("initiationRate_%s", exprName) #'Values[kinit_expr]'
 kpre = "prematureTerminationRate" # 'Values[kpre]'
 krel_base = sprintf("pauseReleaseRate_%s", baselineName) #'Values[krel_base]'
 krel_expr = sprintf("pauseReleaseRate_%s", exprName) #'Values[krel_expr]'
 colNamesToParse = c("gene", kelong, kpre, kinit_base, kinit_expr, krel_base, krel_expr)
 combinedResults = read.table(combinedResultsFile,sep="\t", header=T)
 # if (file.exists(genefile)) { 
#    geneInfoTable = read.table(genefile, sep="\t", header=T)
# } else {
#        stop(sprintf("genefile %s not found", genefile))
# }
 # start reading in the repeated param estimation results
 #rateData = c() 
 rateData = combinedResults[,colNamesToParse]
 rateData = as.data.frame(rateData)
 colnames(rateData) = c("gene", "kelong", "kpre", "kinit_base", "kinit_expr", "krel_base", "krel_expr")
 rateData$kinit_FC = rateData$kinit_expr / rateData$kinit_base
 rateData$krel_FC = rateData$krel_expr / rateData$krel_base
 saveRDS(rateData, sprintf("repeatedParamScanData_%s.rds", fileidentifier))
 dodge = position_dodge(width=0.05)
 # min, max kinit_FC, krel_FC
 minMaxSummary =  group_by(rateData, gene) %>% summarise_at(c("kinit_FC", "krel_FC"), c("min", "max"), na.rm=T)
 lenGene = nrow(minMaxSummary)
 stepSize = 0.25
 minMaxSummary$Ycoord = seq(0.3, 0.3 + stepSize * (lenGene - 1), stepSize)
 minMaxSummary$kinitFCRange = minMaxSummary$kinit_FC_max - minMaxSummary$kinit_FC_min
 minMaxSummary$krelFCRange = minMaxSummary$krel_FC_max - minMaxSummary$krel_FC_min
 fontTheme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 16),
  axis.text.y = element_text(size = 16),
  axis.title.y = element_text(size = 18),
  legend.title = element_text(size=14),
  legend.text = element_text(size=13),
 legend.position="top")
 saveRDS(minMaxSummary, sprintf("minMaxSummary_%s.rds", fileidentifier))
 # line segments
 ggplot(minMaxSummary) + geom_segment(aes(x=kinit_FC_min, xend=kinit_FC_max, y=Ycoord, yend=Ycoord)) + scale_x_continuous(trans=log2_trans(),breaks=trans_breaks('log2', function(x) 2^x),labels=trans_format('log2', math_format(2^.x))) + xlab("segment joining min and max kinit FC for all gene") + ylab("vertical axis to separate intervals") + theme_bw() + fontTheme
 ggsave(sprintf("intervals_kinit_FC_%s.pdf", fileidentifier))
 ggplot(minMaxSummary) + geom_segment(aes(x=krel_FC_min, xend=krel_FC_max, y=Ycoord, yend=Ycoord)) + scale_x_continuous(trans=log2_trans(),breaks=trans_breaks('log2', function(x) 2^x),labels=trans_format('log2', math_format(2^.x))) + xlab("segment joining min and max kinit FC for all gene") + ylab("vertical axis to separate intervals") + theme_bw() + fontTheme
 ggsave(sprintf("intervals_krel_FC_%s.pdf", fileidentifier))


 # width of kinit FC and krel FC
 ggplot(minMaxSummary) + geom_violin(aes(x="width - kinit_FC", y = kinitFCRange)) + geom_point(position = dodge, alpha=0.5, aes(x="width - kinit_FC", y=kinitFCRange)) + xlab("") + ylab(sprintf("(max - min) kinit FC for all genes - %s", exprName)) + theme_bw()
  ggsave(sprintf("Width_kinit_FC_%s.pdf", fileidentifier))
 ggplot(minMaxSummary) + geom_violin(aes(x="width - krel_FC", y = krelFCRange)) + geom_point(position = dodge, alpha=0.5, aes(x="width - krel_FC", y=krelFCRange)) + xlab("") + ylab(sprintf("(max - min) krel FC for all genes - %s", exprName)) + theme_bw()
  ggsave(sprintf("Width_krel_FC_%s.pdf", fileidentifier))
 # save summary
 summaryRateData = group_by(rateData, gene) %>% summarise_at(c("kinit_base", "kinit_expr", "krel_base", "krel_expr", "kinit_FC", "krel_FC","kelong", "kpre"), mean, na.rm=T)
 #summaryRateData$kinit_FC = summaryRateData$kinit_expr / summaryRateData$kinit_base
 #summaryRateData$krel_FC = summaryRateData$krel_expr / summaryRateData$krel_base
 saveRDS(summaryRateData, sprintf("summary_repeatedParamScanData_%s.rds", fileidentifier))
 # savefrequency table
 freqTableGenes = table(rateData$gene)
 saveRDS(freqTableGenes, sprintf("freqTable_repeatedParamScans_%s.rds", fileidentifier))
 uniqueGenes = unique(rateData[,"gene"])
 # start plotting
 uniqGeneNum = min(uniqGeneNum, length(uniqueGenes))
 print(sprintf("plotting data for first %d genes", uniqGeneNum))
 subsetData = rateData[rateData$gene %in% uniqueGenes[1:uniqGeneNum],]
 ## Mean kinit_FC and krel_FC
 ggplot(summaryRateData) + geom_violin(aes(x="Initiation", y=kinit_FC)) + geom_violin(aes(x="Pause Release", y=krel_FC)) +
	       geom_boxplot(aes(x="Initiation", y=kinit_FC), alpha=0.9) + geom_boxplot(aes(x="Pause Release", y=krel_FC),alpha=0.9) +  
	       geom_point(position = "jitter", alpha=ALPHA, aes(x="Initiation", y=kinit_FC)) +
	        geom_point(position = "jitter", alpha=ALPHA, aes(x="Pause Release", y=krel_FC)) +
		scale_y_continuous(trans=log2_trans(),breaks=trans_breaks('log2', function(x) 2^x),labels=trans_format('log2', math_format(.x))) +
		#scale_y_continuous(breaks = seq(0, ceiling(max(max(summaryRateData[,"kinit_FC"]), max(summaryRateData[,"krel_FC"]))), len = 10)) +
	        ylab(sprintf("log2 (Mean Fold Change) for each gene - %s", yLabel)) + xlab("Rates") + theme_bw()  + fontTheme
  ggsave(sprintf("log2_Mean_foldchange_kinit,krel_%s.pdf", fileidentifier))

 
  # log2 TRANSFORMATION - plot y vals
#  ggplot(summaryRateData) + geom_violin(aes(x="kinit", y=kinit_FC)) + 
#	  geom_point(position=dodge, alpha=0.5, aes(x="kinit",y=kinit_FC)) + 
#	  scale_y_continuous(trans=log2_trans(),breaks=trans_breaks('log2', function(x) 2^x),labels=trans_format('log2', math_format(2^.x))) +
#           ylab(sprintf("Mean value of Fold Change for each gene - %s", yLabel)) + xlab("Rates") + theme(text = element_text(size=12)) + theme_bw()
#  ggsave(sprintf("log2_Mean_foldchange_kinit_%s.pdf", fileidentifier))
# 
# ggplot(summaryRateData) + geom_violin(aes(x="krel", y=krel_FC)) + 
#	  geom_point(position=dodge, alpha=0.5, aes(x="krel",y=krel_FC)) + 
#	  scale_y_continuous(trans=log2_trans(),breaks=trans_breaks('log2', function(x) 2^x),labels=trans_format('log2', math_format(2^.x))) +
#           ylab(sprintf("Mean value of Fold Change for each gene - %s", yLabel)) + xlab("Rates") + theme(text = element_text(size=12)) + theme_bw()
#  ggsave(sprintf("log2_Mean_foldchange_krel_%s.pdf", fileidentifier))
#

   ## both kinit fc and krel fc on same plot for one single gene
  #ggplot(subsetData) + geom_violin(aes(x="kinit_fc", y=kinit_FC)) + geom_point(position = position_jitter(width= 0.05), alpha=0.05,  aes(x="kinit_fc", y=kinit_FC)) + geom_violin(aes(x="krel_fc", y=krel_FC)) + geom_point(position = position_jitter(width= 0.05), alpha=0.05,  aes(x="krel_fc", y=kinit_FC)) + ylim(0,NA) + ylab("fold change for TMEM221 over multi param estimation")
  ## KINIT FC
 ggplot(subsetData) + geom_violin(aes(x=gene, y= kinit_FC, group=gene)) + geom_point(position = position_jitter(width= 0.05), alpha=0.05, aes(gene, y=kinit_FC, group=gene)) + ylab(sprintf("kinit foldchange - %s vs %s", exprName, baselineName)) + theme_bw()
ggsave(sprintf("kinit_FC_%s.png", fileidentifier), width=14,dpi=300)
## KREL FC
 ggplot(subsetData) + geom_violin(aes(x=gene, y= krel_FC)) + geom_point(position = position_jitter(width= 0.05), alpha=0.05, aes(gene, y=krel_FC, group=gene)) + ylab(sprintf("krel foldchange - %s vs %s", exprName, baselineName))+ theme_bw()
ggsave(sprintf("krel_FC_%s.png", fileidentifier), width=14,dpi=300)
## KINIT - EXPeriment
 ggplot(subsetData) + geom_violin(aes(x=gene, y= kinit_expr)) + geom_point(position = position_jitter(width= 0.05), alpha=0.05, aes(gene, y=kinit_expr, group=gene)) + ylab(sprintf("kinit - %s", exprName))+ theme_bw()
 ggsave(sprintf("kinit_expr_%s.png", fileidentifier), width=14,dpi=300)
## KINIT - baseline
 ggplot(subsetData) + geom_violin(aes(x=gene, y= kinit_base)) + geom_point(position = position_jitter(width= 0.05), alpha=0.05, aes(gene, y=kinit_base, group=gene)) + ylab(sprintf("kinit - %s", baselineName))+ theme_bw()
 ggsave(sprintf("kinit_base_%s.png", fileidentifier), width=14,dpi=300)
## KREL - EXPeriment
 ggplot(subsetData) + geom_violin(aes(x=gene, y= krel_expr)) + geom_point(position = position_jitter(width= 0.05), alpha=0.05, aes(gene, y=krel_expr, group=gene)) + ylab(sprintf("krel - %s", exprName))+ theme_bw()
 ggsave(sprintf("krel_expr_%s.png", fileidentifier), width=14,dpi=300)
## KREL - baseline
 ggplot(subsetData) + geom_violin(aes(x=gene, y= krel_base)) + geom_point(position = position_jitter(width= 0.05), alpha=0.05, aes(gene, y=krel_base, group=gene)) + ylab(sprintf("krel - %s", baselineName))+ theme_bw()
 ggsave(sprintf("krel_base_%s.png", fileidentifier), width=14,dpi=300)
## KELONG
 ggplot(subsetData) + geom_violin(aes(x=gene, y= kelong)) + geom_point(position = position_jitter(width= 0.05), alpha=0.05, aes(gene, y=kelong, group=gene)) + ylab(sprintf("kelong - %s, %s", exprName,baselineName))+ theme_bw()
 ggsave(sprintf("kelong_%s.png", fileidentifier), width=14,dpi=300)
}



