source("/Users/rudradeepmukherjee/Documents/UConn Health/comparment.model.project/Figures_manuscript/generalPlotter.R")

## please generate this data using the steps mentioned in the vignette - these steps are long, so omitted here for brevity
scanDataFile.entireRange = "/Users/rudradeepmukherjee/Documents/UConn Health/comparment.model.project/datasets, genome/Jonkers_2014_Trp_FP_paper/ParamEstimation_FP_TRP/paramFit_TRP_5kb_body_percentile/calculatedRateData_TRP_125min_repressed.rds"
summmary.TRP.entireRange = summaryFunction(readRDS(scanDataFile.entireRange))
summmary.TRP.entireRange$datasetName = "Triptolide (TRP)"

## Figure 2
summmary.TRP.entireRange$datasetName = "Triptolide (TRP)"
quantile(summmary.TRP.entireRange$B_FC/ summmary.TRP.entireRange$P_FC, probs = c(0.1,0.5,0.9))

data.df = rbind( summmary.TRP.entireRange[,colnames(summary.FP.entireRange)], summary.FP.entireRange)
#data.df$datasetName = factor(data.df$datasetName, levels = c("ZNF143 dTAG", "TBP dTAG"))
data.df$datasetName = factor(data.df$datasetName, levels = c("Triptolide (TRP)","Flavopiridol (FP)"))

#labelData = c(expression(atop("ZNF143 dTAG", "(HEK293T, 30 mins)")), expression(atop("TBP dTAG", "(HAP1, 2 hrs)")))
labelData = c(expression(atop("Triptolide", "(mESCs, 12.5 mins)")), expression(atop("Flavopiridol", "(mESCs, 25 mins)")))

krelFCplotter_multidatasets(data.df, "Fig2_krelFC_TRP_FP_new",  labelData, Y_DOWN_arr = c(rep(1/4,4)), Y_UP_arr = c(rep(2^2.2,4)), H_JUST=0.5,  showPercentLabels = TRUE, showPoints = FALSE, filetype="pdf", WIDTH=8, LOG_trans_break = 10^(seq(-2,1.8,1)), LOG_trans_format = scales::label_number(accuracy = 0.01))

TRP.FP.df = kinitFCvsPfoldchangePlotter_multidatasets(data.df, "Fig2_kinitFC_TRP_FP_new", labelData, Y_DOWN_arr = rep(1/5,3),  Y_UP_arr = rep(2,3), H_JUST_LEFT = 1.2, H_JUST_MIDDLE=0.5, H_JUST_RIGHT=-0.2, showPercentLabels = TRUE,  addKpreEqKrel = T, legend.colors = c("#36454F", "#4169E1", "#D2B48C"), diff_UP = -1.9, diff_DOWN = 0.1, WIDTH= 14, HEIGHT=7.2, filetype="pdf", LOG_trans_break = 10^(seq(-2,1.8,1)), LOG_trans_format = scales::label_number(accuracy = 0.01))

quantile(subset(summmary.TRP.entireRange, WithinRange == TRUE)$Z_factor, probs = c(0.1,0.5,0.9) )

## Figure 3
## TRP - ratio of premature termination to pause release
df = prepareDf.w.Kpre_by_Krel(summmary.TRP.entireRange, kinitFC=0.25)
nrow(subset(df, WithinRange == T))
TRP.df = KpreByKrelPlotter(summmary.TRP.entireRange, "Fig3_TRP_kpre_by_krel_Factor.pdf", plotTitle = "TRP (3326 genes)", kinit_FC=0.25, WIDTH=5)

## plot transcription output with changes in FC(kinit), FC(krel)
transcriptionOutPutvsFC_kinit_krel("Kinit_percentage_vs_B_Foldchange_legend_top.pdf")

## Figure 4
summary.HS.20min.entireRange = summaryFunction(readRDS("Duarte_Drosop_HSF/ParamEstimation/paramFit_HS_20min_HSF_targets/calculatedRateData_HS_20min_HSF_activated_TablS2.rds"))
summary.HS.20min.entireRange$datasetName = "HSF Activated (S2)"
summary.Dex100nm.entireRange$datasetName = "Dex, 100nm"
summary.TBP.entireRange$datasetName = "TBP dTAG"
summary.ZNF143.entireRange$datasetName = "ZNF143 dTAG"
data.df = rbind(summary.TBP.entireRange, summary.ZNF143.entireRange, summary.Dex100nm.entireRange,
                #summary.Celastrol.entireRange, summary.Vihervara_HS.entireRange,
                summary.HS.20min.entireRange)
#data.df$datasetName = factor(data.df$datasetName, levels = c("TBP dTAG", "ZNF143 dTAG", "Dex, 100nm", "HS, 30min (K562)", "HS, 20min (S2)"))
data.df$datasetName = factor(data.df$datasetName, levels = c("TBP dTAG", "ZNF143 dTAG", "Dex, 100nm",
                                                             #"Celastrol (K562)","HS, 30min (K562)",
                                                             "HSF Activated (S2)"))
labelData = c(expression(atop("TBP dTAG", "(HAP1, 2 hrs)")),
              expression(atop("ZNF143 dTAG", "(HEK293T, 30 mins)")) ,
              expression(atop("Dex 100nm", "(A549, 45 mins)")),
              # expression(atop("Celastrol", "(K562, 20 mins)")),
              #expression(atop("Heat Shock", "(K562, 30 mins)")),
              expression(atop("HSF Activated", "(S2, HS 20 mins)")))
              #expression(atop("Heatshock", "(S2, 20 mins)")))

krelFCplotter_multidatasets(data.df, "Fig3_krelFC_TBP_ZNF143_Dex_HS_Droso",  labelData, Y_DOWN_arr = rep(1/5,4), Y_UP_arr = rep(2^3.4,4), H_JUST=0.5, filetype="pdf", showPoints=F, WIDTH=13.5, LOG_trans_break = 10^(seq(-1,1.5,0.5)), LOG_trans_format = scales::label_number(accuracy = 0.01))

##without kpre << krel
Fig3_kinit.df = kinitFCvsPfoldchangePlotter_WoutKrelFaster(data.df, "Fig3_kinitFC_TBP_ZNF143_Dex_HS_Droso_woutKrelFaster", labelData, Y_DOWN_arr = rep(1/5,4),  Y_UP_arr =  rep(2^3.8,4), filetype="pdf", H_JUST_LEFT = 0.5, H_JUST_RIGHT=0.5, showPercentLabels = TRUE, legend.colors = c("#36454F", "#D2B48C"), diff_UP = 19, diff_DOWN = 0.2, WIDTH=13.5, VIOLIN_WIDTH=0.8,  LOG_trans_break =  10^(seq(-1,1.5,0.5)), LOG_trans_format = scales::label_number(accuracy = 0.01))

## Figure 5
# krel
only.Dex.df = kelongVsKrelPlotter_multiDatasets_twoPanels(summary.Dex100nm.entireRange, "Fig4_Dex_Krel_rates", titleText = "GR-activated", WIDTH=4, HEIGHT=8.5, addTrends=T, LOG_trans_break = 10^(seq(-1,2,1)), LOG_trans_format = scales::label_number(accuracy = 0.01))
only.Dex.effective.kreldf = kelongVsEffectiveKrelPlotter_multiDatasets_twoPanels(summary.Dex100nm.entireRange, "Fig4_Dex_Effective_Krel", titleText = "GR-activated", WIDTH=3.5, HEIGHT=8.5, widthMultiplier = 1.15, addTrends=T, LOG_trans_break = 10^(seq(-1,1,1)), Occupancy_LOG_trans_break =  10^(seq(2,5,1)),  LOG_trans_format = scales::label_number(accuracy = 0.01))
# kinit
summary.TBP.entireRange$datasetName = "TBP degradation"
#summary.ZNF143.entireRange$datasetName = "ZNF143 dTAG (HEK293T, 30 mins)"
summary.ZNF143.entireRange$datasetName = "ZNF143 degradation"
data.df = rbind(summary.TBP.entireRange, summary.ZNF143.entireRange)
TBP.ZNF143.kinit.df = plotKinitValues_with_Kpre_Krel(data.df, "Fig_4_ZNF143_TBP_kinit_rates_indvKpre", experiment="TRP", factorMultiplier = 2.6, addTrends= T , ConstantKpre = NULL,globalMaxKpre = FALSE, HEIGHT=9, WIDTH = 9.3, LOG_trans_break = 10^(seq(-2,1,0.5)), LOG_trans_format = scales::label_number(accuracy = 0.01))
# calculate median
only.Dex.df = kelongVsKrelPlotter_multiDatasets_twoPanels(summary.Dex100nm.entireRange, "Fig4_Dex_Krel_rates", titleText = "GR-activated", WIDTH=4, HEIGHT=8.5, addTrends=T, LOG_trans_break = 10^(seq(-1,2,1)), LOG_trans_format = scales::label_number(accuracy = 0.01))
quantile(subset(only.Dex.df, cond=="control")$krel, probs = c(0.1,0.5,0.9))
quantile(subset(only.Dex.df, cond=="treatment")$krel, probs = c(0.1,0.5,0.9))
only.TBP.df = kelongVsKrelPlotter_multiDatasets_twoPanels(summary.TBP.entireRange, "Fig4_TBP_Krel_rates", titleText = "TBP-drag", WIDTH=4, HEIGHT=8.5, addTrends=T, LOG_trans_break = 10^(seq(-1,2,1)), LOG_trans_format = scales::label_number(accuracy = 0.01))
quantile(subset(only.TBP.df, cond=="control")$krel, probs = c(0.1,0.5,0.9))
only.TBP..df = kelongVsEffectiveKrelPlotter_multiDatasets_twoPanels(summary.TBP.entireRange, "Fig4_TBP_Effective_Krel", titleText = "TBP-dtag", WIDTH=5.5, HEIGHT=8.5, widthMultiplier = 0.9, addTrends=T, LOG_trans_break = 2^(seq(-6,2,2)), LOG_trans_format = scales::label_number(accuracy = 0.01))
quantile(subset(only.TBP..df, cond=="control")$effectiveKrel, probs = c(0.1,0.5,0.9))
only.ZNF143.df = kelongVsKrelPlotter_multiDatasets_twoPanels(summary.ZNF143.entireRange, "Fig4_ZNF143_Krel_rates", titleText = "znf143-drag", WIDTH=4, HEIGHT=8.5, addTrends=T, LOG_trans_break = 10^(seq(-1,2,1)), LOG_trans_format = scales::label_number(accuracy = 0.01))
quantile(subset(only.ZNF143.df, cond=="control")$krel, probs = c(0.1,0.5,0.9))
only.ZNF143..df = kelongVsEffectiveKrelPlotter_multiDatasets_twoPanels(summary.ZNF143.entireRange, "Fig4_ZNF143_Effective_Krel", titleText = "ZNF143-dtag", WIDTH=5.5, HEIGHT=8.5, widthMultiplier = 0.9, addTrends=T, LOG_trans_break = 2^(seq(-6,2,2)), LOG_trans_format = scales::label_number(accuracy = 0.01))
merge.df = merge(summary.ZNF143.entireRange[,c("gene", "P_base_mean")], subset(only.ZNF143.df, cond=="control"), by = "gene")
merge.df$effectiveKrel = merge.df$P_base_mean * merge.df$krel
quantile(merge.df$effectiveKrel, probs = c(0.1,0.5,0.9))
# figure 6, halflife
summmary.TRP.entireRange$datasetName = "TRP 12.5mins (mESCs)"
summmary.FP.entireRange$datasetName = "FP 25mins (mESCs)"
summary.ZNF143.entireRange$datasetName = "ZNF143 dTAG (HEK293T)"
summary.TBP.entireRange$datasetName = "TBP dTAG (HAP1)"
summary.Dex100nm.entireRange$datasetName = "Dex 45mins (A549)"
summary.Celastrol.entireRange$datasetName = "Celastrol (K562)"
summary.Vihervara_HS.entireRange$datasetName = "HS 30mins (K562)"
summary.HS.20min.entireRange$datasetName = "HSF Activated (S2)"
labelData = c(expression("TRP 12.5mins (mESCs)"), expression("FP 25mins (mESCs)"), expression("ZNF143 dTAG (HEK293T)"), expression("TBP dTAG (HAP1)") , expression("Dex 45mins (A549)"),
              #expression("Celastrol (K562)"), expression("HS 30mins (K562)"),
              expression("HSF Activated (S2)"))
data.df = rbind(summmary.TRP.entireRange[,colnames(summmary.FP.entireRange)], summmary.FP.entireRange, summary.ZNF143.entireRange, summary.TBP.entireRange, summary.Dex100nm.entireRange,
                    #summary.Celastrol.entireRange, summary.Vihervara_HS.entireRange,
                  summary.HS.20min.entireRange)
data.df$datasetName = factor(data.df$datasetName, levels = labelData)
halflife.df = halfLifePlotter_indv_Kpre_forGenes_multidatasets(data.df, "Fig_Halflife_indv_Kpre_TRP_to_HS", labelData, nROW = 2, factorMultipier = 2.6, useControlOnly = T, log_scale=T, HEIGHT = 9, WIDTH=13.5, LOG_trans_break =  10^(seq(-1,3,1)), LOG_trans_format = scales::label_number(accuracy = 0.1))
halflife.df = halfLifePlotter_indv_Kpre_forGenes_multidatasets(data.df, "Fig_Halflife_indv_Kpre_TRP_to_HS_factor_10", labelData, nROW = 2, factorMultipier = 10, useControlOnly = T, log_scale=T, HEIGHT = 9, WIDTH=13.5,  LOG_trans_break = 10^(seq(-1,3,1)), LOG_trans_format = scales::label_number(accuracy = 0.1))

## supplementary figures begin
# Fig S1
## Dex trends - manually curate gene lists
activated.expectedFC = readLines("activatedLongGenes__gt_150kb_with_ExpectedFoldchange.txt")
wrongTranscript = readLines("wrongTranscript_activated_longGenes.txt")
badSignal = readLines("badSignal_activated_longGenes.txt")
# for activated genes
data.df.control = readRDS("Dex.control.windowSize_10kb_step_1kb.rds")
data.df.Dex = readRDS("Dex.treatment.windowSize_10kb_step_1kb.rds")
plotTrendlineForGeneSets(data.df.control, data.df.Dex, geneList = activated.expectedFC,
                        trendLineGeneList = activated.expectedFC ,
                         plotTitle = "151 Dex-activated genes (len >= 150kb, trendline for 85 genes with expected behavior)", 
                         filename = "Dex_FC_genes_expected_transition.pdf")
## for two genes
plotTrendlineForGeneSets(data.df.control, data.df.Dex, geneList = activated.expectedFC,
                        trendLineGeneList = c("ABTB3", "FLNB"), 
                         plotTitle = "", showLegend = T,
                         filename = "Dex_FC_genes_expected_transition_ABTB3_FLNB.pdf")
#
# for matched genes
matched.unchanged.long.genes = readLines("Dex_matched_unchanged_longgenes_150kb.txt")
data.df.control = readRDS("Dex.Matched.Genes.control.windowSize_10kb_step_1kb.rds")
data.df.Dex = readRDS("Dex.Matched.Genes.treatment.windowSize_10kb_step_1kb.rds")
plotTrendlineForGeneSets(data.df.control, data.df.Dex, geneList = matched.unchanged.long.genes,
                         trendLineGeneList = matched.unchanged.long.genes,
                         plotTitle = "85 Matched Unchanged genes, Dex (len >= 150kb)",
                         filename = "Dex_FC_matched_unchanged_genes.pdf")

## Fig S1 - TRP fit
#We use `fitKinitToNilsonData.py` to find parameters for fits to Fig 6c in Nilson et al 2017 and save the parameters to  `Nilson2017_Kinit_fit_and_deriv.txt`  
# we then use `plotFitsForNilson2017.R` for plotting the fits
cat fitKinitToNilsonData.py
cat plotFitsForNilson2017.R

## Fig S2 - saturation curves for all datasets
source("residualAndFitPlotter.R")
#TRP
rootDir.TRP = "/Users/rudradeepmukherjee/Documents/UConn Health/comparment.model.project/datasets, genome"
paramFileLoc.TRP = sprintf("%s/Jonkers_2014_Trp_FP_paper/TRP_LOG_HYPER_TANGENT_saturation_Parameters.txt", rootDir.TRP)
originalPauseSums.TRP =  read.table(sprintf("%s/Jonkers_2014_Trp_FP_paper/TRP125min_control_allgenes_varlessthanpoint1.txt", rootDir.TRP), sep="\t", header=T)
color.vals = c("#F69994", "#F5807B", "#F26660", "#ED4843", "#E92623")
fitPlotter_multiplePercentiles("TRP", paramFileLoc.TRP, originalPauseSums.TRP, percentileVals = c(70,80,90,95,100), "LOG_HYPER_TANG", "No", plotTitle = "Triptolide (mESCs, no UMIs)", logYPlotter="Yes", showLegend = "Yes", plotMaxSaturation="Yes", outputDir = sprintf("%s/Jonkers_2014_Trp_FP_paper/fittedVals_R", rootDir.TRP), V_JUST = c(1.2,-0.8, -0.7,-1.1, -1.1), H_JUST=c(rep(0.6,2),rep(0.58,2), 0.545), Label_YCoord = 10^2.8, HJUST_LABEL = 0.6, color.list = color.vals)
## FP
rootDir.FP = "/Users/rudradeepmukherjee/Documents/UConn Health/comparment.model.project/datasets, genome"
paramFileLoc.FP = sprintf("%s/Jonkers_2014_Trp_FP_paper/FP_LOG_HYPER_TANGENT_saturation_Parameters.txt", rootDir.FP)
originalPauseSums.FP =  read.table(sprintf("%s/Jonkers_2014_Trp_FP_paper/FP25min_control_allgenes_varlessthanpoint1.txt", rootDir.FP), sep="\t", header=T)
fitPlotter_multiplePercentiles("FP", paramFileLoc.FP, originalPauseSums.FP, percentileVals = c(70,80,90,95,100), "LOG_HYPER_TANG", "No", plotTitle = "Flavopiridol (mESCs, no UMIs)", logYPlotter="Yes", showLegend = "Yes", plotMaxSaturation="Yes", outputDir = sprintf("%s/Jonkers_2014_Trp_FP_paper/fittedVals_R", rootDir.FP), V_JUST = c(1.3,-0.6, -0.7,-1.1, 1.3), H_JUST=c(rep(0.6,3),0.58, 0.56), Label_YCoord = 10^3.4, HJUST_LABEL = 0.63, color.list = color.vals)
## S2 cells
rootDir = "/Users/rudradeepmukherjee/Documents/UConn Health/comparment.model.project/datasets, genome"
paramFileLoc.concat = "/Users/rudradeepmukherjee/Documents/UConn Health/comparment.model.project/datasets, genome/Duarte_Drosop_HSF/Duarte_Droso_HS_LOG_HYPER_TAN_saturation_Parameters_concat.txt"
combn.pause.Droso.HS = tibble(read.table("/Users/rudradeepmukherjee/Documents/UConn Health/comparment.model.project/datasets, genome/Duarte_Drosop_HSF/Duarte_Droso_HS_concat_control_varlesspoint1.txt", sep="\t", header=T))
combn.pause.Droso.HS = subset(combn.pause.Droso.HS, pause_sum > 0 )
combn.pause.Droso.HS = combn.pause.Droso.HS %>% arrange(pause_sum)
fitPlotter_multiplePercentiles("Droso_HS_concat", paramFileLoc.concat, combn.pause.Droso.HS, percentileVals = c(70,80,90,95,99.99), "LOG_HYPER_TANG", "No", plotTitle = "S2 (control + HS, no UMIs)", logYPlotter="Yes", showLegend = "No", plotMaxSaturation="No", outputDir = sprintf("%s/Duarte_Drosop_HSF/fittedVals", rootDir), V_JUST = c(1.2,-0.3, -0.3, -0.3, -0.3), H_JUST=c(rep(0.58,5)), Label_YCoord = 10^4, HJUST_LABEL = 0.6, color.list = color.vals)
# ZNF143
rootDir.ZNF143 = "/Users/rudradeepmukherjee/Documents/UConn Health/comparment.model.project/datasets, genome"
paramFileLoc.ZNF143 = sprintf("%s/ZNF143_Data/ZNF143_dTAG_LOG_HYPER_TANGENT_saturation_Parameters.txt", rootDir.ZNF143)
originalPauseSums.ZNF143 =  read.table(sprintf("%s/ZNF143_Data/ZNF143dTAGv1_baseline_allgenes_varlessthanpoint1.txt", rootDir.ZNF143), sep="\t", header=T)
originalPauseSums.ZNF143 = subset(originalPauseSums.ZNF143 , pause_sum > 0)
originalPauseSums.ZNF143 = originalPauseSums.ZNF143[rank(order(originalPauseSums.ZNF143$pause_sum)),]
fitPlotter_multiplePercentiles("ZNF143", paramFileLoc.ZNF143, originalPauseSums.ZNF143, percentileVals = c(70,80, 90,95,100), "LOG_HYPER_TANG", "No", plotTitle = "ZNF143 dTAGv1 (HEK293T, UMI=8N)", logYPlotter="Yes", plotMaxSaturation="Yes", showLegend = TRUE, outputDir = sprintf("%s/ZNF143_Data/fittedVals_R", rootDir.ZNF143), V_JUST = c(-5,-5.4, -5.6,-6.3, -6), Label_YCoord = 10^3.4, HJUST_LABEL = 0.6, H_JUST=rep(0.5,5), color.list = color.vals, legendPosition = c(0.8,0.2), WIDTH=10)
# Dexamethasone
rootDir.Dex100nm = "/Users/rudradeepmukherjee/Documents/UConn Health/comparment.model.project/datasets, genome"
paramFileLoc.Dex100nm = sprintf("%s/ErinWissink_et_al_2021/Dex_LOG_HYPER_TANGENT_saturation_Parameters.txt", rootDir.Dex100nm)
originalPauseSums.Dex100nm =  read.table(sprintf("%s/ErinWissink_et_al_2021/ParamEstimation/ErinWissink_allgenes_Dex100nm_baseline.txt", rootDir.Dex100nm), sep="\t", header=T)
#originalPauseSums.Dex100nm =  read.table(sprintf("%s/ErinWissink_et_al_2021/Dex100nm_scaled_2071.909116_filtered_var_b_mean_lessthan_0.100000_baseline.txt", rootDir.Dex100nm), sep="\t", header=T)
originalPauseSums.Dex100nm = subset(originalPauseSums.Dex100nm, pause_sum > 0)

fitPlotter_multiplePercentiles("Dex100nm", paramFileLoc.Dex100nm, originalPauseSums.Dex100nm, percentileVals = c(70,80, 90,95,100), "LOG_HYPER_TANG", "No", plotTitle = "Dexamethasone (A549, UMI=6N)", logYPlotter="Yes", showLegend = TRUE, plotMaxSaturation="Yes", outputDir = sprintf("%s/ErinWissink_et_al_2021/fittedVals_R", rootDir.Dex100nm), V_JUST = c(-3.6,-4.2, -4.8,-5.4, -5.9), H_JUST=c(rep(0.5,4),0.49), Label_YCoord = 10^2.9, HJUST_LABEL = 0.6, color.list = color.vals, WIDTH=10)

# TBP dtag v1
rootDir.TBP = "/Users/rudradeepmukherjee/Documents/UConn Health/comparment.model.project/datasets, genome"
paramFileLoc.TBP = sprintf("%s/price2022TBP_paper/TBP_LOG_HYPER_TANGENT_saturation_Parameters.txt", rootDir.TBP)
originalPauseSums.TBP =  read.table(sprintf("%s/price2022TBP_paper/TBP_dTAGv1_control_Allgenes_var_lessthanpoint1.txt", rootDir.TBP), sep="\t", header=T)
originalPauseSums.TBP = originalPauseSums.TBP[rank(order(originalPauseSums.TBP$pause_sum)),]

fitPlotter_multiplePercentiles("TBP", paramFileLoc.TBP, originalPauseSums.TBP, percentileVals = c(70,80, 90,95,100), "LOG_HYPER_TANG", "No", plotTitle = "TBP dTAGv1 (HAP1, UMI=4N)", logYPlotter="Yes", plotMaxSaturation="Yes", showLegend = TRUE, outputDir = sprintf("%s/price2022TBP_paper/fittedVals_R", rootDir.TBP), V_JUST = c(-2,-3.1, -3.8,-4.5, -5) , H_JUST=rep(0.5,5), Label_YCoord = 10^3.4, HJUST_LABEL = 0.6, color.list = color.vals, legendPosition = c(0.8,0.2), WIDTH=10)

## Fig S3 - UMI vs no UMI with pause sems for datasets
## Fig S3A
## these values are from Fig S2
minMaxSaturation.df = rbind(c("TRP", 50, 1154, "Absent"), c("FP", 210, 5212, "Absent"),
                            c("S2 HS", 126, 75721, "Absent"),
                            # c("K562 HS", 36, 62, "Absent"),
                            c("TBP dTAG", 217, 344, "Present"), c("ZNF143 dTAG", 58, 236, "Present"), c("Dex", 43, 105, "Present"))
minMaxSaturation.df = as.data.frame(minMaxSaturation.df)
colnames(minMaxSaturation.df) = c("dataset", "minSaturation", "maxSaturation", "UMI")
minMaxSaturation.df$dataset = factor(minMaxSaturation.df$dataset, c("TRP", "FP", "S2 HS", "Dex", "TBP dTAG", "ZNF143 dTAG"))
minMaxSaturation.df$minSaturation = as.numeric(minMaxSaturation.df$minSaturation)
minMaxSaturation.df$maxSaturation = as.numeric(minMaxSaturation.df$maxSaturation)
minMaxSaturation.df$foldchange = minMaxSaturation.df$maxSaturation / minMaxSaturation.df$minSaturation
ggplot(minMaxSaturation.df, aes(x = UMI, y = foldchange)) +
        #geom_boxplot() +
        geom_point(aes(color = dataset), size=5) +
        scale_y_continuous(trans=log10_trans(),
                breaks=trans_breaks('log10', function(x) 10^x),
                labels=trans_format('log10', math_format(.x)), limits=c(1,NA)) +
        ylab(expression(log[10]~frac("saturation value at 100%", "saturation value at 70%"))) +
        theme_bw() + fontTheme +
        theme(axis.title.x = element_text(size = 30))

# Fig S3B
root.path = "/Users/rudradeepmukherjee/Documents/UConn Health/comparment.model.project/datasets, genome/Duarte_Drosop_HSF"
minus.control = read.table(sprintf("%s/Duarte_control.Hsp70Aa_minus_sizeF.bedGraph", root.path), sep="\t")
minus.HS = read.table(sprintf("%s/Duarte_HS_20min.Hsp70Aa_minus_sizeF.bedGraph", root.path), sep="\t")

plus.control = read.table(sprintf("%s/Duarte_control.Hsp70Aa_plus_sizeF.bedGraph", root.path), sep="\t")
plus.HS = read.table(sprintf("%s/Duarte_HS_20min.Hsp70Aa_plus_sizeF.bedGraph", root.path), sep="\t")

colnames(minus.control) = c("chr", "start", "end", "value")
colnames(minus.HS) = c("chr", "start", "end", "value")
minus.control$condition = "control"
minus.HS$condition = "HS, 20min"

colnames(plus.control) = c("chr", "start", "end", "value")
colnames(plus.HS) = c("chr", "start", "end", "value")
plus.control$condition = "control"
plus.control.nospike = plus.control[c(-527),]
plus.HS.nospike = plus.HS[c(-567),]
data.df = rbind(plus.control.nospike, plus.HS.nospike)
data.df$value = data.df$value / 5  # for Hsp70 gene


ggplot(data.df, aes(x = start, y= value)) + geom_line() +
                      geom_point(shape=1) +
                      geom_vline(xintercept = c(215, 2300), color = "#0071BC", linetype = "dashed", size = 1.5) +
                      geom_vline(xintercept = c(75, 125), color = "#E58C2C", linetype = "dashed", size = 0.5) +
                      facet_wrap(~condition, ncol=1) +
                      ggtitle("Alignment to Drosophila Hsp70 consensus") +
                      xlab("distance from (TSS - 52bp)") +
                      #xlab("distance from gene start") +
                      ylab("signal") +
                      theme_bw() + fontTheme +
                      theme(strip.text = element_text(size=24),
                          strip.background = element_rect(fill = "#E6E6FA", color = "black", linewidth = 1.6),
                      plot.title = element_text(size=23, face=1.2,hjust=0.5))
 ggsave("Hsp70_composite_plus_strand_divide_by_5_RasmussenCoord.pdf")

# Fig S3C
percentile_90th = 182.36
Hsp70_HS_pausesum = 992
rootpath = "/Users/rudradeepmukherjee/Documents/UConn Health/comparment.model.project/datasets, genome/Duarte_Drosop_HSF"
control.pausesums = read.table(sprintf("%s/Duarte_Droso_TableS2_HSF_activated_control.txt", rootpath), sep="\t", header=T)
HS.pausesums = read.table(sprintf("%s/Duarte_Droso_TableS2_HSF_activated_treat.txt", rootpath), sep="\t", header=T)
control.pausesums$cond = "control"
HS.pausesums$cond = "HS, 20min"
ggplot(rbind(control.pausesums, HS.pausesums), aes(x = cond, y = pause_sum)) +
        geom_violin(scale="width") +
        geom_boxplot(width = 0.3) +
        geom_hline(yintercept = c(percentile_90th), color = "#0071BC", linetype = "dashed", size = 1.5) +
        geom_text(aes(x = "HS, 20min", y = percentile_90th*0.65, label = "P, 90th percentile"), color = "#0071BC", size=6,family="mono", hjust=0.7, vjust=0) +
         geom_hline(yintercept = c(Hsp70_HS_pausesum), color = "#E58C2C", linetype = "dashed", size = 1.5) +
        geom_text(aes(x = "HS, 20min", y = Hsp70_HS_pausesum * 1.14, label = "P, Hsp70 consensus (HS)"), color = "#E58C2C", family="mono",size=6, hjust=0.7, vjust=0) +
        scale_y_continuous(trans=log10_trans(),
                breaks=trans_breaks('log10', function(x) 10^x),
                labels=trans_format('log10', math_format(.x)), limits=c(1,NA)) +
    ylab(expression(log[10]~"pause sum")) +
    ggtitle("HSF targets (S2 cells, no UMIs)") +
      theme_bw() + fontTheme +
      theme(axis.title.x = element_blank(),
      plot.title = element_text(size=28))
ggsave("HSF_targets_PauseSum_trend.pdf")

# Fig S4
newSummary.FP.Kpre.constrained = summary.FP.entireRange
newSummary.FP.Kpre.constrained$kinitFCcategory = "Undefined"
newSummary.FP.Kpre.constrained[newSummary.FP.Kpre.constrained$kinit_FC_max < 1, ]$kinitFCcategory = "kinitFC_less_1"
newSummary.FP.Kpre.constrained[newSummary.FP.Kpre.constrained$kinit_FC_max >= 1 & newSummary.FP.Kpre.constrained$kinit_FC_min <= 1, ]$kinitFCcategory = "kinitFC_span_

colorManual.two = c("#1F77B4", "#FF7F0E", "#E7298A", "#1B9E77") # orange, blue, dark palette
df.FP.kpreconstrain = as.data.frame(cbind(c("kinitFCless1", "kinitFCspan1"),
                   c(sum(newSummary.FP.Kpre.constrained$kinitFCcategory == "kinitFC_less_1"),
                     sum(newSummary.FP.Kpre.constrained$kinitFCcategory == "kinitFC_span_1"))))
                               #sum(newSummary.FP$kinitFCcategory == "kinitFC_more_1"))))
colnames(df.FP.kpreconstrain) = c("typeOfData", "countOfData")
df.FP.kpreconstrain$typeOfData = factor(df.FP.kpreconstrain$typeOfData, levels = c("kinitFCless1", "kinitFCspan1"))
df.FP.kpreconstrain$countOfData = as.numeric(df.FP.kpreconstrain$countOfData)
sum_all = sum(df.FP.kpreconstrain$countOfData)
df.FP.kpreconstrain$countOfData = df.FP.kpreconstrain$countOfData / sum_all

colors.New = c("#E7298A", "#1B9E77")
newSummary.FP.Kpre.constrained = summary.FP.entireRange
df.FP.pauserelease = as.data.frame(cbind(c("krelFCless1", "krelFCmore1"),
                   c(sum(newSummary.FP.Kpre.constrained$krel_FC_min > 1),
                     sum(newSummary.FP.Kpre.constrained$krel_FC_max <= 1))))
colnames(df.FP.pauserelease) = c("typeOfData", "countOfData")
df.FP.pauserelease$typeOfData = factor(df.FP.pauserelease$typeOfData, levels = c("krelFCless1", "krelFCmore1"))
df.FP.pauserelease$countOfData = as.numeric(df.FP.pauserelease$countOfData)
sum_all = sum(df.FP.pauserelease$countOfData)
df.FP.pauserelease$countOfData = df.FP.pauserelease$countOfData / sum_all
## fig S4A
ggplot(df.FP.pauserelease, aes(x = "Repressed\n(Flavopiridol)", y = countOfData, fill=typeOfData)) +
 geom_bar(stat="identity", position=position_fill(reverse=F)) +
         ylab("Fraction of Genes") +
         scale_fill_manual(
                           values = colors.New,
                           labels= c("Increase"~italic(K[rel]), "Decrease"~italic(K[rel]))) +
         theme_bw() + fontTheme +
         theme(legend.justification = "center",
               axis.text.x = element_text(size=31),
              axis.text.y = element_text(size=27),
               axis.title.y = element_text(size=32),
              axis.title.x = element_blank(), legend.title = element_blank())
ggsave("CountFrequency_stacked_KrelFC_profile_FP_treat_Kpre_constrained.pdf", height=8.4, width=8) # 6.53 x 8.4
## Fig S4B
newSummary.FP.Kpre.constrained = summary.FP.entireRange
newSummary.FP.Kpre.constrained$krelFCcategory = "Undefined"
newSummary.FP.Kpre.constrained[newSummary.FP.Kpre.constrained$krel_FC_max < 1, ]$krelFCcategory = "krelFC_less_1"
newSummary.FP.Kpre.constrained[newSummary.FP.Kpre.constrained$krel_FC_min >= 1, ]$krelFCcategory = "krelFC_more_1"
ggplot(newSummary.FP.Kpre.constrained,aes(x= krelFCcategory, y = P_FC) ) +
        geom_violin(scale="width") +
        geom_boxplot(width=0.3) +
        geom_hline(yintercept = 1, linetype=2) +
        scale_y_continuous(trans=log10_trans(),
          breaks=trans_breaks('log10', function(x) 10^x),
          labels=trans_format('log10', math_format(.x))) +
        scale_x_discrete(labels = c(expression(atop("decreasing",italic(K[rel]))), expression(atop("increasing",italic(K[rel]))))) +
        ylab(expression(log[10]~"foldchange in Pause Density")) +
                  theme_bw() + fontTheme +
                  theme(axis.title.x = element_blank())

ggsave("FP_treat_FC_P_Vs_krel_change.pdf")

## Fig S4C
FP.kpreConstrain.krelFCless1.PFCless1 = subset(newSummary.FP.Kpre.constrained, P_FC < 1 & krel_FC_max < 1)
FP.kpreConstrain.krelFCless1.PFCmore1 = subset(newSummary.FP.Kpre.constrained,  P_FC >= 1  & krel_FC_max < 1)

FP.kpreConstrain.krelFCmore1.PFCless1 = subset(newSummary.FP.Kpre.constrained, P_FC < 1 & krel_FC_min >= 1)
FP.kpreConstrain.krelFCmore1.PFCmore1 = subset(newSummary.FP.Kpre.constrained, P_FC >= 1 & krel_FC_min >= 1)

# NOTE - All unexpected genes (1480) showing increase in Krel, show decreasing kinit which means pause sum decreases upon FP treatment - which is wrong. Thus, these pause windows are wrong

df.FP.kpreConstrain.krel.P = as.data.frame(cbind(
                                                c("krelFCless1", "krelFCless1", "krelFCmore1", "krelFCmore1"),
                                                c("P_FCless1", "P_FCmore1", "P_FCless1", "P_FCmore1"),
                                         c(nrow(FP.kpreConstrain.krelFCless1.PFCless1), nrow(FP.kpreConstrain.krelFCless1.PFCmore1),
                                         nrow(FP.kpreConstrain.krelFCmore1.PFCless1), nrow(FP.kpreConstrain.krelFCmore1.PFCmore1) )))
colnames(df.FP.kpreConstrain.krel.P) = c("class", "typeOfData", "countOfData")
df.FP.kpreConstrain.krel.P$countOfData = as.numeric(df.FP.kpreConstrain.krel.P$countOfData)
df.FP.kpreConstrain.krel.P[df.FP.kpreConstrain.krel.P$class == "krelFCless1", ]$countOfData = df.FP.kpreConstrain.krel.P[df.FP.kpreConstrain.krel.P$class == "krelFCless1", ]$countOfData / (nrow(FP.kpreConstrain.krelFCless1.PFCless1) + nrow(FP.kpreConstrain.krelFCless1.PFCmore1))
df.FP.kpreConstrain.krel.P[df.FP.kpreConstrain.krel.P$class == "krelFCmore1", ]$countOfData = df.FP.kpreConstrain.krel.P[df.FP.kpreConstrain.krel.P$class == "krelFCmore1", ]$countOfData / (nrow(FP.kpreConstrain.krelFCmore1.PFCless1) + nrow(FP.kpreConstrain.krelFCmore1.PFCmore1))

colors.New = c("#E7298A", "#1B9E77")
ggplot(df.FP.kpreConstrain.krel.P, aes(x=class, y = countOfData, fill=typeOfData)) +
         #geom_bar(stat="identity", position=position_fill(reverse=T)) +
         geom_bar(stat="identity", position=position_stack(reverse=T), width=0.6) +
         #scale_y_continuous(labels = scales::percent) +
         ylab("Fraction of Genes") +
         #ylab(expression("Fraction of Genes (with" ~ italic(K[init]) ~ " decrease)")) +
         #xlab(expression("Genes with decreasing" ~ italic(K[init]) )) +
         #labsForLegend +
         #scale_fill_brewer(palette = "Dark2",
        scale_fill_manual(values = colorManual.two ,
                          #labels= c(expression("Decrease" ~ italic(K[init])), "Uncertain Change")) +
                          labels= c(expression("Decrease Pause Sum"), "Increase Pause Sum")) +
        #scale_fill_manual(values = colorManual.two,
        #              labels= c("Decrease Pause Release", "Increase Pause Release")) +
         #scale_fill_frontiers(labels= c("Decrease Pause Release", "Increase Pause Release")) +
        #labs(fill = "Gene Categories") +
        scale_x_discrete(labels = c(expression(atop("Decrease","in" ~ italic(K[rel]))),
                                    expression(atop("Increase","in" ~ italic(K[rel]))) )) +
        theme_bw() + fontTheme +
         theme(legend.justification = "center",
               axis.ticks.length.x.bottom = unit(4, "mm"),
           axis.title.x= element_blank(),
               axis.text.x = element_text(size=30, angle=90, face=1.3, vjust=0.65, hjust=0.5),
              axis.text.y = element_text(size=26),
               axis.title.y = element_text(size=32), legend.title = element_blank())
ggsave("CountFrequency_KrelFC_FP_trends_kinitFC_PauseSum.pdf", height=8)
# Fig S4D-I can be found from session: https://genome.ucsc.edu/s/rmukherjee/Figure_for_FP_treatment_Jonkers_et_al_mm39

# Fig S5 are from analysis scripts and will be too long to completely describe here.
# please refer to `Nascent RNA methods` on guertinlab (GitHub) to process raw sequencing files
# please refer to TRP `runDeseq2_TRP_12min_matched_internalControl.R`, `run_pauseworkflow_InternalTRP_matched.R`
# please refer to FP `runDeseq2_FP_25min_matched_internalControl.R`, `run_pauseworkflow_InternalFP_matched.R`

## Fig S6 - plots of genes with limits on changes in initiation rate
summmary.TRP.entireRange = prepareDf.w.Kpre_by_Krel(summmary.TRP.entireRange)
## change WithinRange to `TRUE` to get corresponding figure
pdf("TRP_outside_bounds.pdf", width= 7, height=8)
#ggplot(subset(summmary.TRP.entireRange, WithinRange == TRUE)) +
ggplot(subset(summmary.TRP.entireRange, WithinRange == FALSE)) +
        geom_point(aes(x=P_FC, y= Yind), color = "red") +
        geom_point(aes(x=B_FC, y= Yind), color = "blue") +
        geom_segment(aes(x=P_FC, xend = B_FC, y = Yind, yend = Yind), arrow = arrow(length = unit(2,"mm"), ends = "both", type = "open"), alpha = 0.15) +
        scale_x_continuous(trans=log2_trans(),
          breaks=trans_breaks('log2', function(x) 2^x),
          labels=trans_format('log2', math_format(.x))) +
        xlab(expression(log[2]~"bounds on FC("~italic(K[init])~")")) +
        ylab(expression(atop("gene indices with", "FC("~italic(K[init])~") bounds outside 0.25"))) +
        #ylab(expression(atop("gene indices with", "FC("~italic(K[init])~") bounds spanning 0.25"))) +
        theme_bw() + fontTheme
dev.off()
## Fig S7 - body percentiles vs Kpre/Krel
prob.list = c(c(1:10)*0.1)
prob.list.percent = sprintf("%0.0f%%", prob.list*100)
summmary.TRP.entireRange$P_quantile = ">99%"
summmary.TRP.entireRange$B_by_P_control = summmary.TRP.entireRange$B_base_mean / summmary.TRP.entireRange$P_base_mean
summmary.TRP.entireRange$B_by_P_quantile = ">99%"
summmary.TRP.entireRange$B_quantile = ">99%"
i = 1
prevThres = 0
for (thresh in quantile(summmary.TRP.entireRange$P_base_mean, probs = prob.list)) {
  summmary.TRP.entireRange[summmary.TRP.entireRange$P_base_mean <= thresh & summmary.TRP.entireRange$P_base_mean > prevThres,]$P_quantile = prob.list.percent[i]
  i = i + 1
  prevThres = thresh
}

df = summmary.TRP.entireRange %>% subset( WithinRange == TRUE & Z_factor > 0) %>% group_by(B_quantile) %>%
        summarise_at(c(
                       #"Kpre_min", "Kinit_min",
                       "halfLife_base_second",
                       "Z_factor",
                       #"Kinit_by_Kpre", "Kinit_by_Krel", "Kinit_by_Kpre_And_Krel",
                       #"Kpre_timesP",
                       "Krel_timesP",
                       #"Kpre_And_Krel_timesP", "Occupancy",
                       "Kinit_by_Kpre_timesP"),
                       #"Kinit_by_Krel_timesP", "Kinit_by_Kpre_And_Krel_timesP"),
                     "median") %>%
         mutate(
                #Kpre_min_str = sprintf("%.1f", Kpre_min),
               Z_factor_str = sprintf("%.1f", Z_factor),
            #    Kinit_min_str = sprintf("%.1f", Kinit_min),
            #   Kinit_by_Kpre_str = sprintf("%.1f", Kinit_by_Kpre),
            #   Kinit_by_Krel_str = sprintf("%.1f", Kinit_by_Krel),
            #   Kinit_by_Kpre_An_Krel_str = sprintf("%.1f", Kinit_by_Kpre_And_Krel),
            #   Kpre_timesP_str = sprintf("%.1f", Kpre_timesP),
              Krel_timesP_str = sprintf("%.1f", Krel_timesP),
            #   Kpre_And_Krel_timesP_str = sprintf("%.1f", Kpre_And_Krel_timesP),
            #   Occupancy_str = sprintf("%.0f", Occupancy),
               Kinit_by_Kpre_timesP_str = sprintf("%.1f", Kinit_by_Kpre_timesP),
            #   Kinit_by_Krel_timesP_str = sprintf("%.1f", Kinit_by_Krel_timesP),
            #   Kinit_by_Kpre_And_Krel_timesP_str = sprintf("%.1f", Kinit_by_Kpre_And_Krel_timesP),
                halfLife_base_sec_str = sprintf("%.1f", halfLife_base_second))
ggplot(subset(summmary.TRP.entireRange, Z_factor > 0), aes(x=B_quantile, y = Z_factor)) +
     geom_violin(scale="width") +
       geom_boxplot(width=0.35) +
         xlab("percentiles of Pause Sum (untreated)") +
         #xlab("percentiles of Body Density (untreated)") +
         #xlab(expression("percentiles of B/P ("~italic(K[rel])~"), untreated")) +
         ggtitle("mESC cells, Triptolide (repressed)") +
         geom_text(data = df, aes(
                                 x = P_quantile,
                                 #x = B_quantile,
                                 y = Z_factor,
                                 #  y = Kinit_min,
                                #y = Kpre_timesP,
                                #y = Krel_timesP,
                                  #y = Occupancy,
                                 #y = Kinit_by_Kpre_timesP,
                                #y = Kinit_by_Krel_timesP,
                                  # y = Kpre_min,
                                 # y = halfLife_base_second,
                               label = Z_factor_str
                               #label = Occupancy_str
                                #label = Kpre_timesP_str
                                #label = Krel_timesP_str
                               #label = Kinit_by_Kpre_timesP_str
                                #label = Kinit_by_Krel_timesP_str
                                #label = Kinit_min_str
                                #label = Kpre_min_str
                           #       label = halfLife_base_sec_str
                                ),
            vjust = -0.45, size = 5, color = "black") +
scale_x_discrete(labels = c("0-10","10-20","20-30","30-40", "40-50","50-60","60-70","70-80","80-90", "90-100")) +
scale_y_continuous(trans=log10_trans(),
                    guide=guide_axis_logticks(),
                    #n.breaks=6,
                    breaks=10^(seq(-2,2,1)),
                    labels=scales::label_number(accuracy = 0.01)) +
 theme_bw() + fontTheme + theme(plot.title = element_text(size=23, face=1.2,hjust=0.5),
              strip.text = element_text(size=24),
              axis.text.x = element_text(size=28, face=1.5, angle=30, vjust=0.5), #angle=30, vjust=0.9, hjust=1),
              #axis.title.x = element_blank(),
              axis.title.x  = element_text(size=27, face=1.5, hjust = 0.5),
              axis.title.y = element_text(size=32),
              legend.title = element_blank())

## Fig S8 are from analysis scripts and will be too long to completely describe here.
# please refer to `Nascent RNA methods` on guertinlab (GitHub) to process raw sequencing files
# please refer to TBP `runDeseq2_matched.R`, `run_pauseworkflow.R`
# please refer to ZNF143 `runDeseq2.R`, `run_pauseworkflow_sizeFactor.R`
# please refer to Dexamethasone - `runDeseq2_matchIt.R` and `run_pauseworkflow.R` 
# for Heatshock in Drosophila- `runDeseq2_matchIt.R` and `run_pauseworkflow.R`

 ## Fig S9 - multi datasets
summary.HS.20min.entireRange$datasetName = "HSF Activated (S2)"
summary.Dex100nm.entireRange$datasetName = "Dex, 100nm"
summary.TBP.entireRange$datasetName = "TBP dTAG"
summary.ZNF143.entireRange$datasetName = "ZNF143 dTAG"
data.df = rbind(summary.TBP.entireRange, summary.ZNF143.entireRange, summary.Dex100nm.entireRange,
                #summary.Celastrol.entireRange, summary.Vihervara_HS.entireRange,
                summary.HS.20min.entireRange)
#data.df$datasetName = factor(data.df$datasetName, levels = c("TBP dTAG", "ZNF143 dTAG", "Dex, 100nm", "HS, 30min (K562)", "HS, 20min (S2)"))
data.df$datasetName = factor(data.df$datasetName, levels = c("TBP dTAG", "ZNF143 dTAG", "Dex, 100nm",
                                                             #"Celastrol (K562)","HS, 30min (K562)",
                                                             "HSF Activated (S2)"))
labelData = c(expression(atop("TBP dTAG", "(HAP1, 2 hrs)")),
              expression(atop("ZNF143 dTAG", "(HEK293T, 30 mins)")) ,
              expression(atop("Dex 100nm", "(A549, 45 mins)")),
              # expression(atop("Celastrol", "(K562, 20 mins)")),
              #expression(atop("Heat Shock", "(K562, 30 mins)")),
              expression(atop("HSF Activated", "(S2, HS 20 mins)")))
              #expression(atop("Heatshock", "(S2, 20 mins)")))
Fig3_kinit.df = kinitFCvsPfoldchangePlotter_multidatasets(data.df, "Fig3_kinitFC_TBP_ZNF143_Dex_HS_Droso_new", labelData, Y_DOWN_arr = rep(1/5,4),  Y_UP_arr =  rep(2^3.8,4), filetype="pdf", H_JUST_LEFT = 0.8, H_JUST_RIGHT=0.2, showPercentLabels = TRUE, legend.colors = c("#36454F", "#D2B48C"), diff_UP = 19, diff_DOWN = 0.2, WIDTH=20.5, LOG_trans_break = 2^(seq(-4,8,2)), LOG_trans_format = scales::label_number(accuracy = 0.01))

# Fig S10
## U2OS cells, GR targets
## U2OS cells, with percentile
scanDataFile.entireRange = "/Users/rudradeepmukherjee/Documents/UConn Health/comparment.model.project/datasets, genome/ErinWissink_et_al_2021/ParamEstimation/paramFit_U2OS_GRtargets_percentile/calculatedRateData_U2OS.GR_targets_activated.rds"
summary.U2OS.Dex100nm = summaryFunction(readRDS(scanDataFile.entireRange))
summary.U2OS.Dex100nm$datasetName = "Dex, 100nm (U2OS, GR targets)"
krelFoldchangePlotter(summary.U2OS.Dex100nm, "U2OS_GRtargets_krel_FC_violin", "GR targets, Dex (U2OS, 45 mins)", "(Dex, 45 mins)", constraintType = "SteurerVock",  filetype="pdf", pointALPHA=0.15)
kinitFCvsPfoldchangePlotter(summary.U2OS.Dex100nm, "U2OS_GRtargets_kinitFC_vs_P_FC_times_krel", "GR targets, Dex", "(U2OS, 45 mins)",  directionKinit = "up", filetype="pdf", pointALPHA=0.1)
## C7 cells
scanDataFile.entireRange = "/Users/rudradeepmukherjee/Documents/UConn Health/comparment.model.project/datasets, genome/Arun_C7_EtOH_Dex_data/ParamEstimation/paramFit_C7_Activated/calculatedRateData_C7.Dex_1hr.rds"
summary.C7.Dex = summaryFunction(readRDS(scanDataFile.entireRange))
summary.C7.Dex$datasetName = "Dex 1hr (C7)"

krelFoldchangePlotter(summary.C7.Dex, "C7_Activated_krel_FC_violin", "Activated genes, Dex (C7, 1 hr)", "(Dex, 1 hr)", constraintType = "SteurerVock",  filetype="pdf", pointALPHA=0.15)
kinitFCvsPfoldchangePlotter(summary.C7.Dex, "C7_Activated_kinitFC_vs_P_FC_times_krel", "Activated genes, Dex", "(C7, 1 hr)",  directionKinit = "up", filetype="pdf", pointALPHA=0.1)

## Fig S11
summmary.TRP.entireRange$datasetName = "TRP 12.5mins (mESCs)"
summmary.FP.entireRange$datasetName = "FP 25mins (mESCs)"
summary.ZNF143.entireRange$datasetName = "ZNF143 dTAG (HEK293T)"
summary.TBP.entireRange$datasetName = "TBP dTAG (HAP1)"
summary.Dex100nm.entireRange$datasetName = "Dex 45mins (A549)"
summary.Celastrol.entireRange$datasetName = "Celastrol (K562)"
summary.Vihervara_HS.entireRange$datasetName = "HS 30mins (K562)"
summary.HS.20min.entireRange$datasetName = "HSF Activated (S2)"
labelData = c(expression("TRP 12.5mins (mESCs)"), expression("FP 25mins (mESCs)"), expression("ZNF143 dTAG (HEK293T)"), expression("TBP dTAG (HAP1)") , expression("Dex 45mins (A549)"),
              #expression("Celastrol (K562)"), expression("HS 30mins (K562)"),
              expression("HSF Activated (S2)"))
data.df = rbind(summmary.TRP.entireRange[,colnames(summmary.FP.entireRange)], summmary.FP.entireRange, summary.ZNF143.entireRange, summary.TBP.entireRange, summary.Dex100nm.entireRange,
                    #summary.Celastrol.entireRange, summary.Vihervara_HS.entireRange,
                  summary.HS.20min.entireRange)
data.df$datasetName = factor(data.df$datasetName, levels = labelData)

halflife.df = halfLifePlotter_indv_Kpre_forGenes_multidatasets(data.df, "Fig_Halflife_indv_Kpre_TRP_to_HS_factor_10", labelData, nROW = 2, factorMultipier = 10, useControlOnly = T, log_scale=T, HEIGHT = 9, WIDTH=13.5,  LOG_trans_break = 10^(seq(-1,3,1)), LOG_trans_format = scales::label_number(accuracy = 0.1))


# Fig S12
# plot ratio of FC(kinit) / FC(P) to figure out [A* Kpre + B * Krel] / [Kpre + Krel]
## prepare a data object with all different rations of FC(kelong) & FC(Kpre)
summmary.TRP.entireRange$ratio_FC_kinit_by_FC_P = summmary.TRP.entireRange$closestKinitFCtopoint25 /summmary.TRP.entireRange$P_FC
summmary.TRP.entireRange$Kpre_by_Krel_custom = NA
summmary.TRP.entireRange$customType = NA
custom.df = tibble(c())
#rate_changes = c(0.25, 0.5,1,1.5,2)
rate_changes = c(seq(1,2.3,0.25), 1.125)
#rate_changes = seq(1,2.3,0.125)
rate_changes = c(rate_changes, 1/rate_changes)
for (kpre_fc in rate_changes) {
  for (kelong_fc in rate_changes) {
    customTypestr = sprintf("(%0.2f,%0.2f)",kpre_fc,kelong_fc)
    summmary.TRP.entireRange$customType = customTypestr
    summmary.TRP.entireRange$kpreFC = kpre_fc
    summmary.TRP.entireRange$kelongFC = kelong_fc
    summmary.TRP.entireRange$P_FC_limit = kpre_fc *  summmary.TRP.entireRange$P_FC
        summmary.TRP.entireRange$B_FC_limit = kelong_fc *  summmary.TRP.entireRange$B_FC
    ## set closest FC(kinit) to 0.25
     summmary.TRP.entireRange$closestKinitFCtopoint25 = 0.25
     summmary.TRP.entireRange$WithinRange = T
  for (i in c(1:nrow(summmary.TRP.entireRange))) {
   oneEnd = summmary.TRP.entireRange[i,]$P_FC_limit
   otherEnd = summmary.TRP.entireRange[i,]$B_FC_limit
          leftEnd = min(oneEnd, otherEnd)
          rightEnd = max(oneEnd, otherEnd)
        #  print(sprintf("left: %0.1f, right: %0.1f, index: %d", leftEnd, rightEnd, i))
          if (leftEnd > 0.25 || rightEnd < 0.25) {
                   summmary.TRP.entireRange[i,]$WithinRange = FALSE
                 if (leftEnd > 0.25) {
                  summmary.TRP.entireRange[i, ]$closestKinitFCtopoint25 = leftEnd + (rightEnd - leftEnd) * 0.05
                 }
                 if (rightEnd < 0.25) {
                  summmary.TRP.entireRange[i, ]$closestKinitFCtopoint25 = rightEnd - (rightEnd - leftEnd) * 0.05
                 }
          }
          if (leftEnd == 0.25 || rightEnd == 0.25) {
            summmary.TRP.entireRange[i,]$WithinRange = FALSE
            summmary.TRP.entireRange[i,]$boundEqualspoint25 = TRUE
           if (leftEnd == 0.25) {
                  summmary.TRP.entireRange[i, ]$closestKinitFCtopoint25 = leftEnd + (rightEnd - leftEnd) * 0.05
                 }
                 if (rightEnd == 0.25) {
                  summmary.TRP.entireRange[i, ]$closestKinitFCtopoint25 = rightEnd - (rightEnd - leftEnd) * 0.05
                 }
             }
            }

    ## now calculate ratios
    summmary.TRP.entireRange$Kpre_by_Krel_custom = (kelong_fc * summmary.TRP.entireRange$B_FC - summmary.TRP.entireRange$closestKinitFCtopoint25) / (summmary.TRP.entireRange$closestKinitFCtopoint25 - kpre_fc * summmary.TRP.entireRange$P_FC)
    #summmary.TRP.entireRange$Kpre_by_Krel_custom = (kelong_fc * summmary.TRP.entireRange$B_FC - 0.25) / (0.25 - kpre_fc * summmary.TRP.entireRange$P_FC)
    custom.df= rbind(custom.df, summmary.TRP.entireRange)
  }
}

custom.df$label_kpre_FC = sprintf("%0.2f", custom.df$kpreFC)
custom.df$label_kelong_FC = sprintf("%0.2f", custom.df$kelongFC)
custom.df$krel_FC = custom.df$kelongFC * custom.df$B_FC/custom.df$P_FC
saveRDS(custom.df, "TRP.custom.kpre.krel.FCs.UseFC_Kpre_Kelong_toSetClosesKinitFCto_point25.rds")

freq_all_df = custom.df %>% group_by(customType, kpreFC, kelongFC) %>% summarise(count_exceeding_1  = sum(Kpre_by_Krel_custom > 1),
                                                                            count_betwn_0_1 = sum(Kpre_by_Krel_custom > 0 & Kpre_by_Krel_custom <=1),
                                                                           count_negative = sum(Kpre_by_Krel_custom < 0),
                                                           total_count = n()) %>%
                                              mutate(percentGreat1 = 100 * count_exceeding_1 / total_count,
                                                     percentBetwn_0_1 = 100 * count_betwn_0_1 / total_count,
                                                     percentNegative = 100 * count_negative / total_count)
median_kinit_FC =  custom.df %>% group_by(customType, kpreFC, kelongFC) %>% summarise_at(c("closestKinitFCtopoint25"), c("median"))
freq_all_df$label_kpre_FC = sprintf("%0.2f", freq_all_df$kpreFC)
freq_all_df$label_kelong_FC = sprintf("%0.2f", freq_all_df$kelongFC)
median_kinit_FC$label_kpre_FC = sprintf("%0.2f", median_kinit_FC$kpreFC)
median_kinit_FC$label_kelong_FC = sprintf("%0.2f", median_kinit_FC$kelongFC)
## plot the figure
ggplot(freq_all_df, aes(x = label_kpre_FC, y = label_kelong_FC, fill=percentGreat1)) +
       geom_tile(color = "white", size = 1) +
       scale_fill_gradient2(
            low = "#FEE090",   # Yellow for low values
            mid = "#F7F7F7",   # A neutral white for the midpoint
            high = "#A50026",  # Dark red for high values
            midpoint = 50,     # Set the midpoint of the color scale
            name = "% of genes" # Set the title of the legend
       ) +
  labs(
    title = expression("TRP repressed genes, % of genes with" ~ frac(italic(K[pre]), italic(K[rel]^"control")) ~">1"),
    x = expression("Fold change in" ~ italic(K[pre])),
    y = expression("Fold change in" ~ italic(K[elong]))
  ) +
theme_minimal(base_size = 12) + fontTheme +
theme(
    plot.title = element_text(size=22, hjust = 0.5, face = 1.2), # Center the plot title
    axis.title.x = element_text(size=25),
    axis.title.y = element_text(size=25),
    axis.text.x = element_text(angle = 45, hjust = 1),     # Rotate x-axis labels for readability
    panel.grid.major = element_blank(),                     # Remove major grid lines
    panel.grid.minor = element_blank(),                     # Remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1) # Add a border around the plot
  )
ggsave("TRP_heatmap_Kpre_by_Krel_withGenes_KinitFC_span_point25.pdf", height=8, width=9)

# Fig S13
df = summarizeKelongvsKrelKpreKinit(subset(summmary.TRP.entireRange, WithinRange == T), filename = "TRP_vary_kelong", kelongVals = c(1000,2600, 3400,5000)/60, VJUST=-0.5, HJUST=0.5, WIDTH=8, LOG_trans_break = 10^(seq(-2,5,1)),LOG_trans_format = scales::label_number(accuracy = 0.01))
# Fig S14
kinitFCvals.list = c(0.05, 0.1,0.2, 0.25, 0.3, 0.35, 0.45, 0.65, 0.85)
df = summarizeFC_KinitvsKpreByKrel_Kpre_Kinit(summmary.TRP.entireRange, filename = "TRP_vary_KinitFC", kinitFCvals = kinitFCvals.list, WIDTH=12, HEIGHT=9, LOG_trans_break = 10^(seq(-2,6,1)), LOG_trans_format = scales::label_number(accuracy = 0.1))
# Fig S15
rootDir = "/Users/rudradeepmukherjee/Documents/UConn Health/comparment.model.project/datasets, genome/Jonkers_2014_Trp_FP_paper"
controlDens = read.table(sprintf("%s/ParamEstimation_FP_TRP/TRP_control_genes_P_B_bodyOnly4kb.txt", rootDir),  sep="\t", header=T)
treatDens = read.table(sprintf("%s/ParamEstimation_FP_TRP/TRP_treatment_genes_P_B_bodyOnly4kb.txt", rootDir),  sep="\t", header=T)
df = summarizeP_BscalevsKpreByKrel_Kpre_Kinit(controlDens, treatDens, filename = "TRP_vary_P_B_scale", prob.list = c(0.48,0.7,0.90,0.95,0.998), WIDTH=10, LOG_trans_break = 10^(seq(-2,6,1)), LOG_trans_format = scales::label_number(accuracy = 0.01))
## Fig S16
# see Fig S7 how to prepare P_quantile and B_quantile ranges
# for median related to P_quantile
df = summmary.TRP.entireRange %>% subset(WithinRange == TRUE & Z_factor > 0) %>% group_by(P_quantile) %>%
        summarise_at(c("Kpre_min", "Kinit_min", "halfLife_base_second",
                       "Kinit_by_Kpre", "Kinit_by_Krel", "Kinit_by_Kpre_And_Krel",
                       "Kpre_timesP", "Krel_timesP", "Kpre_And_Krel_timesP", "Occupancy",
                       "Kinit_by_Kpre_timesP", "Kinit_by_Krel_timesP", "Kinit_by_Kpre_And_Krel_timesP"), "median") %>%
         mutate(Kpre_min_str = sprintf("%.1f", Kpre_min),
                Kinit_min_str = sprintf("%.1f", Kinit_min),
               Kinit_by_Kpre_str = sprintf("%.1f", Kinit_by_Kpre),
               Kinit_by_Krel_str = sprintf("%.1f", Kinit_by_Krel),
               Kinit_by_Kpre_And_Krel_str = sprintf("%.1f", Kinit_by_Kpre_And_Krel),
               Kpre_timesP_str = sprintf("%.1f", Kpre_timesP),
               Krel_timesP_str = sprintf("%.1f", Krel_timesP),
               Kpre_And_Krel_timesP_str = sprintf("%.1f", Kpre_And_Krel_timesP),
               Occupancy_str = sprintf("%.0f", Occupancy),
               Kinit_by_Kpre_timesP_str = sprintf("%.1f", Kinit_by_Kpre_timesP),
               Kinit_by_Krel_timesP_str = sprintf("%.1f", Kinit_by_Krel_timesP),
               Kinit_by_Kpre_And_Krel_timesP_str = sprintf("%.1f", Kinit_by_Kpre_And_Krel_timesP),
                halfLife_base_sec_str = sprintf("%.1f", halfLife_base_second))
# for median related to B_quantile
         df = summmary.TRP.entireRange %>% subset( WithinRange == TRUE & Z_factor > 0) %>% group_by(B_quantile) %>%
        summarise_at(c(
                       #"Kpre_min", "Kinit_min",
                       "halfLife_base_second",
                       "Z_factor",
                       #"Kinit_by_Kpre", "Kinit_by_Krel", "Kinit_by_Kpre_And_Krel",
                       #"Kpre_timesP",
                       "Krel_timesP",
                       #"Kpre_And_Krel_timesP", "Occupancy",
                       "Kinit_by_Kpre_timesP"),
                       #"Kinit_by_Krel_timesP", "Kinit_by_Kpre_And_Krel_timesP"),
                     "median") %>%
         mutate(
                #Kpre_min_str = sprintf("%.1f", Kpre_min),
               Z_factor_str = sprintf("%.1f", Z_factor),
            #    Kinit_min_str = sprintf("%.1f", Kinit_min),
            #   Kinit_by_Kpre_str = sprintf("%.1f", Kinit_by_Kpre),
            #   Kinit_by_Krel_str = sprintf("%.1f", Kinit_by_Krel),
            #   Kinit_by_Kpre_An_Krel_str = sprintf("%.1f", Kinit_by_Kpre_And_Krel),
            #   Kpre_timesP_str = sprintf("%.1f", Kpre_timesP),
              Krel_timesP_str = sprintf("%.1f", Krel_timesP),
            #   Kpre_And_Krel_timesP_str = sprintf("%.1f", Kpre_And_Krel_timesP),
            #   Occupancy_str = sprintf("%.0f", Occupancy),
               Kinit_by_Kpre_timesP_str = sprintf("%.1f", Kinit_by_Kpre_timesP),
            #   Kinit_by_Krel_timesP_str = sprintf("%.1f", Kinit_by_Krel_timesP),
            #   Kinit_by_Kpre_And_Krel_timesP_str = sprintf("%.1f", Kinit_by_Kpre_And_Krel_timesP),
                halfLife_base_sec_str = sprintf("%.1f", halfLife_base_second))
# toggle betweel P_quantile and B_quantile
ggplot(subset(summmary.TRP.entireRange, WithinRange == TRUE & Z_factor > 0), aes(x=P_quantile, y = Kinit_by_Kpre_timesP)) +
#ggplot(subset(summmary.TRP.entireRange, WithinRange == TRUE & Z_factor > 0), aes(x=B_quantile, y = Kinit_by_Kpre_timesP)) +
   geom_violin(scale="width") +
       geom_boxplot(width=0.35) +
         xlab("percentiles of Pause Sum (untreated)") +
         #xlab("percentiles of Body Density (untreated)") +
         ggtitle("mESC cells, Triptolide (repressed)") +
         geom_text(data = df, aes(
                                 x = P_quantile,
                                 #x = B_quantile,
                                 y = Kinit_by_Kpre_timesP,
                                 label = Kinit_by_Kpre_timesP_str),
                               vjust = -0.45, size = 5, color = "black") +
scale_x_discrete(labels = c("0-10","10-20","20-30","30-40", "40-50","50-60","60-70","70-80","80-90", "90-100")) +
 scale_y_continuous(trans=log10_trans(),
                    guide=guide_axis_logticks(),
                    #n.breaks=6,
                    breaks=10^(seq(-2,2,1)),
                    labels=scales::label_number(accuracy = 0.01)) +
   ylab(expression(frac(italic(K[init]), italic(K[pre]) ~ "x P"))) +
    theme_bw() + fontTheme + theme(plot.title = element_text(size=23, face=1.2,hjust=0.5),
              strip.text = element_text(size=24),
              axis.text.x = element_text(size=28, face=1.5, angle=30, vjust=0.5), #angle=30, vjust=0.9, hjust=1),
              #axis.title.x = element_blank(),
              axis.title.x  = element_text(size=27, face=1.5, hjust = 0.5),
              axis.title.y = element_text(size=32),
              legend.title = element_blank())
ggsave("TRP_quantiles_P_vs_kinit_by_KprexP.pdf", width=12)
ggsave("TRP_quantiles_B_vs_kinit_by_KprexP.pdf", width=12)









