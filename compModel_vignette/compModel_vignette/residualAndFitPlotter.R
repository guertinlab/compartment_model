library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(scales)
options(error = traceback)

fontTheme = theme(
  axis.title.x = element_text(size = 20),
  #axis.text.x = element_text(size = 4),
  #axis.text.x = element_text(size = 11),
  #axis.text.x = element_text(size=28, angle = 30, hjust=1),
  axis.text.x = element_text(size=27),
  axis.text.y = element_text(size = 25),
  axis.title.y = element_text(size = 28),
  legend.title = element_text(size=25),
  legend.text = element_text(size=24),
  panel.border = element_rect(colour="black", linewidth=2),
  panel.grid.major = element_blank(),
 panel.grid.minor = element_blank(),
 axis.ticks.length = unit(2, "mm")
 #axis.line = element_line(colour = "black"),
 #axis.line.x = element_line(colour = "black"),
 #axis.line.y = element_line(colour = "black"),
  #legend.position="top"
 #axis.line.x.top = element_line(colour = "black"),
 #axis.line.y.right = element_line(colour = "black")
)

generalized_tanh <- function(x, S1, S2, B) {
  A <- (S2 - S1) / 2
  C <- (S2 + S1) / 2
  return (A * tanh(B * x) + C)
}

percentilePlotter <- function(datasetName, originalRankedPauseSums, percentileVal, treatment = "(control)", plotTitle=NA, logYPlotter = "No", logXPlotter = "No", showLegend = FALSE, outputDir=NULL) {
   originalRankedPauseSums = originalRankedPauseSums[rank(order(originalRankedPauseSums$pause_sum)), ]
   N_row = nrow(originalRankedPauseSums)
   xs <- 1:(N_row)
   #ysForSaturation <- originalRankedPauseSums[xsForSaturation, 'pause_sum']
   data_fitting <- data.frame(x = xs, y = originalRankedPauseSums$pause_sum, type = "Original Data for Fitting")
   plot_data <- data_fitting
   maxSaturation = originalRankedPauseSums[round(N_row * percentileVal/100),]$pause_sum
   p <- ggplot(plot_data, aes(x = x, y = y, color = type)) +
    geom_line(size=1.8) + 
    ggtitle(plotTitle) 
   p <- p + geom_hline(yintercept = maxSaturation, linetype = "dashed", size=1.3) + #linetype = sprintf("Max saturation (fitted curve) = %0.2f", maxSaturation))
            annotate("text", x=N_row / 2,y=maxSaturation,label = sprintf("Percentile (%0.1f) = %0.0f", percentileVal, maxSaturation), vjust = -1, color="black", family = "mono", size=8)
   xLabel = paste("Rank of empirical pause sums", treatment)
   yLabel = paste("Empirical pause sums", treatment)
   if (logYPlotter == "Yes") {
     p <- p + scale_y_continuous(trans=log10_trans(),
                breaks=trans_breaks('log10', function(x) 10^x),
                labels=trans_format('log10', math_format(.x)))
      p <- p + labs(y=bquote(log[10](.(yLabel))))
    }
   else {
      p <- p + ylab(yLabel)
   }
   if (logXPlotter == "Yes") {
      p <- p + scale_x_continuous(trans=log10_trans(),
                breaks=trans_breaks('log10', function(x) 10^x),
                labels=trans_format('log10', math_format(.x))) 
 p <- p + labs(x=bquote(log[10](.(xLabel))))
} else {
   p <- p + xlab(xLabel)
}
p <- p + theme_bw() + fontTheme
if (!is.na(plotTitle)) {
 p <- p + labs(title = plotTitle) +
         theme(plot.title = element_text(size = 28, face=1.5, hjust=0.5))
}  

p <- p + theme(
    text = element_text(size = 22),
     axis.text = element_text(size=23),
      axis.title.y = element_text(size=30),
      axis.title.x = element_text(size=28),
    #legend.title = element_blank()
  )
Width = 15
if (showLegend == FALSE) {
        p <- p + theme(legend.position = "None")
Width = 10
}
   ggsave(filename = sprintf('%s/%s_threshPercentile_%0.2f_%s_logX_%s_logY.pdf', outputDir, datasetName, percentileVal, logXPlotter, logYPlotter), plot = p, width = Width, height = 8)
          
} 

percentilePlotter_multiPercentile <- function(datasetName, originalRankedPauseSums, percentileVals, treatment = "(control)", plotTitle=NA, logYPlotter = "No", logXPlotter = "No", showLegend = FALSE, outputDir=NULL) {
   originalRankedPauseSums = originalRankedPauseSums[rank(order(originalRankedPauseSums$pause_sum)), ]
   N_row = nrow(originalRankedPauseSums)
   xs <- 1:(N_row)
   #ysForSaturation <- originalRankedPauseSums[xsForSaturation, 'pause_sum']
   data_fitting <- data.frame(x = xs, y = originalRankedPauseSums$pause_sum, type = "Original Data for Fitting")
   plot_data <- data_fitting
   p <- ggplot(plot_data, aes(x = x, y = y, color = type)) +
    geom_line(size=1.8) + 
    ggtitle(plotTitle) 
   for (percentileVal in percentileVals) {
           maxSaturation = quantile(originalRankedPauseSums$pause_sum, probs =c( percentileVal/100))[[1]]
      p <- p + geom_hline(yintercept = maxSaturation, linetype = "dashed", size=1.3)  #linetype = sprintf("Max saturation (fitted curve) = %0.2f", maxSaturation))
           p <- p + annotate("text", x=N_row / 2,y=maxSaturation,label = sprintf("Percentile (%0.1f) = %0.0f", percentileVal, maxSaturation), vjust = -1, color="black", family = "mono", size=8)
   }
           xLabel = paste("Rank of empirical pause sums", treatment)
   yLabel = paste("Empirical pause sums", treatment)
   if (logYPlotter == "Yes") {
     p <- p + scale_y_continuous(trans=log10_trans(),
                breaks=trans_breaks('log10', function(x) 10^x),
                labels=trans_format('log10', math_format(.x)))
      p <- p + labs(y=bquote(log[10](.(yLabel))))
    }
   else {
      p <- p + ylab(yLabel)
   }
   if (logXPlotter == "Yes") {
      p <- p + scale_x_continuous(trans=log10_trans(),
                breaks=trans_breaks('log10', function(x) 10^x),
                labels=trans_format('log10', math_format(.x))) 
 p <- p + labs(x=bquote(log[10](.(xLabel))))
} else {
   p <- p + xlab(xLabel)
}
p <- p + theme_bw() + fontTheme
if (!is.na(plotTitle)) {
 p <- p + labs(title = plotTitle) +
         theme(plot.title = element_text(size = 28, face=1.5, hjust=0.5))
}  

p <- p + theme(
    text = element_text(size = 22),
     axis.text = element_text(size=23),
      axis.title.y = element_text(size=30),
      axis.title.x = element_text(size=28),
    #legend.title = element_blank()
  )
Width = 15
if (showLegend == FALSE) {
        p <- p + theme(legend.position = "None")
Width = 10
}
   ggsave(filename = sprintf('%s/%s_multiPercentiles_%s_logX_%s_logY.pdf', outputDir, datasetName, logXPlotter, logYPlotter), plot = p, width = Width, height = 8)
          
} 
saturationRankPlotter <- function(datasetName, paramFileloc, outputDir, plotTitle = "S2 cells", typeofFit="LOG_HYPER_TANG", WIDTH=7) { 
        paramInfo <- read.table(paramFileloc, sep="\t", header=T)

        p <- ggplot(paramInfo, aes(x = paramInfo['rankIndex'][[1]], y = paramInfo['S2'][[1]])) +
               geom_point() + 
                xlab('rank fitted') + 
               ylab(expression(log[10]~"saturation")) +
               ggtitle(plotTitle) + 
                theme_bw() + fontTheme + theme(
    text = element_text(size = 22),
     axis.text = element_text(size=23),
      axis.title.y = element_text(size=30),
      axis.title.x = element_text(size=28),
     plot.title = element_text(size = 25, face=1.5, hjust=0.5))
      ggsave(filename = sprintf('%s/%s_%s_SaturationVsRank.pdf', outputDir, datasetName, typeofFit), plot = p, width = WIDTH, height = 8)

}
## R version of functions from `residualAndFitPlotter.py`
fitPlotter_multiplePercentiles_logPlotter <- function(datasetName, paramFileloc, residualfileloc, originalRankedPauseSums, percentileVals = c(90,95,100), typeofFit="LOG_HYPER_TANG", is75percentile="No" ,startIndex=1, endIndex=3, logYPlotter="Yes", logXPlotter="No", plotMaxSaturation="No", plotTitle=NA,plotRestOfFunction="Yes", outputDir=NULL, showLegend=TRUE, plotRestOfData="Yes", treatment = "(control)", V_JUST = c(1.2,1.2,1.2), H_JUST = c(0.5,0.5,0.5), color.list = c( "#EBD3F8", "#AD49E1", "#7A1CAC"), legendPosition = c(0.8,0.2), WIDTH=10) {
   originalRankedPauseSums = originalRankedPauseSums[rank(order(originalRankedPauseSums$pause_sum)),]
     N_row = nrow(originalRankedPauseSums)
     ecdf.pause_sum = ecdf(originalRankedPauseSums$pause_sum)
    paramInfo <- read.csv(paramFileloc, sep="\t", header=T)
    fitDf <- read.csv(residualfileloc, sep="\t", header=T)
    if (typeofFit == "RATIONAL_FUNC") {
        paramInfo <- paramInfo %>% filter(A / P > 0)
    }
    levelList = c("Original Data")
    ## extract saturation values for percentiles
    fittedData.df = as.data.frame(c())
    saturationVal.df = as.data.frame(c())
    for (percentileval in percentileVals) {
      indexToConsider = round(N_row * percentileval / 100.0) 
      parameterSet = subset(paramInfo, rankIndex == indexToConsider)
      if (typeofFit == "LOG_HYPER_TANG") {
      S1 <- parameterSet["S1"][[1]]
      S2 <- parameterSet["S2"][[1]]
      B <- parameterSet["B"][[1]]
      maxSaturation = S2
      percentileMaxSaturation = ecdf.pause_sum(maxSaturation) * 100.0
      xsForSaturation = 1:(indexToConsider)
      fittedSaturation <- generalized_tanh(xsForSaturation, S1, S2, B)
      #print(length(fittedSaturation))
      #typestring =  sprintf("Fitted Data, percentile %0.1f", percentileval)
      typestring =  sprintf("%0.0f percentile", percentileval)
      #saturationTypeString =  sprintf("Saturation value = %0.0f (%0.0f percentile, fit at %0.1f percentile)", maxSaturation, percentileMaxSaturation, percentileval)
      saturationTypeString =  sprintf("Saturation = %0.0f (%0.0f percentile)", maxSaturation, percentileMaxSaturation)

      fittedData.df = rbind(fittedData.df, data.frame(x = xsForSaturation, y = fittedSaturation, type = typestring))
      saturationVal.df = rbind(saturationVal.df, data.frame(percentile = percentileval, saturationVal = maxSaturation, 
                          rankIndex = indexToConsider, type = saturationTypeString))      
         levelList = c(levelList,typestring) 
      }
    }
    originData <- data.frame(x = (1:N_row), y = log(originalRankedPauseSums$pause_sum, base=10), type = "Original Data")
    #print(ncol(originData))
    #print(ncol(fittedData.df))
    plot_data <- rbind(fittedData.df, originData)
    plot_data$type = factor(plot_data$type, levels = levelList)
    p <- ggplot(plot_data, aes(x = x, y = y, color = type)) +
        geom_line(linewidth=1.8) + 
        scale_color_manual(values = c("black", color.list)) 
        #ggtitle(plotTitle) 
   #color.list = c( "#EBD3F8", "#AD49E1", "#7A1CAC")  
   ind = 1
   if (plotMaxSaturation == "Yes") {
           for (percentileval in percentileVals) {
      subset.saturationVal.df = subset(saturationVal.df, percentile == percentileval)
      p <- p + geom_hline(yintercept = subset.saturationVal.df$saturationVal , linetype = "dashed", size=1.3, color = color.list[ind]) + #linetype = sprintf("Max saturation (fitted curve) = %0.2f", maxSaturation))
            annotate("text", x=N_row / 2,y= subset.saturationVal.df$saturationVal,label = subset.saturationVal.df$type, hjust = H_JUST[ind],vjust = V_JUST[ind], family = "mono", size=8, color =  color.list[ind], fontface="bold") 
      ind = ind + 1 
   } 
} 
    
 xLabel = paste("Rank of pause sums", treatment)
yLabel = paste("Pause sums", treatment)
 p <- p + labs(y=bquote(log[10](.(yLabel))))
#if (logYPlotter == "Yes") {
#  p <- p + scale_y_continuous(trans=log10_trans(),
#                breaks=trans_breaks('log10', function(x) 10^x),
#                labels=trans_format('log10', math_format(.x)))
# p <- p + labs(y=bquote(log[10](.(yLabel))))
#}
#else {
#   p <- p + ylab(yLabel)
#}
if (logXPlotter == "Yes") {
  p <- p + scale_x_continuous(trans=log10_trans(),
                breaks=trans_breaks('log10', function(x) 10^x),
                labels=trans_format('log10', math_format(.x))) 
 p <- p + labs(x=bquote(log[10](.(xLabel))))
} else {
   p <- p + xlab(xLabel)
}
p <- p + theme_bw() + fontTheme
if (!is.na(plotTitle)) {
 p <- p + labs(title = plotTitle) +
         theme(plot.title = element_text(size = 28, face=1.5, hjust=0.5))
}  

p <- p + theme(
    text = element_text(size = 22),
     axis.text = element_text(size=23),
      axis.title.y = element_text(size=30),
      axis.title.x = element_text(size=28),
    legend.title = element_blank(),
    legend.text =  element_text(size=30),
    legend.position = legendPosition
  )
#Width = 20
if (showLegend == FALSE) {
        p <- p + theme(legend.position = "None")
#Width = 14
}
   ggsave(filename = sprintf('%s/%s_%s_multiPercentiles_logY.pdf', outputDir, datasetName, typeofFit), plot = p, width = WIDTH, height = 8)
     
}


## R version of functions from `residualAndFitPlotter.py`
fitPlotter_multiplePercentiles <- function(datasetName, paramFileloc, originalRankedPauseSums, percentileVals = c(90,95,100), typeofFit="LOG_HYPER_TANG", is75percentile="No" ,startIndex=1, endIndex=3, logYPlotter="Yes", logXPlotter="No", plotMaxSaturation="No", plotTitle=NA,plotRestOfFunction="Yes", outputDir=NULL, showLegend=TRUE, plotRestOfData="Yes", treatment = "(control)", V_JUST = c(1.2,1.2,1.2), H_JUST = c(0.5,0.5,0.5), color.list = c( "#EBD3F8", "#AD49E1", "#7A1CAC"), legendPosition = c(0.8,0.2), WIDTH=10) {
   originalRankedPauseSums = originalRankedPauseSums[rank(order(originalRankedPauseSums$pause_sum)),]
     N_row = nrow(originalRankedPauseSums)
     ecdf.pause_sum = ecdf(originalRankedPauseSums$pause_sum)
    paramInfo <- read.csv(paramFileloc, sep="\t", header=T)
    if (typeofFit == "RATIONAL_FUNC") {
        paramInfo <- paramInfo %>% filter(A / P > 0)
    }
    levelList = c("Original Data")
    ## extract saturation values for percentiles
    fittedData.df = as.data.frame(c())
    saturationVal.df = as.data.frame(c())
    for (percentileval in percentileVals) {
      indexToConsider = round(N_row * percentileval / 100.0) 
      parameterSet = subset(paramInfo, rankIndex == indexToConsider)
      if (typeofFit == "LOG_HYPER_TANG") {
      S1 <- parameterSet["S1"][[1]]
      S2 <- parameterSet["S2"][[1]]
      B <- parameterSet["B"][[1]]
      maxSaturation = 10^S2
      percentileMaxSaturation = ecdf.pause_sum(maxSaturation) * 100.0
      xsForSaturation = 1:(indexToConsider)
      fittedSaturation <- 10^(generalized_tanh(xsForSaturation, S1, S2, B))
      #print(length(fittedSaturation))
      #typestring =  sprintf("Fitted Data, percentile %0.1f", percentileval)
      typestring =  sprintf("%0.0f percentile", percentileval)
      #saturationTypeString =  sprintf("Saturation value = %0.0f (%0.0f percentile, fit at %0.1f percentile)", maxSaturation, percentileMaxSaturation, percentileval)
      saturationTypeString =  sprintf("Saturation = %0.0f (%0.0f percentile)", maxSaturation, percentileMaxSaturation)

      fittedData.df = rbind(fittedData.df, data.frame(x = xsForSaturation, y = fittedSaturation, type = typestring))
      saturationVal.df = rbind(saturationVal.df, data.frame(percentile = percentileval, saturationVal = maxSaturation, 
                          rankIndex = indexToConsider, type = saturationTypeString))      
         levelList = c(levelList,typestring) 
      }
    }
    originData <- data.frame(x = (1:N_row), y = originalRankedPauseSums$pause_sum, type = "Original Data")
    #print(ncol(originData))
    #print(ncol(fittedData.df))
    plot_data <- rbind(fittedData.df, originData)
    plot_data$type = factor(plot_data$type, levels = levelList)
    p <- ggplot(plot_data, aes(x = x, y = y, color = type)) +
        geom_line(linewidth=1.8) + 
        scale_color_manual(values = c("black", color.list)) 
        #ggtitle(plotTitle) 
   #color.list = c( "#EBD3F8", "#AD49E1", "#7A1CAC")  
   ind = 1
   for (percentileval in percentileVals) {
      subset.saturationVal.df = subset(saturationVal.df, percentile == percentileval)
      p <- p + geom_hline(yintercept = subset.saturationVal.df$saturationVal , linetype = "dashed", size=1.3, color = color.list[ind]) + #linetype = sprintf("Max saturation (fitted curve) = %0.2f", maxSaturation))
            annotate("text", x=N_row / 2,y= subset.saturationVal.df$saturationVal,label = subset.saturationVal.df$type, hjust = H_JUST[ind],vjust = V_JUST[ind], family = "mono", size=8, color =  color.list[ind], fontface="bold") 
      ind = ind + 1 
   } 
xLabel = paste("Rank of pause sums", treatment)
yLabel = paste("Pause sums", treatment)
if (logYPlotter == "Yes") {
  p <- p + scale_y_continuous(trans=log10_trans(),
                breaks=trans_breaks('log10', function(x) 10^x),
                labels=trans_format('log10', math_format(.x)))
 p <- p + labs(y=bquote(log[10](.(yLabel))))
}
else {
   p <- p + ylab(yLabel)
}
if (logXPlotter == "Yes") {
  p <- p + scale_x_continuous(trans=log10_trans(),
                breaks=trans_breaks('log10', function(x) 10^x),
                labels=trans_format('log10', math_format(.x))) 
 p <- p + labs(x=bquote(log[10](.(xLabel))))
} else {
   p <- p + xlab(xLabel)
}
p <- p + theme_bw() + fontTheme
if (!is.na(plotTitle)) {
 p <- p + labs(title = plotTitle) +
         theme(plot.title = element_text(size = 28, face=1.5, hjust=0.5))
}  

p <- p + theme(
    text = element_text(size = 22),
     axis.text = element_text(size=23),
      axis.title.y = element_text(size=30),
      axis.title.x = element_text(size=28),
    legend.title = element_blank(),
    legend.text =  element_text(size=30),
    legend.position = legendPosition
  )
#Width = 20
if (showLegend == FALSE) {
        p <- p + theme(legend.position = "None")
#Width = 14
}
   ggsave(filename = sprintf('%s/%s_%s_multiPercentiles_%s_showLegend_%s_logX_%s_logY_Rplotter.pdf', outputDir, datasetName, typeofFit, showLegend, logXPlotter, logYPlotter), plot = p, width = WIDTH, height = 8)
     
}

## R version of functions from `residualAndFitPlotter.py`
fitPlotter <- function(datasetName, paramFileloc, residualfileloc, originalRankedPauseSums, typeofFit="GeneralLogistic", is75percentile="No" ,startIndex=1, endIndex=3, logYPlotter="Yes", logXPlotter="No", plotMaxSaturation="No", plotTitle=NA,plotRestOfFunction="Yes", outputDir=NULL, showLegend="False",plotRestOfData="Yes", treatment = "(control)") {
    N_row = nrow(originalRankedPauseSums) 
    paramInfo <- read.csv(paramFileloc, sep="\t", header=T)
    fitDf <- read.csv(residualfileloc, sep="\t", header=T)
    if (typeofFit == "RATIONAL_FUNC") {
        paramInfo <- paramInfo %>% filter(A / P > 0)
    }
    
    for (i in seq_len(nrow(paramInfo))) {
        #row <- paramInfo[i,]
        indexThresh <- fitDf[i, 'rankIndex']
        
        if (typeofFit %in% c("QUADRATIC", "LOG_QUADRATIC", "EXPONENTIAL")) {
            mean_residual <- fitDf[i, 'residualSaturation'] / (nrow(originalRankedPauseSums) - indexThresh + 1)
        } else {
            mean_residual <- fitDf[i, 'residualSaturation'] / (indexThresh - 1)
        }
        
        if (indexThresh >= startIndex && indexThresh <= endIndex && !is.nan(mean_residual)) {
            xsForSaturation <- 1:(indexThresh)
            ysForSaturation <- originalRankedPauseSums[xsForSaturation, 'pause_sum']
            restofXs <- indexThresh:nrow(originalRankedPauseSums)
            restofYs <- originalRankedPauseSums[restofXs, 'pause_sum']
           maxP <- originalRankedPauseSums[indexThresh, 'pause_sum'] 
            if (typeofFit == "LOG_HYPER_TANG") {
                    row = subset(paramInfo, rankIndex == indexThresh)
                    S1 <- row["S1"][[1]]
                    S2 <- row["S2"][[1]]
                    B <- row["B"][[1]]
                 fittedSaturation <- 10^(generalized_tanh(xsForSaturation, S1, S2, B))
                 restOfFuncEval <- 10^(generalized_tanh(restofXs, S1, S2, B))
                 maxSaturation <- 10^S2
           }
            # Add other fit types here
           data_fitting <- data.frame(x = xsForSaturation, y = ysForSaturation, type = "Original Data for Fitting")
       data_fitted <- data.frame(x = xsForSaturation, y = fittedSaturation, type = "Fitted Data" ) #sprintf("Fitted Data - %s, up to rank %d", typeofFit, indexThresh))
       data_rest <- data.frame(x = restofXs, y = restofYs, type = "Original Data (not fitted)")#type = "Rest of original data (not fitted)")
       data_func_eval <- data.frame(x = restofXs, y = restOfFuncEval, type = sprintf("Function evaluation (on data excluded for fitting)"))#type = sprintf("Rest of function evaluation (%s)", typeofFit))
 plot_data <- rbind(data_fitting, data_fitted)
if (plotRestOfData == "Yes") {
  plot_data <- rbind(plot_data, data_rest)
}
if (plotRestOfFunction == "Yes") {
  plot_data <- rbind(plot_data, data_func_eval)
}

# Plot using ggplot2
p <- ggplot(plot_data, aes(x = x, y = y, color = type)) +
  geom_line(size=1.8) + 
 ggtitle(plotTitle) 
# labs(
 #   x = 'Rank of pause sum in original data',
 #   y = 'Pause sums',
 #   title = sprintf("Mean residual: %0.3f, max Pause Sum: %0.2f", mean_residual, maxP)
 # ) 
 if (plotMaxSaturation == "Yes") {
  p <- p + geom_hline(yintercept = maxSaturation, linetype = "dashed", size=1.3) + #linetype = sprintf("Max saturation (fitted curve) = %0.2f", maxSaturation))
            annotate("text", x=N_row / 2,y=maxSaturation,label = sprintf("Saturation value (fitted curve) = %0.2f", maxSaturation), vjust = -1, color="black", family = "mono", size=8)
          #+
 #         scale_linetype_manual(name="Reference", values = c(sprintf("Max saturation (fitted curve) = %0.2f", maxSaturation) = "dashed"))
}
xLabel = paste("Rank of empirical pause sums", treatment)
yLabel = paste("Empirical pause sums", treatment)
if (logYPlotter == "Yes") {
  p <- p + scale_y_continuous(trans=log10_trans(),
                breaks=trans_breaks('log10', function(x) 10^x),
                labels=trans_format('log10', math_format(.x)))
 p <- p + labs(y=bquote(log[10](.(yLabel))))
}
else {
   p <- p + ylab(yLabel)
}
if (logXPlotter == "Yes") {
  p <- p + scale_x_continuous(trans=log10_trans(),
                breaks=trans_breaks('log10', function(x) 10^x),
                labels=trans_format('log10', math_format(.x))) 
 p <- p + labs(x=bquote(log[10](.(xLabel))))
} else {
   p <- p + xlab(xLabel)
}
p <- p + theme_bw() + fontTheme
if (!is.na(plotTitle)) {
 p <- p + labs(title = plotTitle) +
         theme(plot.title = element_text(size = 28, face=1.5, hjust=0.5))
}  

p <- p + theme(
    text = element_text(size = 22),
     axis.text = element_text(size=23),
      axis.title.y = element_text(size=30),
      axis.title.x = element_text(size=28),
    #legend.title = element_blank()
  )
Width = 20
if (showLegend == "False") {
        p <- p + theme(legend.position = "None")
Width = 14
}
   ggsave(filename = sprintf('%s/%s_%s_fitted_vs_actual_uptorank_%d_%s_showLegend_%s_75percentile_%s_logX_%s_logY_%s_maxSaturation_%s_restOfFuncEval_%s_restOfData_Rplotter.pdf', outputDir, datasetName, typeofFit, indexThresh, showLegend, is75percentile, logXPlotter, logYPlotter, plotMaxSaturation, plotRestOfFunction, plotRestOfData), plot = p, width = Width, height = 8)
        }
    }
}


