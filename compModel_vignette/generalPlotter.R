library(ggplot2)
library(dplyr)
library(ggsci)
library(RColorBrewer)
library(scales)
# https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html
source("/Users/rudradeepmukherjee/Documents/UConn Health/comparment.model.project/datasets, genome/plotting_composites_lattice.R")
source("/Users/rudradeepmukherjee/Documents/UConn Health/comparment.model.project/datasets, genome/residualAndFitPlotter.R")

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

summaryFunction <- function(rateData) {
    rateData$P_base = rateData$kinit_base / (rateData$kpre + rateData$krel_base)
    rateData$P_expr = rateData$kinit_expr / (rateData$kpre + rateData$krel_expr)
    rateData$B_base = ( rateData$krel_base * rateData$P_base) / (rateData$kelong)
    rateData$B_expr = ( rateData$krel_expr * rateData$P_expr) / (rateData$kelong)
    rateData$kinit_FC = rateData$kinit_expr / rateData$kinit_base
    rateData$krel_FC = rateData$krel_expr / rateData$krel_base
  summary.rateData <- group_by(rateData, gene) %>%  summarise_at(c("kinit_base", "kinit_expr", "krel_base", "krel_expr", "kpre", "kelong", "kinit_FC", "krel_FC", "P_base", "P_expr", "B_base", "B_expr"), c("min", "max", "mean", "var"))
  summary.rateData$P_FC = summary.rateData$P_expr_mean / summary.rateData$P_base_mean
  summary.rateData$B_FC = summary.rateData$B_expr_mean / summary.rateData$B_base_mean

  return(summary.rateData)
}

plotBandP_foldchange <- function(summaryObj, filename, condnIdentifier = "treatment", WIDTH=7, filetype="pdf") {

         p<- ggplot(summaryObj[,c("gene", "P_FC", "B_FC")]) +
                 geom_violin(aes(x="Pause FC", y = P_FC), scale="width") +
                 geom_violin(aes(x="Body FC", y = B_FC), scale="width") +
                 geom_boxplot(aes(x="Pause FC", y = P_FC), width = 0.3) +
                 geom_boxplot(aes(x="Body FC", y = B_FC), width = 0.3) +
                 geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
                 scale_y_continuous(trans=log2_trans(), breaks=trans_breaks('log2', function(x) 2^x), labels=trans_format('log2', math_format(.x))) + 
                 ylab(expression(log[2]~"values")) + 
                 ggtitle(condnIdentifier) + 
                 theme_bw() + fontTheme +
                 theme(
              plot.title = element_text(size=23, face=1.2,hjust=0.5),
              #strip.text = element_text(size=24),
              axis.text.x = element_text(size=28, face=1.5), #angle=30, vjust=0.9, hjust=1),
              axis.title.x = element_blank(),
              #axis.title.x  = element_text(size=27, face=1.5, hjust = 0.5),
              axis.title.y = element_text(size=32),
              legend.title = element_blank()
              )
    ggsave(paste0(filename,'.',filetype), width=WIDTH, plot=p)
}

# each of control.df and treat.df are rate.csv files from Unimod run without scale factors and steric hindrance - https://github.com/CshlSiepelLab/UniMod
siepel_pauserelease_plotter <- function(control.df, treat.df, cond1="control", cond2="treat", kelong_vals_permin = c(2000, 1800,2700,3600), filetype="pdf",WIDTH=10) {
     control.df$type = cond1
     treat.df$type = cond2
     data.df = rbind(control.df, treat.df)
     plot.df = as.data.frame(c())
     for (kelong_val in kelong_vals_permin) {
       data.df$kelong = kelong_val
       data.df$krel_betaorg_permin = data.df$kelong *  data.df$beta_org
       data.df$krel_betaadp_permin = data.df$kelong *  data.df$beta_adp 
       plot.df = rbind(plot.df, data.df)
     }
     plot.df.org = plot.df[,c("gene_id", "krel_betaorg_permin", "kelong", "type")]
     plot.df.adp = plot.df[,c("gene_id", "krel_betaadp_permin", "kelong", "type")]
     plot.df.org$estimate = "original"
     plot.df.adp$estimate = "variable pause region"
     colnames(plot.df.org) = c("gene", "krel_permin", "kelong", "type", "estimate")
     colnames(plot.df.adp) = c("gene", "krel_permin", "kelong", "type", "estimate")
    plot.df.long = rbind(plot.df.org, plot.df.adp)
    plot.df.long$kelong = sprintf("%s",plot.df.long$kelong)
     p <- ggplot(plot.df.long, aes(x = kelong, y = krel_permin, color = type)) +
            geom_violin(position = "dodge", scale="width") +
            geom_boxplot(position=position_dodge(width = 0.9), width=0.3, alpha = 0.4) + 
            facet_wrap(~estimate)
    p <- p + scale_y_continuous(trans=log10_trans(),
          breaks=trans_breaks('log10', function(x) 10^x),
          labels=trans_format('log10', math_format(.x)))
    p <- p + ylab(expression(log[10] ~ "Pause escape"~beta~"x"~zeta~"(minute"^"-1"~")")) +
            xlab(expression(italic(k[elong])~"bp/min")) + 
            theme_bw() + fontTheme +
            theme(
            #  plot.title = element_text(size=23, face=1.2,hjust=0.5),
              strip.text = element_text(size=24),
              axis.text.x = element_text(size=28, face=1.5), #angle=30, vjust=0.9, hjust=1),
              #axis.title.x = element_blank(),
              axis.title.x  = element_text(size=27, face=1.5, hjust = 0.5),
              axis.title.y = element_text(size=32),
              legend.title = element_blank()
        )
    filename = sprintf("%s_%s_krel_permin_siepel", cond2, cond1)
    ggsave(paste0(filename,'.',filetype), width=WIDTH, plot=p)
    return(plot.df.long)
}

plotKinitValues_multiDatasets <- function(summaryObj, filename, NCOL=1, baselineName = "control", condnShortHand = "treatment", kinit_scale="min", WIDTH=10, HEIGHT=10, maxKpre = 10, filetype="pdf") {
        df.control = summaryObj[,c("datasetName", "gene","B_base_mean")]
        colnames(df.control) = c("datasetName", "gene", "B")
        df.control$type = baselineName
        df.treat = summaryObj[,c("datasetName", "gene","B_expr_mean")]
        colnames(df.treat) = c("datasetName", "gene", "B")
        df.treat$type = condnShortHand
        df1 = as.data.frame(c())
        for (kelong_val in c(30,45,60)) {
           df.control$kelong = kelong_val
           df.treat$kelong = kelong_val
           df1 = rbind(df1, df.control, df.treat)
        }
        df1$limitType = sprintf("%s bp/sec", df1$kelong)
        df1$kinit = df1$B * df1$kelong
        df.P.control = summaryObj[,c("datasetName","gene","P_base_mean")]
        df.P.treat = summaryObj[,c("datasetName","gene","P_expr_mean")]
        colnames(df.P.control) = c("datasetName","gene", "P")
        colnames(df.P.treat) = c("datasetName","gene","P")
        df.P.control$type = baselineName
        df.P.treat$type = condnShortHand
        df2 = rbind(df.P.control, df.P.treat)
        df2$kinit = maxKpre * df2$P 
        df2$limitType = "max(kpre) x P"
        df = rbind(df1[,c("datasetName", "gene", "kinit","type", "limitType")], df2[,c("datasetName", "gene", "kinit","type", "limitType")])
     df$limitType = factor(df$limitType, levels = c("30 bp/sec", "45 bp/sec", "60 bp/sec", "max(kpre) x P"))
        if (kinit_scale == "min") {
          df$kinit = df$kinit * 60
        }
        p <- ggplot(df, aes(x = limitType, y = kinit, color = type)) +
        geom_violin(position = "dodge", scale="width") +
        geom_boxplot(position=position_dodge(width = 0.9), width=0.3, alpha = 0.4) +
        scale_y_continuous(trans=log10_trans(),
          breaks=trans_breaks('log10', function(x) 10^x),
          labels=trans_format('log10', math_format(.x))) + 
                facet_wrap(~datasetName, ncol=NCOL)  
        if (kinit_scale == "min") {
            p <- p + ylab(expression(log[10] ~ italic(k[init])~"values (minute"^"-1"~")"))
        }
        if (kinit_scale == "sec") {
            p <- p + ylab(expression(log[10] ~ italic(k[init])~"values (second"^"-1"~")"))
        }

           p <- p + scale_x_discrete(labels = c("30",expression(atop("45","(bp/sec)")),"60", expression(atop("max("~italic(k[pre])~")","x"~italic("P"[normalized]))))) + 
        xlab(expression(atop(italic(k[elong])~"x"~italic("B"[normalized]), "(low"~italic(k[pre])~")" ))) +
        theme_bw() + fontTheme +
         #ggtitle(plotTitle) + 
         theme(
              #plot.title = element_text(size=23, face=1.2,hjust=0.5),
              legend.position="top",
               strip.text= element_text(size=26),
               axis.text.x = element_text(size=28, face=1.5), #angle=30, vjust=0.9, hjust=1),
              #axis.title.x = element_blank(),
              axis.title.x  = element_text(size=27, face=1.5, hjust = 0.2, vjust = 1),
              axis.title.y = element_text(size=32),
              legend.title = element_blank(),
              legend.text = element_text(size=24)
        )
        #filename = sprintf("%s_kinitValuesLimits", condnShortHand)
        if (filetype == "png") {
           ggsave(paste0(filename,'.',filetype), plot=p, width= WIDTH, height=HEIGHT,dpi=300)
        } else {
                ggsave(paste0(filename,'.',filetype), width= WIDTH, height=HEIGHT, plot=p)
        }
}

plotKinitValues_multiDatasets_twoPanels <- function(summaryObj, filename, NCOL=2, baselineName = "control", condnShortHand = "treatment", kinit_scale="min", WIDTH=10, HEIGHT=10, maxKpre = 1, filetype="pdf") {
        df.control = summaryObj[,c("datasetName", "gene","B_base_mean")]
        colnames(df.control) = c("datasetName", "gene", "B")
        df.control$type = baselineName
        df.treat = summaryObj[,c("datasetName", "gene","B_expr_mean")]
        colnames(df.treat) = c("datasetName", "gene", "B")
        df.treat$type = condnShortHand
        df1 = as.data.frame(c())
        #for (kelong_val in c(30,45,60)) {
        for (kelong_val in c(2000 / 60)) {
           df.control$kelong = kelong_val
           df.treat$kelong = kelong_val
           df1 = rbind(df1, df.control, df.treat)
        }
        #df1$limitType = sprintf("%s bp/sec", df1$kelong)
        df1$limitType = "2000 bp/min"
        df1$kinit = df1$B * df1$kelong
        df.P.control = summaryObj[,c("datasetName","gene","P_base_mean")]
        df.P.treat = summaryObj[,c("datasetName","gene","P_expr_mean")]
        colnames(df.P.control) = c("datasetName","gene", "P")
        colnames(df.P.treat) = c("datasetName","gene","P")
        df.P.control$type = baselineName
        df.P.treat$type = condnShortHand
        df2 = rbind(df.P.control, df.P.treat)
        df2$kinit = maxKpre * df2$P 
        df2$limitType = "large[kpre] x P"
        df = rbind(df1[,c("datasetName", "gene", "kinit","type", "limitType")], df2[,c("datasetName", "gene", "kinit","type", "limitType")])
     #df$limitType = factor(df$limitType, levels = c("30 bp/sec", "45 bp/sec", "60 bp/sec", "large[kpre] x P"))
     df$limitType = factor(df$limitType, levels = c("2000 bp/min", "large[kpre] x P"))

        if (kinit_scale == "min") {
          df$kinit = df$kinit * 60
        }
        p <- ggplot(df, aes(x = limitType, y = kinit, color = type)) +
        geom_violin(position = "dodge", scale="width") +
        geom_boxplot(position=position_dodge(width = 0.9), width=0.3, alpha = 0.4) +
        scale_y_continuous(trans=log10_trans(),
          breaks=trans_breaks('log10', function(x) 10^x),
          labels=trans_format('log10', math_format(.x))) + 
                facet_wrap(~datasetName, ncol=NCOL)  
        if (kinit_scale == "min") {
            p <- p + ylab(expression(log[10] ~ italic(K[init])~"(RNAP x minute"^"-1"~")"))
        }
        if (kinit_scale == "sec") {
            p <- p + ylab(expression(log[10] ~ italic(K[init])~"(RNAP x second"^"-1"~")"))
        }

      #  p <- p + scale_x_discrete(labels = c("30",expression(atop("45","(bp/sec)")),"60", 
      #                                          expression(atop("large("~italic(k[pre])~")","x"~italic("P"[normalized]))))) + 
      #  xlab(expression(atop(italic(k[elong])~"x"~italic("B"[normalized]), "(low"~italic(k[pre])~")" ))) +
        p <- p + scale_x_discrete(labels = c(expression(atop(italic(k[elong])~"x"~italic("B"[normalized]),
                                                             atop("(low"~italic(k[pre])~")", italic(k[elong])~"= 2000 bp/min"))), 
                                                expression(atop("large("~italic(k[pre])~")","x"~italic("P"[normalized]))))) + 
                  theme_bw() + fontTheme +
         #ggtitle(plotTitle) + 
         theme(
              #plot.title = element_text(size=23, face=1.2,hjust=0.5),
              legend.position="top",
               axis.text.x = element_text(size=28, face=1.5), #angle=30, vjust=0.9, hjust=1),
              #axis.title.x = element_blank(),
              #axis.title.x  = element_text(size=27, face=1.5, hjust = 0.2, vjust = 1),
              strip.background = element_rect(fill = "#E6E6FA", color = "black", linewidth = 1.6),
    strip.text = element_text(color = "#6A5ACD", face = "bold", size = 21),
               axis.title.x  = element_blank(),
              axis.title.y = element_text(size=32),
              legend.title = element_blank(),
              legend.text = element_text(size=24)
        )
        #filename = sprintf("%s_kinitValuesLimits", condnShortHand)
        if (filetype == "png") {
           ggsave(paste0(filename,'.',filetype), plot=p, width= WIDTH, height=HEIGHT,dpi=300)
        } else {
                ggsave(paste0(filename,'.',filetype), width= WIDTH, height=HEIGHT, plot=p)
        }
       return(df)
}

plotKinitvsKpre <- function(summaryObj, condnShortHand, kelong_minute = 2000, baselineName = "control", kinit_scale="min", WIDTH=7, filetype = "pdf") {
   df.control = summaryObj[,c("gene", "P_base", "B_base", "kpre")]
   df.treat  = summaryObj[,c("gene", "P_expr", "B_expr", "kpre")]
   colnames(df.control) = c("gene", "P", "B", "kpre")
   colnames(df.treat) = c("gene", "P", "B", "kpre")
   df.control$cond = baselineName
   df.treat$cond = condnShortHand
   df = rbind(df.control, df.treat)
   df$kelong = kelong_minute / 60
   df$kinit = df$kpre * df$P + df$B * df$kelong
   if (kinit_scale == "min") {
      df$kinit = df$kinit * 60
   }
   df$kpre_str = ""
   df[df$kpre == 0.001,]$kpre_str = "0.001"
   df[df$kpre == 0.01,]$kpre_str = "0.01"
   df[df$kpre == 0.1,]$kpre_str = "0.1"
   df[df$kpre == 1,]$kpre_str = "1"
   print(head(df))
   df$cond = factor(df$cond, levels = c(baselineName, condnShortHand))
   #df$kpre = factor(df$kpre, levels = unique(df$kpre))
   p <- ggplot(df, aes(x=kpre_str, y=kinit, color=cond)) +
           geom_violin(  position = "dodge", scale="width") + 
           geom_boxplot( position=position_dodge(width = 0.9), width=0.3, alpha = 0.4) +
                  #facet_wrap(~cond, ncol=1) +
        #   scale_x_continuous(trans=log10_trans(),
        #  breaks=trans_breaks('log10', function(x) 10^x),
        #  labels=trans_format('log10', math_format(.x))) +
           scale_y_continuous(trans=log10_trans(),
          breaks=trans_breaks('log10', function(x) 10^x),
          labels=trans_format('log10', math_format(.x))) 
           if (kinit_scale == "min") {
            p <- p + ylab(expression(log[10] ~ italic(k[init])~"values (minute"^"-1"~")"))
        } else {
            p <- p + ylab(expression(log[10] ~ italic(k[init])~"values ("~kinit_scale^"-1"~")"))
          }
           #p <- p + xlab(expression(log[10] ~ italic(k[pre]) ~"(second"^"-1"~")")) + 
           p <- p + xlab(expression(italic(k[pre]) ~"(second"^"-1"~")")) + 
                   theme_bw() + fontTheme + 
                   theme(legend.title = element_blank(),
                         strip.text = element_text(size=24))

     filename = sprintf("%s_kinitvskpre", condnShortHand)
        if (filetype == "png") {
           ggsave(paste0(filename,'.',filetype), plot=p, width=WIDTH, dpi=300)
        } else {
                ggsave(paste0(filename,'.',filetype), width=WIDTH, plot=p)
        }
}

plotKinitValuesForAllGenes <- function(summaryObj, condnShortHand, plotTitle, baselineName = "control", kinit_scale="min", maxKpre = 10, filetype="pdf") {
        df.control = summaryObj[,c("gene","B_base_mean")]
        colnames(df.control) = c("gene", "B")
        df.control$type = baselineName
        df.treat = summaryObj[,c("gene","B_expr_mean")]
        colnames(df.treat) = c("gene", "B")
        df.treat$type = condnShortHand
        df1 = as.data.frame(c())
        for (kelong_val in c(30,45,60)) {
           df.control$kelong = kelong_val
           df.treat$kelong = kelong_val
           df1 = rbind(df1, df.control, df.treat)
        }
        df1$limitType = sprintf("%s bp/sec", df1$kelong)
        df1$kinit = df1$B * df1$kelong
        df.P.control = summaryObj[,c("gene","P_base_mean")]
        df.P.treat = summaryObj[,c("gene","P_expr_mean")]
        colnames(df.P.control) = c("gene", "P")
        colnames(df.P.treat) = c("gene","P")
        df.P.control$type = baselineName
        df.P.treat$type = condnShortHand
        df2 = rbind(df.P.control, df.P.treat)
        df2$kinit = maxKpre * df2$P 
        df2$limitType = "max(kpre) x P"
        df = rbind(df1[,c("gene", "kinit","type", "limitType")], df2[,c("gene", "kinit","type", "limitType")])
     df$limitType = factor(df$limitType, levels = c("30 bp/sec", "45 bp/sec", "60 bp/sec", "max(kpre) x P"))
        if (kinit_scale == "min") {
          df$kinit = df$kinit * 60
        }
        p <- ggplot(df) +
        geom_violin(aes(x = limitType, y = kinit, color = type), position = "dodge", scale="width") +
        geom_boxplot(aes(x = limitType, y = kinit, color = type), position=position_dodge(width = 0.9), width=0.3, alpha = 0.4) +
        scale_y_continuous(trans=log10_trans(),
          breaks=trans_breaks('log10', function(x) 10^x),
          labels=trans_format('log10', math_format(.x)))  
          if (kinit_scale == "min") {
            p <- p + ylab(expression(log[10] ~ italic(k[init])~"values (minute"^"-1"~")"))
        }
        if (kinit_scale == "sec") {
            p <- p + ylab(expression(log[10] ~ italic(k[init])~"values (second"^"-1"~")"))
        }

           p <- p + scale_x_discrete(labels = c("30",expression(atop("45","(bp/sec)")),"60", expression(atop("max("~italic(k[pre])~")","x"~italic("P"[normalized]))))) + 
        xlab(expression(atop(italic(k[elong])~"x"~italic("B"[normalized]), "(low"~italic(k[pre])~")" ))) +
        theme_bw() + fontTheme +
         ggtitle(plotTitle) + 
         theme(
              plot.title = element_text(size=23, face=1.2,hjust=0.5),
              axis.text.x = element_text(size=28, face=1.5), #angle=30, vjust=0.9, hjust=1),
              #axis.title.x = element_blank(),
              axis.title.x  = element_text(size=27, face=1.5, hjust = 0.2, vjust = 1),
              axis.title.y = element_text(size=32),
              legend.title = element_blank()
        )
        filename = sprintf("%s_kinitValuesLimits", condnShortHand)
        if (filetype == "png") {
           ggsave(paste0(filename,'.',filetype), plot=p, width=10, dpi=300)
        } else {
                ggsave(paste0(filename,'.',filetype), width=10, plot=p)
        }
}

plotKinitLimits <- function(summaryObj, condnShortHand, gene.name) {
  summaryObj.gene = subset(summaryObj, gene == gene.name)
  df = c()
  df.kinit = summaryObj.gene[,c("B_base_mean", "P_base_mean")]
  df.kinit$kelong = 45
  for (kpre_val in c(1e-6, 1e-5,1e-4,1e-3,1e-2,1e-1,1,10)) {
    df.kinit$kpre = kpre_val
    df = rbind(df, df.kinit)
  }
  df$kinit = df$kpre * df$P_base_mean  + df$kelong * df$B_base_mean
 btimeselongintercept =  45 * summaryObj.gene$B_base_mean
 kpreTimesP = 10 * summaryObj.gene$P_base_mean
  p <- ggplot(df, aes(x = kpre, y = kinit)) +
         geom_point(color="blue", size=2) + 
          geom_line(color="orange", linewidth=1.2) +
        geom_hline(yintercept = btimeselongintercept, linetype = "dashed", color = "red") +  
        annotate("text", x = 1e-4, y = btimeselongintercept, label = 'B x Kelong',
            hjust = 1.1, vjust = -0.5, color = "red", size = 5) +
        geom_hline(yintercept = kpreTimesP, linetype = "dashed", color = "red") +  
        annotate("text", x = 10, y = kpreTimesP, label = bquote(italic("10"~(k[pre])~"x P") == .(kpreTimesP)),
            hjust = 1.1, vjust = 1, color = "red", size = 5) + 
       scale_y_continuous(trans=log10_trans(),
          breaks=trans_breaks('log10', function(x) 10^x),
          labels=trans_format('log10', math_format(.x))) +
        scale_x_continuous(trans=log10_trans(),
          breaks=trans_breaks('log10', function(x) 10^x),
          labels=trans_format('log10', math_format(.x))) +
        ylab(expression(log[10] ~ italic(k[init]))) +
        xlab(expression(italic(k[elong]))) +
                  theme_bw() + fontTheme
  filename = sprintf("%s_kinitvskpre_gene_%s", condnShortHand, gene.name)
  if (filetype == "png") {
        ggsave(paste0(filename,'.',filetype), plot=p, width=8, dpi=300)
  } else {
        ggsave(paste0(filename,'.',filetype), width=8, plot=p)
    }

}

printMedianAndQuartileRange <- function(dataVals, dataIdentiifer, quartileProbs = c(0.25,0.50,0.75)) {
  percentiles <- quantile(dataVals, probs = quartileProbs)
return(percentiles)
}
getGeometricMean <- function(dataVals) {
    # Calculate geometric mean
    geometric_mean <- exp(mean(log(dataVals)))
    return(geometric_mean)
}
# Calculate geometric standard deviation
getGeometricSD <- function(dataVals) {
        geometric_sd <- exp(sd(log(dataVals)))
        return(geometric_sd)
}
## 95 percent confidence interval
getGeometricConfidenceIntvl <- function(dataVals) {
        geometric_mean <- exp(mean(log(dataVals)))

        # Standard error of the log-transformed dataVals
        stderr_log <- sd(log(dataVals)) / sqrt(length(dataVals))

        # Confidence interval on the log scale
        conf_interval_log <- mean(log(dataVals)) + c(-1, 1) * qnorm(0.975) * stderr_log

        # Convert back to the original scale
        conf_interval <- exp(conf_interval_log)

        return(conf_interval)
}
getArithConfInterval <- function(dataVals) {
   stderr_arith = sd(dataVals) / sqrt(length(dataVals))
   conf_interval <- mean(dataVals) + c(-1, 1) * qnorm(0.975) * stderr_arith
   return(conf_interval)
}
printGeometricSummaryOfValues <- function(dataVals, identifier, scale="minute") {
   if (!is.na(scale) && scale == "minute") {
      dataVals = dataVals * 60
   print(sprintf("geometric summary of %s, timescale in minute", identifier))
   } else {
           if (is.na(scale)) {
             print(sprintf("arithmetic summary of %s, timescale N/A", identifier))
         } else {
             print(sprintf("arithmetic summary of %s, timescale %s", identifier, scale))

           }
   }
   iqrvals = printMedianAndQuartileRange(dataVals)
   print(sprintf("median: %0.1f, IQR = %0.1f - %0.1f = %0.1f", iqrvals[2][[1]], iqrvals[3][[1]], iqrvals[1][[1]], iqrvals[3][[1]] - iqrvals[1][[1]]))
   geom_mean = getGeometricMean(dataVals)
   geom_sd = getGeometricSD(dataVals)
   geom_conf_interv = getGeometricConfidenceIntvl(dataVals)
   print(sprintf("geometric mean: %0.1f, sd: %0.1f, 95%% conf: %0.1f - %0.1f", geom_mean, geom_sd, geom_conf_interv[1], geom_conf_interv[2])) 
}

printArithmeticSummaryOfValues <- function(dataVals, identifier, scale="minute") {
   if (!is.na(scale) && scale == "minute") {
      dataVals = dataVals * 60
      print(sprintf("arithmetic summary of %s, timescale in minute", identifier))
   } else {
           if (is.na(scale)) {
             print(sprintf("arithmetic summary of %s, timescale N/A", identifier))
         } else {
             print(sprintf("arithmetic summary of %s, timescale %s", identifier, scale))
           
           }
   }
   iqrvals = printMedianAndQuartileRange(dataVals)
   print(sprintf("median: %0.1f, IQR = %0.1f - %0.1f = %0.1f", iqrvals[2][[1]], iqrvals[3][[1]], iqrvals[1][[1]], iqrvals[3][[1]] - iqrvals[1][[1]]))
   geom_mean = mean(dataVals)
   geom_sd = sd(dataVals)
   geom_conf_interv = getArithConfInterval(dataVals)
   print(sprintf("arithmetic mean: %0.1f, sd: %0.1f, 95%% conf: %0.1f - %0.1f", geom_mean, geom_sd, geom_conf_interv[1], geom_conf_interv[2])) 
}

plotBtimesKelong <- function(summaryRates, condnShortHand, quartileProbs = c(0, 0.25,0.5,0.75,1),baseline="control",filetype= "pdf") {
        df.elong.30 = summaryRates[,c("B_base_mean", "P_base_mean")]
        df.elong.30$elongation = 30
        df.elong.45 = summaryRates[,c("B_base_mean", "P_base_mean")]
        df.elong.45$elongation = 45
        df.elong.60 = summaryRates[,c("B_base_mean", "P_base_mean")]
        df.elong.60$elongation = 60
        df.elongtimesB_base = rbind(df.elong.30, df.elong.45, df.elong.60)
        df.elongtimesB_base$type = baseline
        df.elong.30 = summaryRates[,c("B_expr_mean", "P_expr_mean")]
        df.elong.30$elongation = 30
        df.elong.45 = summaryRates[,c("B_expr_mean", "P_expr_mean")]
        df.elong.45$elongation = 45
        df.elong.60 = summaryRates[,c("B_expr_mean", "P_expr_mean")]
        df.elong.60$elongation = 60
        df.treat = rbind(df.elong.30, df.elong.45, df.elong.60)
        df.treat$type = condnShortHand
        colnames(df.elongtimesB_base) = c("B", "P","elongation", "type")
        colnames(df.treat) = c("B", "P", "elongation", "type")
        df.elongtimesB = rbind(df.elongtimesB_base, df.treat)
        df.elongtimesB$product = df.elongtimesB$B * df.elongtimesB$elongation
        df.elongtimesB$elongationLabel = sprintf("%s bp/sec", df.elongtimesB$elongation)
       for (cond in c(baseline, condnShortHand)) {
         for (elong in c(30,45,60) ) {
              subsetData = subset(df.elongtimesB, elongation == elong & type == cond )$product
              percentiles <- quantile(subsetData, probs = quartileProbs)
              print(sprintf("quartiles, %s elong: %d", cond, elong))
              print(percentiles)
              print(sprintf("geometric mean, %f", exp(mean(log(subsetData)))))
           } 
         }    
        p <- ggplot(df.elongtimesB) +
        geom_violin(aes(x = elongationLabel, y = product, color = type), position = "dodge", scale="width") +
        geom_boxplot(aes(x = elongationLabel, y = product, color = type), position=position_dodge(width = 0.9), width=0.3, alpha = 0.4) +
        scale_y_continuous(trans=log10_trans(),
          breaks=trans_breaks('log10', function(x) 10^x),
          labels=trans_format('log10', math_format(.x)), limits=c(NA,1)) +
                  ylab(expression(log[10] ~ italic(k[elong])~"x Body density")) +
        xlab(expression(italic(k[elong]))) +
                  theme_bw() + fontTheme +
        theme(
              #legend.position="none",
              axis.text.x = element_text(size=28, face=1.5), #angle=30, vjust=0.9, hjust=1),
              #axis.text.x = element_blank(),
              #axis.title.x = element_blank(),
              axis.title.x = element_text(size=30),
              axis.title.y = element_text(size=32),
              legend.title = element_blank()
        )
        filename = sprintf("%s_Btimeskelong", condnShortHand)
        if (filetype == "png") {
           ggsave(paste0(filename,'.',filetype), plot=p, width=10, dpi=300)
        } else {
                ggsave(paste0(filename,'.',filetype), width=10, plot=p)
        }
        ## plot [ B x Kelong + kpre * P ]
        df.kpre  = as.data.frame(c())
                df.elongtimesB$elongationLabel = sprintf("%s", df.elongtimesB$elongation)
        for (kpre_val in c(1e-6,1e-5,1e-4, 1e-3,1e-2,0.1,1,10)) {
          df.elongtimesB$kpre = kpre_val
          df.elongtimesB$kinit = df.elongtimesB$product + df.elongtimesB$P * df.elongtimesB$kpre
          df.kpre = rbind(df.kpre, df.elongtimesB)
        }
        p <- ggplot(df.kpre, aes(x = elongationLabel,  y = kinit, color = type)) +
                 geom_violin( position = "dodge", scale="width") + 
                 geom_boxplot(position=position_dodge(width = 0.9), width=0.3, alpha = 0.4) +
                scale_y_continuous(trans=log10_trans(),
          breaks=trans_breaks('log10', function(x) 10^x),
          labels=trans_format('log10', math_format(.x))) +
                  ylab(expression(log[10] ~ italic(k[init]))) +
                 facet_wrap(~ kpre, labeller = label_bquote(italic(k[pre]) == .(kpre)), nrow=2) +
                  xlab(expression(italic(k[elong])~"(bp/sec)")) +
                  theme_bw() + fontTheme +
        theme(
              #legend.position="none",
              axis.text.x = element_text(size=28, face=1.5), #angle=30, vjust=0.9, hjust=1),
              #axis.text.x = element_blank(),
              #axis.title.x = element_blank(),
              axis.title.x = element_text(size=30),
              axis.title.y = element_text(size=32),
              strip.text = element_text(size=24),
              legend.title = element_blank()
        )
     filename = sprintf("%s_kinitvskpre_kelong", condnShortHand)
        if (filetype == "png") {
           ggsave(paste0(filename,'.',filetype), plot=p, width=12, dpi=300)
        } else {
                ggsave(paste0(filename,'.',filetype), width=12, plot=p)
        }
       p <- ggplot(df.kpre, aes(x = elongationLabel,  y = kinit / kpre, color = type)) +
                 geom_violin( position = "dodge", scale="width") + 
                 geom_boxplot(position=position_dodge(width = 0.9), width=0.3, alpha = 0.4) +
                scale_y_continuous(trans=log10_trans(),
          breaks=trans_breaks('log10', function(x) 10^x),
          labels=trans_format('log10', math_format(.x))) +
                  ylab(expression(log[10] ~ frac(italic(k[init]), italic(k[pre])))) +
                 facet_wrap(~ kpre, labeller = label_bquote(italic(k[pre]) == .(kpre)), nrow=2) +
                  xlab(expression(italic(k[elong])~"(bp/sec)")) +
                  theme_bw() + fontTheme +
        theme(
              #legend.position="none",
              axis.text.x = element_text(size=28, face=1.5), #angle=30, vjust=0.9, hjust=1),
              #axis.text.x = element_blank(),
              #axis.title.x = element_blank(),
              axis.title.x = element_text(size=30),
              axis.title.y = element_text(size=32),
              strip.text = element_text(size=24),
              legend.title = element_blank()
        )
     filename = sprintf("%s_kinitBYkpre_kelong", condnShortHand)
        if (filetype == "png") {
           ggsave(paste0(filename,'.',filetype), plot=p, width=12, dpi=300)
        } else {
                ggsave(paste0(filename,'.',filetype), width=12, plot=p)
        }

}

makeQQplots <- function(dataVals, filename, title = "QQplot", lineCol = "blue") {
pdf(paste0(filename, ".pdf")) 
par(cex.axis = 1.5, cex.lab = 1.8, cex.main = 2.2)
qqnorm(dataVals, 
       main = title)
qqline(dataVals, col = "red")
dev.off()
}

absoluteRatePlotter <- function(summaryObj, filename, condnIdentifier, condnShortHand, constraintType="SteurerVock", filetype="png", pointALPHA=0.09) {

        p <- ggplot(summaryObj) +
#ggplot(newSummary.TRP.Kpre.constrained) +
    #geom_violin(aes(x="release_base_mean",y= Release_base_mean), scale="width") +
     #   geom_violin(aes(x="release_expr_mean",y= Release_expr_mean), scale="width") +
        geom_violin(aes(x="kinit_base_mean",y= kinit_base_mean), scale="width") +
        geom_violin(aes(x="kinit_expr_mean",y= kinit_expr_mean), scale="width") +
        geom_violin(aes(x="kpre_mean",y= kpre_mean), scale="width") +
        geom_violin(aes(x="krel_base_mean",y= krel_base_mean), scale="width") +
        geom_violin(aes(x="krel_expr_mean",y= krel_expr_mean), scale="width") +
      #  geom_violin(aes(x="terminal_base_mean",y= Terminal_base_mean), scale="width") +
      #  geom_violin(aes(x="terminal_expr_mean",y= Terminal_expr_mean), scale="width") +

      #  geom_point(aes(x="release_base_mean",y= Release_base_mean), position="jitter", alpha=pointALPHA) +
       # geom_point(aes(x="release_expr_mean",y= Release_expr_mean), position="jitter", alpha=pointALPHA) +
        geom_point(aes(x="kinit_base_mean",y= kinit_base_mean), position="jitter", alpha=pointALPHA) +
        geom_point(aes(x="kinit_expr_mean",y= kinit_expr_mean), position="jitter", alpha=pointALPHA) +
        geom_point(aes(x="kpre_mean",y= kpre_mean), position="jitter", alpha=pointALPHA) +
        geom_point(aes(x="krel_base_mean",y= krel_base_mean), position="jitter", alpha=pointALPHA) +
        geom_point(aes(x="krel_expr_mean",y= krel_expr_mean), position="jitter", alpha=pointALPHA) +
        #geom_point(aes(x="terminal_base_mean",y= Terminal_base_mean), position="jitter", alpha=pointALPHA) +
        #geom_point(aes(x="terminal_expr_mean",y= Terminal_expr_mean), position="jitter", alpha=pointALPHA) +

        #geom_boxplot(aes(x="release_base_mean",y= Release_base_mean), alpha=0.8, width=0.36, linewidth=0.3) +
        #geom_boxplot(aes(x="release_expr_mean",y= Release_expr_mean), alpha=0.8, width=0.36, linewidth=0.3)+
        geom_boxplot(aes(x="kinit_base_mean",y= kinit_base_mean), alpha=0.8, width=0.36, linewidth=0.3) +
        geom_boxplot(aes(x="kinit_expr_mean",y= kinit_expr_mean), alpha=0.8, width=0.36, linewidth=0.3) +
        geom_boxplot(aes(x="kpre_mean",y= kpre_mean), alpha=0.8, width=0.36, linewidth=0.3) +
        geom_boxplot(aes(x="krel_base_mean",y= krel_base_mean), alpha=0.8, width=0.36, linewidth=0.3) +
        geom_boxplot(aes(x="krel_expr_mean",y= krel_expr_mean), alpha=0.8, width=0.36, linewidth=0.3) +
        #geom_boxplot(aes(x="terminal_base_mean",y= Terminal_base_mean), alpha=0.8, width=0.36, linewidth=0.3) +
        #geom_boxplot(aes(x="terminal_expr_mean",y= Terminal_expr_mean), alpha=0.8, width=0.36, linewidth=0.3)+
                  scale_y_continuous(trans=log10_trans(),
                breaks=trans_breaks('log10', function(x) 10^x),
                labels=trans_format('log10', math_format(.x))) 
        p <- p +
                  scale_x_discrete(labels = c(expression(atop(italic(K[init]), "(control)")), bquote(atop(italic(K[init]), .(condnShortHand))),
                                     expression(italic(K[pre])),expression(atop(italic(K[rel]), "(control)")),bquote(atop(italic(K[rel]), .(condnShortHand))))) + 
                  ylab(expression(log[10] ~ "Mean Values (1 / second)")) 
         if (constraintType == "SteurerVock") {
           p <- p + ggtitle(bquote(.(condnIdentifier)~","~italic(k[pre])~"= 8.67 x"~10^{"-6"}~"to"~"0.023 / second"))
                  # ggtitle(expression({{condnIdentifier}}~italic(k[pre])~"= 8.67 x"~10^{"-6"}~"to"~"0.023 / second")) +
         ## condnIdentifier
         } else {
            p <- p + ggtitle(bquote(.(condnIdentifier)~","~italic(k[pre])~"="~10^{-8}~"to"~"10 / second")) 
         }
         #ylab(expression(log[10] ~ "Mean Raw Values - All Genes")) +
         p <- p + theme_bw() + fontTheme +
         theme(plot.title = element_text(size=23, face=1.2,hjust=0.5), axis.title.x = element_blank(),
         axis.text.x = element_text(size=32, face=1, vjust=0.9))
        if (filetype == "png") {
           ggsave(paste0(filename,'.',filetype), plot=p, height=8, width=11.7, dpi=300)
        } else {
                ggsave(paste0(filename,'.',filetype), plot=p, height=8, width=11.7)
        }

}

krelFoldchangePlotter_twoSeparateViolins <- function(summaryObj, filename, condnIdentifier, condnShortHand, constraintType="SteurerVock", filetype="png", pointALPHA=0.2) {
 summaryObj$krelFCcategory = "Undefined"
summaryObj[summaryObj$krel_FC_max < 1, ]$krelFCcategory = "krelFC_less_1"
summaryObj[summaryObj$krel_FC_min >= 1, ]$krelFCcategory = "krelFC_more_1"

        df.krel.kinit = as.data.frame(cbind(c("krelFC_less_1", "krelFC_more_1"),
                             c(nrow(subset(summaryObj, krelFCcategory == "krelFC_less_1")) ,
                   nrow(subset(summaryObj, krelFCcategory == "krelFC_more_1")))))
  colnames(df.krel.kinit) = c("typeOfData", "countOfData")
  df.krel.kinit$countOfData = as.numeric(df.krel.kinit$countOfData)
  sum_all = sum(df.krel.kinit$countOfData)

  df.krel.kinit$percentLabel = sprintf("%0.2f%%",100.0*(df.krel.kinit$countOfData / sum_all))
YCOORD = max(summaryObj$krel_FC_min) + 2
        
        p <- ggplot(summaryObj) +
        geom_violin(aes(x= krelFCcategory, y=krel_FC_min), width=0.5) +
        geom_boxplot(aes(x= krelFCcategory, y=krel_FC_min), width=0.3) +
        geom_point(aes(x= krelFCcategory, y=krel_FC_min), position="jitter", alpha=pointALPHA) +
        geom_label(data=df.krel.kinit, aes(label=percentLabel,x = typeOfData, y = YCOORD),size=6,size.unit="mm") +
       scale_y_continuous(trans=log2_trans(), breaks=trans_breaks('log2', function(x) 2^x), labels=trans_format('log2', math_format(.x))) +
       scale_x_discrete(labels = c(expression(atop("Decrease in",italic(k[rel]))), expression(atop("Increase in",italic(k[rel]))) )) +
              ggtitle(condnIdentifier) +
       ylab(expression(log[2] ~ "Foldchange (FC) in" ~ italic(K[rel]))) +
        theme_bw() + fontTheme +
        theme(legend.position="none",
              plot.title = element_text(size=28, face=1.5,  hjust = 0.5),
              axis.text.x = element_text(size=26, face=1.5), #vjust=0.9, hjust=1),
              #axis.text.x = element_blank(),
              axis.title.x = element_blank(), axis.title.y = element_text(size=32), legend.title = element_blank())

         if (filetype == "png") {
           ggsave(paste0(filename,'.',filetype), plot=p, dpi=300, width=7)
        } else {
                ggsave(paste0(filename,'.',filetype), plot=p, width=7)
        }

}
krelFC_barplots_multiDatasets <- function(summaryObjMultiData, filename, labelDatasets,filetype="pdf", WIDTH=8, HEIGHT=9,bar.colors = c("#2ca02c", "#d62728"), REVERSE=TRUE, addTotalLabel = TRUE, VJUST = 0.5, HJUST=0.8,VJUST_X = 0.9, textsize=10) {
if (any(is.na(unique(summaryObjMultiData$datasetName)))) {
              print(sprintf("atlease one `datasetName` in the dataframe submitted is NA. Can't proceed."))
              return
     }
 summaryObjMultiData$krelFCcategory = "Undefined"
 summaryObjMultiData[summaryObjMultiData$krel_FC_max < 1, ]$krelFCcategory = "krelFC_less_1"
 summaryObjMultiData[summaryObjMultiData$krel_FC_min >= 1, ]$krelFCcategory = "krelFC_more_1"
 summary_data <- summaryObjMultiData %>%
        group_by(datasetName, krelFCcategory) %>%
        summarize(numGenes = n(),
                  ) %>%
        ungroup() %>%
        group_by(datasetName) %>%
  mutate(totalGenes = sum(numGenes)) %>%
  mutate(proportion = numGenes / totalGenes) %>%
  #select(-totalGenes) %>%
   ungroup()
   summary_data$krelFCcategory = factor(summary_data$krelFCcategory, levels = c("krelFC_more_1", "krelFC_less_1"))
   p <- ggplot(summary_data, aes(x = datasetName, y = proportion, fill = krelFCcategory)) +
       geom_bar(stat = "identity", position = position_fill(reverse = REVERSE)) +
         scale_fill_manual(labels = c(expression("Increase" ~ italic(k[rel])), expression("Decrease" ~ italic(k[rel]))), values = bar.colors) +
  labs(y = "Proportion of genes")
  if (addTotalLabel) {
          datasets = levels(summaryObjMultiData$datasetName)
          for (dataset in datasets) {
                  totalgenes = subset(summary_data, datasetName == dataset)$totalGenes[1]
                  p <- p +  geom_text(x=dataset, y = 1, vjust = VJUST, position = "identity", family = "sans",label = totalgenes, size=textsize)
          }
          }   
p <- p + scale_x_discrete(labels = labelDatasets) + 
        theme_bw() + fontTheme +
        theme(
              #legend.position="none",
              #plot.title = element_text(size=28, face=1.5,  hjust = 0.5),
              axis.text.x = element_text(size=26, face=1.2, angle = 30, hjust=HJUST, vjust=VJUST_X), #vjust=0.9, hjust=1),
              #axis.text.x = element_blank(),
              axis.title.x = element_blank(), 
              axis.title.y = element_text(size=32), 
              legend.title = element_blank())

     if (filetype == "png") {
           ggsave(paste0(filename,'.',filetype), plot=p, dpi=300, width=WIDTH, height= HEIGHT)
        } else {
                ggsave(paste0(filename,'.',filetype), plot=p, width=WIDTH, height= HEIGHT)
        }
 
}
# please factor order summaryObjMultiData$datasetName to the order you desire from left to right
# labelDatasets will contain the same order of names
krelFCplotter_multidatasets <- function(summaryObjMultiData, filename, labelDatasets, Y_DOWN_arr = c(1/4.5,1/4.5), Y_UP_arr = c(4, 4, 4), filetype="png", H_JUST = 0.5, showPercentLabels = TRUE, showPoints = TRUE, pointALPHA=0.2, WIDTH=10, VIOLIN_WIDTH = 0.5) {
  if (any(is.na(unique(summaryObjMultiData$datasetName)))) {
              print(sprintf("atlease one `datasetName` in the dataframe submitted is NA. Can't proceed."))
              return 
     }
 summaryObjMultiData$krelFCcategory = "Undefined"
 summaryObjMultiData[summaryObjMultiData$krel_FC_max < 1, ]$krelFCcategory = "krelFC_less_1"
 summaryObjMultiData[summaryObjMultiData$krel_FC_min >= 1, ]$krelFCcategory = "krelFC_more_1"
 summary_data <- summaryObjMultiData %>% 
        group_by(datasetName, krelFCcategory) %>%  
        summarize(numGenes = n(),
                  ) %>%
        ungroup() %>%
        group_by(datasetName) %>%
  mutate(totalGenes = sum(numGenes)) %>%
  mutate(proportion = numGenes / totalGenes) %>%
  select(-totalGenes) %>%
 ungroup() 
summary_data$percentLabel = "undefined"
#summary_data[ summary_data$krelFCcategory == "krelFC_less_1",]$percentLabel = sprintf("decrease = %0.1f%%", summary_data[ summary_data$krelFCcategory == "krelFC_less_1",]$proportion * 100) 
#summary_data[ summary_data$krelFCcategory == "krelFC_more_1",]$percentLabel = sprintf("increase = %0.1f%%", summary_data[ summary_data$krelFCcategory == "krelFC_more_1",]$proportion * 100) 
summary_data[ summary_data$krelFCcategory == "krelFC_less_1",]$percentLabel = sprintf("decrease = %0.0f%% (%d)", summary_data[ summary_data$krelFCcategory == "krelFC_less_1",]$proportion * 100, summary_data[ summary_data$krelFCcategory == "krelFC_less_1",]$numGenes) 
summary_data[ summary_data$krelFCcategory == "krelFC_more_1",]$percentLabel = sprintf("increase = %0.0f%% (%d)", summary_data[ summary_data$krelFCcategory == "krelFC_more_1",]$proportion * 100, summary_data[ summary_data$krelFCcategory == "krelFC_more_1",]$numGenes) 
  
  #print(summary_data)
   p<- ggplot(summaryObjMultiData, aes(x = datasetName, y = krel_FC_mean)) +
          geom_violin(scale="width", width=VIOLIN_WIDTH) + 
          geom_boxplot(width=0.25) + 
          geom_hline(yintercept = 1, linetype = "dashed", color = "#4d4d4d")
  if (showPoints) {
           p <- p + geom_point(position="jitter", alpha=pointALPHA)
  }
  DOWNCOORD = min(summaryObjMultiData$krel_FC_mean) * 0.85
  UPCOORD = max(summaryObjMultiData$krel_FC_mean) * 1.11
  
  ## add dataset proportion labels 
  datasets =  levels(summaryObjMultiData$datasetName)
   p <- p + scale_y_continuous(trans=log2_trans(), breaks=trans_breaks('log2', function(x) 2^x), labels=trans_format('log2', math_format(.x)), limits=c(NA,NA)) 
  if (showPercentLabels) {
          ind = 1
          for (dataset in datasets) {
            subset.summary = subset(summary_data, datasetName == dataset)
           print(subset.summary)
           #print(subset(subset.summary, krelFCcategory == "krelFC_more_1")$percentLabel)
           print(subset(subset.summary, krelFCcategory == "krelFC_more_1")$datasetName)
           downcoord = DOWNCOORD
           upcoord = UPCOORD
           if ((ind %% 2) != 1) {
              #downcoord = downcoord * 8/0.85
              #upcoord = upcoord * 0.3 / 1.11
              p <- p + geom_label(data=subset(subset.summary, krelFCcategory == "krelFC_less_1"), hjust=H_JUST, aes(label=percentLabel,x = datasetName, y = downcoord * 4/0.85),size=6,size.unit="mm")  
               p <- p + geom_label(data=subset(subset.summary, krelFCcategory == "krelFC_more_1"), hjust=H_JUST, aes(label=percentLabel,x = datasetName, y = upcoord * 0.3 / 1.11),size=6,size.unit="mm") 
            } else {
           print(downcoord)
           print(upcoord)
            #print(Y_DOWN_arr[ind])
           #print(Y_UP_arr[ind]) 
          p <- p + geom_label(data=subset(subset.summary, krelFCcategory == "krelFC_less_1"), hjust=H_JUST, aes(label=percentLabel,x = datasetName, y = downcoord),size=6,size.unit="mm")  
           p <- p + geom_label(data=subset(subset.summary, krelFCcategory == "krelFC_more_1"), hjust=H_JUST, aes(label=percentLabel,x = datasetName, y = upcoord),size=6,size.unit="mm") 
           }
          ind = ind + 1
          }   
  }
        p <- p + ylab(expression(log[2] ~ "Foldchange (FC) in" ~ italic(K[rel]))) +
        scale_x_discrete(labels = labelDatasets) + 
        theme_bw() + fontTheme +
        theme(legend.position="none",
              #plot.title = element_text(size=28, face=1.5,  hjust = 0.5),
              axis.text.x = element_text(size=26, face=1.5), #vjust=0.9, hjust=1),
              #axis.text.x = element_blank(),
              axis.title.x = element_blank(), 
              axis.title.y = element_text(size=32), 
              legend.title = element_blank())

     if (filetype == "png") {
           ggsave(paste0(filename,'.',filetype), plot=p, dpi=300, width=WIDTH)
        } else {
                ggsave(paste0(filename,'.',filetype), plot=p, width=WIDTH)
        }
 
}
krelFoldchangePlotter_twoDatasetsVertically <- function(summaryObj, filename, condnIdentifier, condnShortHand, constraintType="SteurerVock", filetype="png", H_JUST = 0.5, pointALPHA=0.2, UP=32, DOWN=1/4.5, WIDTH=6,showPoints = FALSE) {
 summaryObj$krelFCcategory = "Undefined"
summaryObj[summaryObj$krel_FC_max < 1, ]$krelFCcategory = "krelFC_less_1"
summaryObj[summaryObj$krel_FC_min >= 1, ]$krelFCcategory = "krelFC_more_1"
 summary_data <- summaryObj %>% 
        group_by(datasetName, krelFCcategory) %>%  
        summarize(numGenes = n(),
                  ) %>%
        ungroup() %>%
        group_by(datasetName) %>%
  mutate(totalGenes = sum(numGenes)) %>%
  mutate(proportion = numGenes / totalGenes) %>%
  select(-totalGenes) %>%
 ungroup() 
summary_data$percentLabel = "undefined"
#summary_data[ summary_data$krelFCcategory == "krelFC_less_1",]$percentLabel = sprintf("decrease = %0.1f%%", summary_data[ summary_data$krelFCcategory == "krelFC_less_1",]$proportion * 100) 
#summary_data[ summary_data$krelFCcategory == "krelFC_more_1",]$percentLabel = sprintf("increase = %0.1f%%", summary_data[ summary_data$krelFCcategory == "krelFC_more_1",]$proportion * 100) 
summary_data[ summary_data$krelFCcategory == "krelFC_less_1",]$percentLabel = sprintf("decrease = %0.0f%% (%d)", summary_data[ summary_data$krelFCcategory == "krelFC_less_1",]$proportion * 100, summary_data[ summary_data$krelFCcategory == "krelFC_less_1",]$numGenes) 
summary_data[ summary_data$krelFCcategory == "krelFC_more_1",]$percentLabel = sprintf("increase = %0.0f%% (%d)", summary_data[ summary_data$krelFCcategory == "krelFC_more_1",]$proportion * 100, summary_data[ summary_data$krelFCcategory == "krelFC_more_1",]$numGenes) 
 
## not pursued further
}

krelFoldchangePlotter <- function(summaryObj, filename, condnIdentifier, condnShortHand, constraintType="SteurerVock", filetype="png", H_JUST = 0.5, pointALPHA=0.2, UP=32, DOWN=1/4.5, WIDTH=6,showPoints = FALSE) {
 summaryObj$krelFCcategory = "Undefined"
summaryObj[summaryObj$krel_FC_max < 1, ]$krelFCcategory = "krelFC_less_1"
summaryObj[summaryObj$krel_FC_min >= 1, ]$krelFCcategory = "krelFC_more_1"

        df.krel.kinit = as.data.frame(cbind(c("krelFC_less_1", "krelFC_more_1"),
                             c(nrow(subset(summaryObj, krelFCcategory == "krelFC_less_1")) ,
                   nrow(subset(summaryObj, krelFCcategory == "krelFC_more_1")))))
  colnames(df.krel.kinit) = c("typeOfData", "countOfData")
  df.krel.kinit$countOfData = as.numeric(df.krel.kinit$countOfData)
  sum_all = sum(df.krel.kinit$countOfData)
 df.krel.kinit$countOfData = df.krel.kinit$countOfData
  df.krel.kinit$percentLabel = ""
  df.krel.kinit[df.krel.kinit$typeOfData == "krelFC_less_1", ]$percentLabel = sprintf("decrease = %0.1f%% (%d)",100.0*(df.krel.kinit[df.krel.kinit$typeOfData == "krelFC_less_1",]$countOfData/sum_all), df.krel.kinit[df.krel.kinit$typeOfData == "krelFC_less_1", ]$countOfData)
  df.krel.kinit[df.krel.kinit$typeOfData == "krelFC_more_1", ]$percentLabel = sprintf("increase = %0.1f%% (%d)",100.0*(df.krel.kinit[df.krel.kinit$typeOfData == "krelFC_more_1",]$countOfData/sum_all), df.krel.kinit[df.krel.kinit$typeOfData == "krelFC_more_1", ]$countOfData)
UP_COORD = max(summaryObj$krel_FC_min) + 2
DOWN_COORD = min(summaryObj$krel_FC_min) * 0.85
        p <- ggplot(summaryObj, aes(x= "krel", y=krel_FC_mean)) +
        geom_violin( size=1.2,width=0.3) +
        geom_boxplot( size=1.2,width=0.2)
     if (showPoints) {
        p <- p + geom_point( position="jitter", alpha=pointALPHA) 
     }
       p <- p + geom_label(data=subset(df.krel.kinit, typeOfData == "krelFC_less_1"), hjust=H_JUST, aes(label=percentLabel,x = "krel", y = DOWN_COORD),size=6,size.unit="mm") +
        geom_label(data=subset(df.krel.kinit, typeOfData == "krelFC_more_1"), hjust=H_JUST, aes(label=percentLabel,x = "krel", y = UP_COORD),size=6,size.unit="mm") +
        scale_y_continuous(trans=log2_trans(), breaks=trans_breaks('log2', function(x) 2^x), labels=trans_format('log2', math_format(.x))) +
       scale_x_discrete(labels = c(expression(atop("Pause release","("~italic(k[rel])~")" )) )) +
              ggtitle(condnIdentifier) +
       ylab(expression(log[2] ~ "Foldchange (FC) in" ~ italic(K[rel]))) +
        theme_bw() + fontTheme +
        theme(legend.position="none",
              plot.title = element_text(size=28, face=1.5,  hjust = 0.5),
              axis.text.x = element_text(size=26, face=1.5), #vjust=0.9, hjust=1),
              #axis.text.x = element_blank(),
              axis.title.x = element_blank(), axis.title.y = element_text(size=32), legend.title = element_blank())

         if (filetype == "png") {
           ggsave(paste0(filename,'.',filetype), plot=p, dpi=300, width=WIDTH)
        } else {
                ggsave(paste0(filename,'.',filetype), plot=p, width=WIDTH)
        }
}

kelongVsKrelPlotter <- function(summaryObj, filename, condnIdentifier, condnShortHand, constraintType="SteurerVock", filetype="png", krel_scale="min", pointALPHA=0.2) {
## absolute values of krel by using kelong values
condn.kelong.krel.df = summaryObj[,c("gene", "P_base_min", "B_base_min", "P_expr_min", "B_expr_min")]
condn.30.kelong.krel.df = condn.kelong.krel.df
condn.30.kelong.krel.df$elongation=30
condn.45.kelong.krel.df = condn.kelong.krel.df
condn.45.kelong.krel.df$elongation=45
condn.60.kelong.krel.df = condn.kelong.krel.df
condn.60.kelong.krel.df$elongation=60
condn.kelong.krel.df = rbind(condn.30.kelong.krel.df, condn.45.kelong.krel.df, condn.60.kelong.krel.df)
condn.kelong.krel.df$krel_base  = (condn.kelong.krel.df$elongation * condn.kelong.krel.df$B_base_min) / condn.kelong.krel.df$P_base_min
condn.kelong.krel.df$krel_expr  = (condn.kelong.krel.df$elongation * condn.kelong.krel.df$B_expr_min) / condn.kelong.krel.df$P_expr_min
condn.kelong.krel.df.baseline =  condn.kelong.krel.df[,c("gene", "elongation", "krel_base")]
condn.kelong.krel.df.baseline$cond = "control"
condn.kelong.krel.df.expr =  condn.kelong.krel.df[,c("gene", "elongation", "krel_expr")]
condn.kelong.krel.df.expr$cond = condnShortHand
colnames(condn.kelong.krel.df.baseline) = c("gene", "elongation", "krel", "cond")
colnames(condn.kelong.krel.df.expr) = c("gene", "elongation", "krel", "cond")
condn.kelong.krel.df = rbind(condn.kelong.krel.df.baseline, condn.kelong.krel.df.expr)
condn.kelong.krel.df$elongationType = sprintf("%d bp/sec", condn.kelong.krel.df$elongation)
condn.kelong.krel.df$cond = factor(condn.kelong.krel.df$cond, levels = c("control", condnShortHand))
if (krel_scale == "min") { # minute
    condn.kelong.krel.df$krel = condn.kelong.krel.df$krel * 60
}
p <- ggplot(condn.kelong.krel.df) +
        geom_violin(aes(x= elongationType, y= krel, color=cond), position="dodge", scale="width") +
        geom_boxplot(aes(x= elongationType, y= krel, color=cond), position=position_dodge(width = 0.9), alpha = 0.4, width=0.3) +
        scale_y_continuous(trans=log10_trans(),
          breaks=trans_breaks('log10', function(x) 10^x),
          labels=trans_format('log10', math_format(.x)))
        if (krel_scale == "min") {
            p <- p + ylab(expression(log[10] ~ italic(k[rel])~"values (minute"^"-1"~")"))
        }
        if (krel_scale == "sec") {
            p <- p + ylab(expression(log[10] ~ italic(k[rel])~"values (second"^"-1"~")"))
        }
 
           p <- p + xlab(expression(italic(k[elong]))) +
                  theme_bw() + fontTheme +
        theme(
              #legend.position="none",
              axis.text.x = element_text(size=28, face=1.5), #angle=30, vjust=0.9, hjust=1),
              #axis.text.x = element_blank(),
              #axis.title.x = element_blank(),
              axis.title.x = element_text(size=30),
              axis.title.y = element_text(size=32),
              legend.title = element_blank(),
              strip.text = element_text(size=25,face=1.2)
        )
 if (filetype == "png") {
           ggsave(paste0(filename,'.',filetype), plot=p, width=10, dpi=300)
        } else {
                ggsave(paste0(filename,'.',filetype), width=10, plot=p)
        }

}

kelongVsKrelPlotter_facet_dataset <- function(summaryObj, filename, condnShortHand = "treatment", constraintType="SteurerVock", filetype="pdf", krel_scale="min", pointALPHA=0.2, WIDTH=10, HEIGHT=8) {
## absolute values of krel by using kelong values
condn.kelong.krel.df = summaryObj[,c("datasetName","gene", "P_base_min", "B_base_min", "P_expr_min", "B_expr_min")]
condn.30.kelong.krel.df = condn.kelong.krel.df
condn.30.kelong.krel.df$elongation=30
condn.45.kelong.krel.df = condn.kelong.krel.df
condn.45.kelong.krel.df$elongation=45
condn.60.kelong.krel.df = condn.kelong.krel.df
condn.60.kelong.krel.df$elongation=60
condn.kelong.krel.df = rbind(condn.30.kelong.krel.df, condn.45.kelong.krel.df, condn.60.kelong.krel.df)
condn.kelong.krel.df$krel_base  = (condn.kelong.krel.df$elongation * condn.kelong.krel.df$B_base_min) / condn.kelong.krel.df$P_base_min
condn.kelong.krel.df$krel_expr  = (condn.kelong.krel.df$elongation * condn.kelong.krel.df$B_expr_min) / condn.kelong.krel.df$P_expr_min
condn.kelong.krel.df.baseline =  condn.kelong.krel.df[,c("datasetName","gene", "elongation", "krel_base")]
condn.kelong.krel.df.baseline$cond = "control"
condn.kelong.krel.df.expr =  condn.kelong.krel.df[,c("datasetName","gene", "elongation", "krel_expr")]
condn.kelong.krel.df.expr$cond = condnShortHand
colnames(condn.kelong.krel.df.baseline) = c("datasetName","gene", "elongation", "krel", "cond")
colnames(condn.kelong.krel.df.expr) = c("datasetName", "gene", "elongation", "krel", "cond")
condn.kelong.krel.df = rbind(condn.kelong.krel.df.baseline, condn.kelong.krel.df.expr)
condn.kelong.krel.df$elongationType = sprintf("%d bp/sec", condn.kelong.krel.df$elongation)
condn.kelong.krel.df$cond = factor(condn.kelong.krel.df$cond, levels = c("control", condnShortHand))
if (krel_scale == "min") { # minute
    condn.kelong.krel.df$krel = condn.kelong.krel.df$krel * 60
}
p <- ggplot(condn.kelong.krel.df) +
        geom_violin(aes(x= elongationType, y= krel, color=cond), position="dodge", scale="width") +
        geom_boxplot(aes(x= elongationType, y= krel, color=cond), position=position_dodge(width = 0.9), alpha = 0.4, width=0.3) +
        scale_y_continuous(trans=log10_trans(),
          breaks=trans_breaks('log10', function(x) 10^x),
          labels=trans_format('log10', math_format(.x))) +
                  facet_wrap(~datasetName, nrow = 2, scale = "free_y")
        if (krel_scale == "min") {
            p <- p + ylab(expression(log[10] ~ italic(k[rel])~"values (minute"^"-1"~")"))
        }
        if (krel_scale == "sec") {
            p <- p + ylab(expression(log[10] ~ italic(k[rel])~"values (second"^"-1"~")"))
        }
 
           p <- p + xlab(expression(italic(k[elong]))) +
                  theme_bw() + fontTheme +
        theme(
              #legend.position="none",
              axis.text.x = element_text(size=28, face=1.5), #angle=30, vjust=0.9, hjust=1),
              #axis.text.x = element_blank(),
              #axis.title.x = element_blank(),
              axis.title.x = element_text(size=30),
              axis.title.y = element_text(size=32),
              legend.title = element_blank(),
              legend.position="top",
              strip.text = element_text(size=25,face=1.2)
        )
 if (filetype == "png") {
           ggsave(paste0(filename,'.',filetype), plot=p, width= WIDTH, height=HEIGHT, dpi=300)
        } else {
                ggsave(paste0(filename,'.',filetype), width= WIDTH, height=HEIGHT, plot=p)
        }

}

kelongVsKrelPlotter_multiDatasets <- function(summaryObj, filename, labelDatasets, condnShortHand="treatment", kelong_bp_per_sec=30, filetype="png", krel_scale="min", WIDTH=22) {
## absolute values of krel by using kelong values
condn.kelong.krel.df = summaryObj[,c("datasetName","gene", "P_base_min", "B_base_min", "P_expr_min", "B_expr_min")]
condn.kelong.krel.df = condn.kelong.krel.df
condn.kelong.krel.df$elongation= kelong_bp_per_sec
condn.kelong.krel.df$krel_base  = (condn.kelong.krel.df$elongation * condn.kelong.krel.df$B_base_min) / condn.kelong.krel.df$P_base_min
condn.kelong.krel.df$krel_expr  = (condn.kelong.krel.df$elongation * condn.kelong.krel.df$B_expr_min) / condn.kelong.krel.df$P_expr_min
condn.kelong.krel.df.baseline =  condn.kelong.krel.df[,c("datasetName","gene", "elongation", "krel_base")]
condn.kelong.krel.df.baseline$cond = "control"
condn.kelong.krel.df.expr =  condn.kelong.krel.df[,c("datasetName","gene", "elongation", "krel_expr")]
condn.kelong.krel.df.expr$cond = condnShortHand
colnames(condn.kelong.krel.df.baseline) = c("datasetName","gene", "elongation", "krel", "cond")
colnames(condn.kelong.krel.df.expr) = c("datasetName","gene", "elongation", "krel", "cond")
condn.kelong.krel.df = rbind(condn.kelong.krel.df.baseline, condn.kelong.krel.df.expr)
condn.kelong.krel.df$elongationType = sprintf("%d bp/sec", condn.kelong.krel.df$elongation)
condn.kelong.krel.df$cond = factor(condn.kelong.krel.df$cond, levels = c("control", condnShortHand))
if (krel_scale == "min") { # minute
    condn.kelong.krel.df$krel = condn.kelong.krel.df$krel * 60
}
p <- ggplot(condn.kelong.krel.df, aes(x=datasetName, y = krel,color=cond)) +
        geom_violin(position="dodge", scale="width") +
        geom_boxplot(position=position_dodge(width = 0.9), alpha = 0.4, width=0.3) +
        scale_y_continuous(trans=log10_trans(),
          breaks=trans_breaks('log10', function(x) 10^x),
          labels=trans_format('log10', math_format(.x)))
        if (krel_scale == "min") {
            p <- p + ylab(expression(log[10] ~ italic(k[rel])~"values (minute"^"-1"~")"))
        }
        if (krel_scale == "sec") {
            p <- p + ylab(expression(log[10] ~ italic(k[rel])~"values (second"^"-1"~")"))
        }
 
       #    p <- p + xlab(expression(italic(k[elong]))) +
       p <- p + scale_x_discrete(labels = labelDatasets) + 
                  theme_bw() + fontTheme +
        ggtitle(bquote(italic(k[elong])==.(kelong_bp_per_sec * 60)~"bp / min")) +
             theme(
              #legend.position="none",
              axis.text.x = element_text(size=28, face=1.5), #angle=30, vjust=0.9, hjust=1),
              #axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              #axis.title.x = element_text(size=30),
              axis.title.y = element_text(size=32),
              legend.title = element_blank(),
              plot.title = element_text(size=26,face=1.2, hjust=0.5)
        )
 if (filetype == "png") {
           ggsave(paste0(filename,'.',filetype), plot=p, width=WIDTH, dpi=300)
        } else {
                ggsave(paste0(filename,'.',filetype), width=WIDTH, plot=p)
        }

}

kelongVsKrelPlotter_multiDatasets_twoPanels <- function(summaryObj, filename, condnShortHand="treatment", kelong_bp_per_min=2000, filetype="pdf", krel_scale="min", HEIGHT=8, WIDTH=22) {
## absolute values of krel by using kelong values
condn.kelong.krel.df = summaryObj[,c("datasetName","gene", "P_base_min", "B_base_min", "P_expr_min", "B_expr_min")]
condn.kelong.krel.df = condn.kelong.krel.df
condn.kelong.krel.df$elongation= kelong_bp_per_min / 60
condn.kelong.krel.df$krel_base  = (condn.kelong.krel.df$elongation * condn.kelong.krel.df$B_base_min) / condn.kelong.krel.df$P_base_min
condn.kelong.krel.df$krel_expr  = (condn.kelong.krel.df$elongation * condn.kelong.krel.df$B_expr_min) / condn.kelong.krel.df$P_expr_min
condn.kelong.krel.df.baseline =  condn.kelong.krel.df[,c("datasetName","gene", "elongation", "krel_base")]
condn.kelong.krel.df.baseline$cond = "control"
condn.kelong.krel.df.expr =  condn.kelong.krel.df[,c("datasetName","gene", "elongation", "krel_expr")]
condn.kelong.krel.df.expr$cond = condnShortHand
colnames(condn.kelong.krel.df.baseline) = c("datasetName","gene", "elongation", "krel", "cond")
colnames(condn.kelong.krel.df.expr) = c("datasetName","gene", "elongation", "krel", "cond")
condn.kelong.krel.df = rbind(condn.kelong.krel.df.baseline, condn.kelong.krel.df.expr)
#condn.kelong.krel.df$elongationType = sprintf("%d bp/sec", condn.kelong.krel.df$elongation)
condn.kelong.krel.df$cond = factor(condn.kelong.krel.df$cond, levels = c("control", condnShortHand))
if (krel_scale == "min") { # minute
    condn.kelong.krel.df$krel = condn.kelong.krel.df$krel * 60
}
p <- ggplot(condn.kelong.krel.df, aes(x=cond, y = krel,color=cond)) +
        geom_violin(position="dodge", scale="width") +
        geom_boxplot(position=position_dodge(width = 0.9), alpha = 0.4, width=0.3) +
        scale_y_continuous(trans=log10_trans(),
          breaks=trans_breaks('log10', function(x) 10^x),
          labels=trans_format('log10', math_format(.x)))
        if (krel_scale == "min") {
            p <- p + ylab(expression(log[10] ~ italic(k[rel])~"values (minute"^"-1"~")"))
        }
        if (krel_scale == "sec") {
            p <- p + ylab(expression(log[10] ~ italic(k[rel])~"values (second"^"-1"~")"))
        }

        p <- p  + facet_wrap(~ datasetName, nrow = 1)
       #    p <- p + xlab(expression(italic(k[elong]))) +
      # p <- p + scale_x_discrete(labels = labelDatasets) + 
        p <- p + theme_bw() + fontTheme +
        ggtitle(bquote(italic(k[elong])==.(kelong_bp_per_min)~"bp / min")) +
             theme(
              legend.position="none",
              axis.text.x = element_text(size=28, face=1.5), #angle=30, vjust=0.9, hjust=1),
              #axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              #axis.title.x = element_text(size=30),
              axis.title.y = element_text(size=32),
              legend.title = element_blank(),
             # strip.text = element_text(size = 24),
             strip.background = element_rect(fill = "#E6E6FA", color = "black", linewidth = 1.6),
    strip.text = element_text(color = "#6A5ACD", face = "bold", size = 24), 
              plot.title = element_text(size=26,face=1.2, hjust=0.5)
        )
 if (filetype == "png") {
           ggsave(paste0(filename,'.',filetype), plot=p, width=WIDTH, height=HEIGHT, dpi=300)
        } else {
                ggsave(paste0(filename,'.',filetype), width=WIDTH, height=HEIGHT, plot=p)
        }
  return(condn.kelong.krel.df)
}

kelongVsEffectiveKrelPlotter_multiDatasets_twoPanels <- function(summaryObj, filename, condnShortHand="treatment", kelong_bp_per_min=2000, filetype="pdf", krel_scale="min", WIDTH=22, HEIGHT=10) {
## absolute values of krel by using kelong values
condn.kelong.krel.df = summaryObj[,c("datasetName","gene", "P_base_min", "B_base_min", "P_expr_min", "B_expr_min")]
condn.kelong.krel.df = condn.kelong.krel.df
condn.kelong.krel.df$elongation= kelong_bp_per_min / 60
condn.kelong.krel.df$krel_base  = (condn.kelong.krel.df$elongation * condn.kelong.krel.df$B_base_min) / condn.kelong.krel.df$P_base_min
condn.kelong.krel.df$krel_expr  = (condn.kelong.krel.df$elongation * condn.kelong.krel.df$B_expr_min) / condn.kelong.krel.df$P_expr_min
condn.kelong.krel.df.baseline =  condn.kelong.krel.df[,c("datasetName","gene", "elongation", "krel_base", "P_base_min")]
condn.kelong.krel.df.baseline$cond = "control"
condn.kelong.krel.df.expr =  condn.kelong.krel.df[,c("datasetName","gene", "elongation", "krel_expr", "P_expr_min")]
condn.kelong.krel.df.expr$cond = condnShortHand
colnames(condn.kelong.krel.df.baseline) = c("datasetName","gene", "elongation", "krel", "P", "cond")
colnames(condn.kelong.krel.df.expr) = c("datasetName","gene", "elongation", "krel", "P", "cond")
condn.kelong.krel.df = rbind(condn.kelong.krel.df.baseline, condn.kelong.krel.df.expr)
#condn.kelong.krel.df$elongationType = sprintf("%d bp/sec", condn.kelong.krel.df$elongation)
condn.kelong.krel.df$cond = factor(condn.kelong.krel.df$cond, levels = c("control", condnShortHand))
if (krel_scale == "min") { # minute
    condn.kelong.krel.df$krel = condn.kelong.krel.df$krel * 60
}
condn.kelong.krel.df$effectiveKrel =  condn.kelong.krel.df$krel * condn.kelong.krel.df$P 
p <- ggplot(condn.kelong.krel.df, aes(x=cond, y = effectiveKrel,color=cond)) +
        geom_violin(position="dodge", scale="width") +
        geom_boxplot(position=position_dodge(width = 0.9), alpha = 0.4, width=0.3) +
        scale_y_continuous(trans=log10_trans(),
          breaks=trans_breaks('log10', function(x) 10^x),
          labels=trans_format('log10', math_format(.x)))
        if (krel_scale == "min") {
            p <- p + ylab(expression(log[10] ~"Effective"~italic(k[rel])~"(RNAP x minute"^"-1"~")"))
        }
        if (krel_scale == "sec") {
            p <- p + ylab(expression(log[10] ~"Effective"~italic(k[rel])~"(RNAP x second"^"-1"~")"))
        }

        p <- p  + facet_wrap(~ datasetName, nrow = 1)
       #    p <- p + xlab(expression(italic(k[elong]))) +
      # p <- p + scale_x_discrete(labels = labelDatasets) + 
        p <- p + theme_bw() + fontTheme +
        ggtitle(bquote(italic(k[elong])==.(kelong_bp_per_min)~"bp / min")) +
             theme(
              legend.position="none",
              axis.text.x = element_text(size=28, face=1.5), #angle=30, vjust=0.9, hjust=1),
              #axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              #axis.title.x = element_text(size=30),
              axis.title.y = element_text(size=32),
              legend.title = element_blank(),
             # strip.text = element_text(size = 24),
             strip.background = element_rect(fill = "#E6E6FA", color = "black", linewidth = 1.6),
    strip.text = element_text(color = "#6A5ACD", face = "bold", size = 24), 
              plot.title = element_text(size=26,face=1.2, hjust=0.5)
        )
 if (filetype == "png") {
           ggsave(paste0(filename,'.',filetype), plot=p, height=HEIGHT, width=WIDTH, dpi=300)
        } else {
                ggsave(paste0(filename,'.',filetype), height=HEIGHT, width=WIDTH, plot=p)
        }
  return(condn.kelong.krel.df)
}

kinitFC_barplots_multiDatasets <- function(summaryObjMultiData, filename, labelDatasets,filetype="pdf", WIDTH=8, HEIGHT=9,bar.colors = c("#2ca02c", "#d62728"), type.colors =c("#36454F", "#D2B48C") ,REVERSE=TRUE, addTotalLabel = TRUE, VJUST = 0.5, HJUST=0.8,VJUST_X = 0.9, textsize=10) {
        condn.kinitFC.krel = summaryObjMultiObj[,c("datasetName", "gene", "P_base_min", "P_expr_min", "krel_FC_mean")]
        condn.kinitFC.krel$P_foldchange = condn.kinitFC.krel$P_expr_min /  condn.kinitFC.krel$P_base_min
        condn.kinitFC.krel$P_FC_times_krelFC = condn.kinitFC.krel$P_foldchange * condn.kinitFC.krel$krel_FC_mean
        condn.kinitFC.krel.P_FC = condn.kinitFC.krel[,c("datasetName", "gene", "P_foldchange")]
        condn.kinitFC.krel.P_FC_times_krelFC = condn.kinitFC.krel[,c("datasetName", "gene", "P_FC_times_krelFC")]
        condn.kinitFC.krel.P_FC$type = "P_foldchange"
        condn.kinitFC.krel.P_FC_times_krelFC$type = "P_FC_timeskrelFC"
        colnames(condn.kinitFC.krel.P_FC) = c("datasetName","gene", "kinit_FC", "type")
        colnames(condn.kinitFC.krel.P_FC_times_krelFC) = c("datasetName","gene", "kinit_FC", "type")
        condn.kinitFC.krel = rbind(condn.kinitFC.krel.P_FC, condn.kinitFC.krel.P_FC_times_krelFC)
        condn.kinitFC.krel$type = factor(condn.kinitFC.krel$type, levels = c("P_foldchange", "P_FC_timeskrelFC"))

         condn.kinitFC.krel$kinitFCcategory = "Undefined"
         condn.kinitFC.krel[condn.kinitFC.krel$kinit_FC < 1, ]$kinitFCcategory = "kinitFC_less_1"
         condn.kinitFC.krel[condn.kinitFC.krel$kinit_FC >= 1, ]$kinitFCcategory = "kinitFC_more_1"

        summary_data <- condn.kinitFC.krel %>%
                group_by(datasetName, type, kinitFCcategory) %>%
                summarize(numGenes = n()) %>%
              mutate(total = sum(numGenes),
                proportion = numGenes / total) %>%
                #select(-total) %>%
                ungroup()
       print(summary_data)
       write.table(summary_data, sprintf("%s.txt", filename), row.names=F, sep="\t", quote=F)
        custom_labels <- c(
     "P_foldchange" = bquote(italic(k[pre])~">>"~italic(k[rel])),
     "P_FC_timeskrelFC" = bquote(italic(k[pre])~"<<"~italic(k[rel]))
     )
   summary_data$krelFCcategory = factor(summary_data$kinitFCcategory , levels = c("kinitFC_less_1", "kinitFC_more_1"))
   p <- ggplot(data=summary_data, aes(x = datasetName, y = proportion, color=type, fill = krelFCcategory, color=type)) +
           geom_bar(stat = "identity")#  position = position_stack(reverse=REVERSE))
  #     geom_bar(data=subset(summary_data, type == "P_foldchange"), stat = "identity", position = position_stack(reverse=REVERSE), just=0)
 #p <- p + geom_bar(data=subset(summary_data, type == "P_FC_timeskrelFC"), stat = "identity", position = position_stack(reverse=REVERSE), just=1)
      p <- p +   scale_fill_manual(labels = c(expression("Increase" ~ italic(k[init])), expression("Decrease" ~ italic(k[init]))), values = bar.colors) +
         #facet_wrap(~type, labeller = as_labeller(custom_labels))
         #scale_color_manual(labels = c(expression(italic(k[pre])~">>"~italic(k[rel])), expression(italic(k[pre])~"<<"~italic(k[rel]))), values= type.colors)
  labs(y = "Proportion of genes")
  if (addTotalLabel) {
          datasets = levels(summaryObjMultiData$datasetName)
          for (dataset in datasets) {
                  totalgenes = subset(summary_data, datasetName == dataset)$total[1]
                  p <- p +  geom_text(x=dataset, y = 1, vjust = VJUST, position = "identity", family = "sans",label = totalgenes, size=textsize)
          }
          }   
       p <- p + scale_x_discrete(labels = labelDatasets) + 
        theme_bw() + fontTheme +
        theme(
              #legend.position="none",
              #plot.title = element_text(size=28, face=1.5,  hjust = 0.5),
              axis.text.x = element_text(size=26, face=1.2, angle = 30, hjust=HJUST, vjust=VJUST_X), #vjust=0.9, hjust=1),
           strip.background = element_rect(fill = "#E6E6FA", color = "black", linewidth = 1.6),
    strip.text = element_text(color = "#6A5ACD", face = "bold", size = 24), 
              #axis.text.x = element_blank(),
              axis.title.x = element_blank(), 
              axis.title.y = element_text(size=32), 
              legend.title = element_blank())

     if (filetype == "png") {
           ggsave(paste0(filename,'.',filetype), plot=p, dpi=300, width=WIDTH, height= HEIGHT)
        } else {
                ggsave(paste0(filename,'.',filetype), plot=p, width=WIDTH, height= HEIGHT)
        }
 
}

# summaryObjMultiObj has a `datasetName` that is ordered
kinitFCvsPfoldchangePlotter_multidatasets <- function(summaryObjMultiObj, filename, labelDatasets, Y_DOWN_arr = c(1/4.5,1/4.5), Y_UP_arr = c(4, 4, 4), filetype="png", H_JUST_LEFT = 0.5, H_JUST_RIGHT=0.5, diff_UP = -3, diff_DOWN= 0.2, pointALPHA=0.2, showPercentLabels = FALSE, showPoints = FALSE, WIDTH=16, legend.colors = c("#B3CDE3", "#6497B1")) {
## kinit FC extremes in violin
        condn.kinitFC.krel = summaryObjMultiObj[,c("datasetName", "gene", "P_base_min", "P_expr_min", "krel_FC_mean")]
        condn.kinitFC.krel$P_foldchange = condn.kinitFC.krel$P_expr_min /  condn.kinitFC.krel$P_base_min
        condn.kinitFC.krel$P_FC_times_krelFC = condn.kinitFC.krel$P_foldchange * condn.kinitFC.krel$krel_FC_mean
        condn.kinitFC.krel.P_FC = condn.kinitFC.krel[,c("datasetName", "gene", "P_foldchange")]
        condn.kinitFC.krel.P_FC_times_krelFC = condn.kinitFC.krel[,c("datasetName", "gene", "P_FC_times_krelFC")]
        condn.kinitFC.krel.P_FC$type = "P_foldchange"
        condn.kinitFC.krel.P_FC_times_krelFC$type = "P_FC_timeskrelFC"
        colnames(condn.kinitFC.krel.P_FC) = c("datasetName","gene", "kinit_FC", "type")
        colnames(condn.kinitFC.krel.P_FC_times_krelFC) = c("datasetName","gene", "kinit_FC", "type")
        condn.kinitFC.krel = rbind(condn.kinitFC.krel.P_FC, condn.kinitFC.krel.P_FC_times_krelFC)
        condn.kinitFC.krel$type = factor(condn.kinitFC.krel$type, levels = c("P_foldchange", "P_FC_timeskrelFC"))

         condn.kinitFC.krel$kinitFCcategory = "Undefined"
         condn.kinitFC.krel[condn.kinitFC.krel$kinit_FC < 1, ]$kinitFCcategory = "kinitFC_less_1"
         condn.kinitFC.krel[condn.kinitFC.krel$kinit_FC >= 1, ]$kinitFCcategory = "kinitFC_more_1"
         
        summary_data <- condn.kinitFC.krel %>% 
                group_by(datasetName, type, kinitFCcategory) %>%  
                summarize(numGenes = n()) %>%
              mutate(total = sum(numGenes),
                proportion = numGenes / total) %>%
                #select(-total) %>%
                ungroup()
       print(summary_data) 
      summary_data$percentLabel = "undefined"
      summary_data[ summary_data$kinitFCcategory == "kinitFC_less_1",]$percentLabel = sprintf("decrease = %0.0f%%", summary_data[ summary_data$kinitFCcategory == "kinitFC_less_1",]$proportion * 100) 
      summary_data[ summary_data$kinitFCcategory == "kinitFC_more_1",]$percentLabel = sprintf("increase = %0.0f%%", summary_data[ summary_data$kinitFCcategory == "kinitFC_more_1",]$proportion * 100)  
      #summary_data[ summary_data$kinitFCcategory == "kinitFC_less_1",]$percentLabel = sprintf("decrease = %0.1f%% (%d)", summary_data[ summary_data$kinitFCcategory == "kinitFC_less_1",]$proportion * 100, summary_data[ summary_data$kinitFCcategory == "kinitFC_less_1",]$numGenes) 
      #summary_data[ summary_data$kinitFCcategory == "kinitFC_more_1",]$percentLabel = sprintf("increase = %0.1f%% (%d)", summary_data[ summary_data$kinitFCcategory == "kinitFC_more_1",]$proportion * 100, summary_data[ summary_data$kinitFCcategory == "kinitFC_more_1",]$numGenes)
 
      #print(summary_data)
   p <- ggplot(condn.kinitFC.krel, aes(x = datasetName, y = kinit_FC, color = type)) +
        geom_violin(position="dodge",scale="width") + 
        geom_boxplot(position = position_dodge(width=0.9), width=0.3) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "#4d4d4d")
   if (showPoints) {
       p <- p + geom_point(position="jitter", alpha=pointALPHA)
    }
  p <- p + scale_y_continuous(trans=log2_trans(), breaks=trans_breaks('log2', function(x) 2^x), labels=trans_format('log2', math_format(.x))) 
  ## add dataset proportion labels 
  datasets =  levels(summary_data$datasetName)
     kinit_FC_more_limit = max(condn.kinitFC.krel$kinit_FC)*1.5 
     kinit_FC_less_limit = min(condn.kinitFC.krel$kinit_FC) * 1.15
  
if (showPercentLabels) {
  ind = 1
  for (dataset in datasets) {
    subset.summary = subset(summary_data, datasetName == dataset)
   print(subset.summary)
   print(subset(subset.summary, kinitFCcategory == "kinitFC_more_1")$percentLabel)
   print(subset(subset.summary, kinitFCcategory == "kinitFC_more_1")$datasetName)
   print(Y_UP_arr[ind]) 
  print(Y_DOWN_arr[ind] +diff_DOWN)
   print(Y_UP_arr[ind] + diff_UP) 
#  p <- p + geom_label(data=subset(subset.summary, kinitFCcategory == "kinitFC_less_1" & type == "P_foldchange"), hjust= H_JUST_LEFT, aes(label=percentLabel,x = datasetName, y = Y_DOWN_arr[ind]),size=6,size.unit="mm")  
#   p <- p + geom_label(data=subset(subset.summary, kinitFCcategory == "kinitFC_more_1" & type == "P_foldchange" ), hjust=H_JUST_LEFT, aes(label=percentLabel,x = datasetName, y = Y_UP_arr[ind]),size=6,size.unit="mm") 
# p <- p + geom_label(data=subset(subset.summary, kinitFCcategory == "kinitFC_less_1" & type == "P_FC_timeskrelFC"), hjust= H_JUST_RIGHT, aes(label=percentLabel,x = datasetName, y = Y_DOWN_arr[ind]+ diff_DOWN),size=6,size.unit="mm")  
#   p <- p + geom_label(data=subset(subset.summary, kinitFCcategory == "kinitFC_more_1" & type == "P_FC_timeskrelFC" ), hjust=H_JUST_RIGHT, aes(label=percentLabel,x = datasetName, y = Y_UP_arr[ind]+diff_UP),size=6,size.unit="mm") 
  p <- p + geom_label(data=subset(subset.summary, kinitFCcategory == "kinitFC_less_1" & type == "P_foldchange"), hjust= H_JUST_LEFT, aes(label=percentLabel,x = datasetName, y = kinit_FC_less_limit),size=6,size.unit="mm")  
   p <- p + geom_label(data=subset(subset.summary, kinitFCcategory == "kinitFC_more_1" & type == "P_foldchange" ), hjust=H_JUST_LEFT, aes(label=percentLabel,x = datasetName, y = kinit_FC_more_limit/3),size=6,size.unit="mm") 
 p <- p + geom_label(data=subset(subset.summary, kinitFCcategory == "kinitFC_less_1" & type == "P_FC_timeskrelFC"), hjust= H_JUST_RIGHT, aes(label=percentLabel,x = datasetName, y = kinit_FC_less_limit*3.15),size=6,size.unit="mm")  
   p <- p + geom_label(data=subset(subset.summary, kinitFCcategory == "kinitFC_more_1" & type == "P_FC_timeskrelFC" ), hjust=H_JUST_RIGHT, aes(label=percentLabel,x = datasetName, y = kinit_FC_more_limit/6.1),size=6,size.unit="mm") 
  
 ind = ind + 1
  }   
}
       #p <- p + annotate("text",hjust = 0.7, label = "FC (Pause Sum)", x = datasets[1], y = kinit_FC_more_limit*0.7, color=legend.colors[1], size=8)      
       #p <- p + annotate("text", hjust = 0.45, label = expression("FC(Pause Sum) x FC("~italic(k[rel])~")"), x = datasets[1], y= kinit_FC_more_limit*0.5, color= legend.colors[2], size=8) 
     if (showPercentLabels) { 
         p <- p + annotate("text",hjust = 0.7, label = expression(italic(k[pre])~">>"~italic(k[rel])), x = datasets[1], y = kinit_FC_more_limit*2.2, color=legend.colors[1], size=11, family="sans",fontface = "bold")      
         p <- p + annotate("text", hjust = 0.7, label = expression(italic(k[pre])~"<<"~italic(k[rel])), x = datasets[1], y= kinit_FC_more_limit, color= legend.colors[2], size=11, family="sans",fontface = "bold") 
      } else {
         p <- p + annotate("text",hjust = 0.4, label = expression(italic(k[pre])~">>"~italic(k[rel])), x = datasets[1], y = kinit_FC_more_limit/4.5, color=legend.colors[1], size=14, family="sans",fontface = "bold")      
         p <- p + annotate("text", hjust = 0.4, label = expression(italic(k[pre])~"<<"~italic(k[rel])), x = datasets[1], y= kinit_FC_more_limit/15.2, color= legend.colors[2], size=14, family="sans",fontface = "bold") 
    
      
     } 
       p <- p + scale_color_manual(values = legend.colors) +
                      # labels = c("FC (Pause Sum)",
                      #              expression(atop("FC(Pause Sum) x", "Foldchange("~italic(k[rel])~")")))) + 
        ylab(expression(log[2] ~ "Foldchange (FC) in" ~ italic(K[init]))) +
        scale_x_discrete(labels = labelDatasets) + 
        theme_bw() + fontTheme +
        theme(
              legend.position="none",
              #plot.title = element_text(size=28, face=1.5,  hjust = 0.5),
              axis.text.x = element_text(size=26, face=1.5), #vjust=0.9, hjust=1),
              #axis.text.x = element_blank(),
              axis.title.x = element_blank(), 
              axis.title.y = element_text(size=32), 
              legend.title = element_blank())

     if (filetype == "png") {
           ggsave(paste0(filename,'.',filetype), plot=p, dpi=300, width=WIDTH)
        } else {
                ggsave(paste0(filename,'.',filetype), plot=p, width=WIDTH)
        }
 

}

kinitFCvsPfoldchangePlotter <- function(summaryObj, filename, condnIdentifier, condnShortHand, directionKinit="up", filetype="png", pointALPHA=0.2) {
## kinit FC extremes in violin
condn.kinitFC.krel = summaryObj[,c("gene", "P_base_min", "P_expr_min", "krel_FC_mean")]
condn.kinitFC.krel$P_foldchange = condn.kinitFC.krel$P_expr_min /  condn.kinitFC.krel$P_base_min
condn.kinitFC.krel$P_FC_times_krelFC = condn.kinitFC.krel$P_foldchange * condn.kinitFC.krel$krel_FC_mean
condn.kinitFC.krel.P_FC = condn.kinitFC.krel[,c("gene", "P_foldchange")]
condn.kinitFC.krel.P_FC_times_krelFC = condn.kinitFC.krel[,c("gene", "P_FC_times_krelFC")]
condn.kinitFC.krel.P_FC$type = "P_foldchange"
condn.kinitFC.krel.P_FC_times_krelFC$type = "P_FC_timeskrelFC"
colnames(condn.kinitFC.krel.P_FC) = c("gene", "kinit_FC", "type")
colnames(condn.kinitFC.krel.P_FC_times_krelFC) = c("gene", "kinit_FC", "type")
condn.kinitFC.krel = rbind(condn.kinitFC.krel.P_FC, condn.kinitFC.krel.P_FC_times_krelFC)
condn.kinitFC.krel$type = factor(condn.kinitFC.krel$type, levels = c("P_foldchange", "P_FC_timeskrelFC"))

# frequency of decreasing kinit
if (directionKinit == "down") {
n.p_fc.kinitFCdirection = nrow(subset(condn.kinitFC.krel.P_FC, kinit_FC < 1))
n.p_fctimeskrelFC.kinitFCdirection = nrow(subset(condn.kinitFC.krel.P_FC_times_krelFC, kinit_FC < 1))
} else {
n.p_fc.kinitFCdirection = nrow(subset(condn.kinitFC.krel.P_FC, kinit_FC > 1))
n.p_fctimeskrelFC.kinitFCdirection = nrow(subset(condn.kinitFC.krel.P_FC_times_krelFC, kinit_FC > 1))
}
condn.kinitFC.df = c()
condn.kinitFC.df$direction_prop = n.p_fc.kinitFCdirection / nrow(condn.kinitFC.krel.P_FC)
condn.kinitFC.df$type = "P_foldchange"
condn.kinitFC.df = as.data.frame(condn.kinitFC.df)
condn.kinitFC.df$direction_prop = as.numeric(condn.kinitFC.df$direction_prop)

if (directionKinit == "down") {
   condn.kinitFC.df$textLabel = sprintf("decrease = %0.1f%% (%d)", 100 * condn.kinitFC.df$direction_prop, n.p_fc.kinitFCdirection)
} else {
   condn.kinitFC.df$textLabel = sprintf("increase = %0.1f%% (%d)", 100 * condn.kinitFC.df$direction_prop, n.p_fc.kinitFCdirection)
}

condn.kinitFC.krelFC.df = c()
condn.kinitFC.krelFC.df$direction_prop = n.p_fctimeskrelFC.kinitFCdirection / nrow(condn.kinitFC.krel.P_FC_times_krelFC)
condn.kinitFC.krelFC.df$type = "P_FC_timeskrelFC"
condn.kinitFC.krelFC.df = as.data.frame(condn.kinitFC.krelFC.df)
condn.kinitFC.krelFC.df$direction_prop = as.numeric(condn.kinitFC.krelFC.df$direction_prop)
if (directionKinit == "down") {
   condn.kinitFC.krelFC.df$textLabel = sprintf("decrease = %0.1f%% (%d)", 100 * condn.kinitFC.krelFC.df$direction_prop, n.p_fctimeskrelFC.kinitFCdirection)
} else {
   condn.kinitFC.krelFC.df$textLabel = sprintf("increase = %0.1f%% (%d)", 100 * condn.kinitFC.krelFC.df$direction_prop, n.p_fctimeskrelFC.kinitFCdirection)
}
condn.kinitFC.krelFC.df = rbind(condn.kinitFC.krelFC.df, condn.kinitFC.df)
Y_COORD = max(condn.kinitFC.krel$kinit_FC) + 2
## count frequency of pfoldchange and P_FC times krel_FC
p<- ggplot(condn.kinitFC.krel, aes(x = type, y = kinit_FC)) + 
        geom_violin( width=1, size=1.2) +
        geom_boxplot( alpha=0.4,width=0.2, size=1.2) +
        geom_point( position="jitter",alpha = pointALPHA) +
        #geom_text(data=Dex100nm.kinitFC.krelFC.df, aes(label=textLabel,x = type, y = 8, size=2, size.unit="cm")) +
        geom_label(data=condn.kinitFC.krelFC.df, aes(label=textLabel,x = type, y = Y_COORD),size=6,size.unit="mm") +
        scale_y_continuous(trans=log2_trans(),
          breaks=trans_breaks('log2', function(x) 2^x),
          labels=trans_format('log2', math_format(.x))) +
                  ylab(expression(log[2] ~ "Foldchange (FC) in"~italic(k[init]))) +
          scale_x_discrete(labels = c(expression(italic(k[pre])~">>"~italic(k[rel])),
                                    expression(italic(k[pre])~"<<"~italic(k[rel])))) +
        #scale_x_discrete(labels = c(expression(atop("FC","(Pause Sum)")),
         #                           expression(atop("FC(Pause Sum) x","Foldchange("~italic(k[rel])~")")))) +
               ggtitle(paste(condnIdentifier, condnShortHand)) +
                  theme_bw() + fontTheme +
        theme(
              plot.title = element_text(size=28, face=1.5,  hjust = 0.5),
              legend.position="none",
              axis.text.x = element_text(size=28, face=1.5), #angle=30, vjust=0.9, hjust=1),
              #axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              #axis.title.x = element_text(size=30),
              axis.title.y = element_text(size=32),
              legend.title = element_blank()
        )
 if (filetype == "png") {
           ggsave(paste0(filename,'.',filetype), plot=p, width=8, dpi=300)
        } else {
                ggsave(paste0(filename,'.',filetype), width=8, plot=p)
        }

}

