@
library(ggplot2)
library(dplyr)
library(ggsci)
library(RColorBrewer)
library(scales)
# https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html
source("plotting_composites_lattice.R")

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

#scanDataDir = "/Users/rudradeepmukherjee/Documents/UConn Health/comparment.model.project/datasets, genome/Jonkers_2014_Trp_FP_paper/ParamEstimation_FP_TRP/InternalControl_TRP12min_repressed_KPRE_CONSTRAINED_ALLPARAMS/figures"
scanDataDir = "InternalControl_TRP12min_repressed_KpreConstrained/"
rateData.Kpre.constrained = readRDS(sprintf("%s/repeatedParamScanData_TRP_12min_repressed_Steurer_constrained.rds", scanDataDir))
rateData.Kpre.constrained$P_base = rateData.Kpre.constrained$kinit_base / (rateData.Kpre.constrained$kpre + rateData.Kpre.constrained$krel_base)
rateData.Kpre.constrained$P_expr = rateData.Kpre.constrained$kinit_expr / (rateData.Kpre.constrained$kpre + rateData.Kpre.constrained$krel_expr)
rateData.Kpre.constrained$B_base = ( rateData.Kpre.constrained$krel_base * rateData.Kpre.constrained$P_base) / (rateData.Kpre.constrained$kelong)
rateData.Kpre.constrained$B_expr = ( rateData.Kpre.constrained$krel_expr * rateData.Kpre.constrained$P_expr) / (rateData.Kpre.constrained$kelong)
rateData.Kpre.constrained$Release_base = rateData.Kpre.constrained$P_base * rateData.Kpre.constrained$krel_base
rateData.Kpre.constrained$Release_expr = rateData.Kpre.constrained$P_expr * rateData.Kpre.constrained$krel_expr
rateData.Kpre.constrained$Terminal_base = rateData.Kpre.constrained$kpre * rateData.Kpre.constrained$P_base
rateData.Kpre.constrained$Terminal_expr = rateData.Kpre.constrained$kpre * rateData.Kpre.constrained$P_expr

newSummary.TRP.Kpre.constrained = group_by(rateData.Kpre.constrained, gene) %>%  summarise_at(c("Terminal_base", "Terminal_expr", "Release_base", "Release_expr", "kinit_base", "kinit_expr", "krel_base", "krel_expr", "kpre", "kelong", "kinit_FC", "krel_FC", "P_base", "P_expr", "B_base", "B_expr"), c("min", "max", "mean", "var"))

newSummary.TRP.Kpre.constrained$IndexDisp = newSummary.TRP.Kpre.constrained$kinit_FC_var /  newSummary.TRP.Kpre.constrained$kinit_FC_mean
newSummary.TRP.Kpre.constrained$IndexDisp_Krel = newSummary.TRP.Kpre.constrained$krel_FC_var /  newSummary.TRP.Kpre.constrained$krel_FC_mean

newSummary.TRP.Kpre.constrained = newSummary.TRP.Kpre.constrained[order(newSummary.TRP.Kpre.constrained$IndexDisp, decreasing=T),]
newSummary.TRP.Kpre.constrained$kinitFCcategory = "Undefined"
newSummary.TRP.Kpre.constrained[newSummary.TRP.Kpre.constrained$kinit_FC_max < 1, ]$kinitFCcategory = "kinitFC_less_1"
newSummary.TRP.Kpre.constrained[newSummary.TRP.Kpre.constrained$kinit_FC_max >= 1 & newSummary.TRP.Kpre.constrained$kinit_FC_min <= 1, ]$kinitFCcategory = "kinitFC_span_1"

newSummary.TRP.Kpre.constrained.kinitfc.less1 = subset(newSummary.TRP.Kpre.constrained, kinit_FC_max < 1)
newSummary.TRP.Kpre.constrained.kinitfc.span1 = subset(newSummary.TRP.Kpre.constrained, kinit_FC_min <= 1 & kinit_FC_max >= 1)

colNames = c("gene", "kinit_FC_min", "kinit_FC_max", "krel_FC_min","krel_FC_max","IndexDisp")

newSummary.TRP.Kpre.constrained.kinitfc.span1[, colNames]
min(newSummary.TRP.Kpre.constrained.kinitfc.less1$kinit_FC_min)
max(newSummary.TRP.Kpre.constrained.kinitfc.less1$kinit_FC_max)
min(newSummary.TRP.Kpre.constrained.kinitfc.less1$krel_FC_min)
max(newSummary.TRP.Kpre.constrained.kinitfc.less1$krel_FC_max)

gene.list.Init = c("Nup214", "Fgfr2", "Clhc1", "Dmrta2")
colorManual.one = c("#FF7F0E", "#1F77B4", "#E7298A", "#1B9E77") # blue, orange, dark palette

subset.rateData.TRP.kpreConstrain = subset(rateData.Kpre.constrained, gene %in% gene.list.Init)
subset.rateData.TRP.kpreConstrain$kinitFCcategory = NA

for (gene.name in gene.list.Init) {
	subset.rateData.TRP.kpreConstrain[subset.rateData.TRP.kpreConstrain$gene == gene.name, ]$kinitFCcategory = subset(newSummary.TRP.Kpre.constrained, gene == gene.name)$kinitFCcategory 
}
subset.rateData.TRP.kpreConstrain$kinitFCcategory = factor(subset.rateData.TRP.kpreConstrain$kinitFCcategory, levels = c("kinitFC_less_1", "kinitFC_span_1"))

subset.rateData.TRP.kpreConstrain$gene = factor(subset.rateData.TRP.kpreConstrain$gene, levels = gene.list.Init)

# kinit FC changes
colorManual.one = c("#1F77B4", "#FF7F0E", "#E7298A", "#1B9E77") # orange, blue, dark palette
ggplot(subset.rateData.TRP.kpreConstrain) + 
	geom_violin(aes(x=gene, y=kinit_FC, fill=kinitFCcategory), scale="width") +
        geom_point(aes(x=gene,y=kinit_FC), position="jitter", alpha=0.3 , size=0.16) +
	scale_y_continuous(trans=log2_trans(), breaks=trans_breaks('log2', function(x) 2^x), labels=trans_format('log2', math_format(.x)), limits = c(NA,1)) +
        #ylab(expression(atop(log[2] ~ "Foldchange in" ~ K[init] ~ ", TRP treatment")))  + 
        ylab(expression(log[2] ~ "Foldchange in" ~ italic(K[init])))  +
        scale_fill_manual(values = colorManual.one,
                           labels= c("Decrease Initiation Rate", "Uncertain Change")) + 
	theme_bw() + fontTheme +
        theme(legend.position="none",axis.text.x = element_text(size=32, face=1.5, angle=30, vjust=0.9, hjust=1),
               axis.title.x = element_blank(), axis.title.y = element_text(size=32), legend.title = element_blank())

ggsave("log2_TRP_kinitFC_distribution_Kpre_constrained.pdf", height=8.4, width = 6)

# kinit FC index of dispersion

ggplot() +
        geom_violin(data=newSummary.TRP.Kpre.constrained, aes(x=kinitFCcategory, y = IndexDisp, fill = kinitFCcategory), scale="width") +
        geom_boxplot(data=newSummary.TRP.Kpre.constrained, aes(x=kinitFCcategory, y = IndexDisp),alpha = 0.1, width=0.36, linewidth=0.3, outlier.shape=NA) +

        geom_point(data=subset(newSummary.TRP.Kpre.constrained, gene %in% gene.list.Init), aes(x=kinitFCcategory, y = IndexDisp, fill = kinitFCcategory), size=2) +
       geom_text(data=subset(newSummary.TRP.Kpre.constrained, gene %in% gene.list.Init), aes(x=kinitFCcategory, y = IndexDisp, label=gene, fontface=2), nudge_x=-0.12, nudge_y=0.3, size=7) +
       #scale_y_continuous(trans=log2_trans(), breaks=trans_breaks('log2', function(x) 2^x), labels=trans_format('log2', math_format(.x))) +
       scale_y_continuous(trans=log10_trans(), breaks=trans_breaks('log10', function(x) 10^x), labels=trans_format('log10', math_format(.x)), limits = c(NA,1)) +
       scale_fill_d3() +
      scale_color_d3() +
                  theme_bw() + fontTheme +
        #ggtitle("Gene Categories - Triptolide 12.5 min") +
        theme(plot.title = element_text(size=23),
              axis.text.x = element_text(size=31),
              axis.text.y = element_text(size=26),
               axis.title.y = element_text(size=29),
              axis.title.x = element_blank(),
             # axis.text.x = element_text(size=28, angle=30, hjust=0.1, vjust=0.22),
                legend.position = "None") +
                  labs(fill = "Gene Categories") +
        #ylab("Index of Dispersion - TRP treatment")
        #ylab("log2 Index of Dispersion - TRP treatment")
        #ylab(expression(atop(log[10] ~ frac("Variance", "Mean") ~ "Foldchange in Initiation, TRP treatment", "(valid parameter sets, TRP treatment)"))) +
        ylab(expression(log[10] ~ frac("Variance (Foldchange" ~ italic(K[init]) ~")", "Mean (Foldchange" ~ italic(K[init]) ~")"))) +
        scale_x_discrete(labels = c("Decrease\nInitiation", "Uncertain\nChange"))

ggsave("log10_TRP_kinitFC_indexOfDispersion_Kpre_constrained.pdf", width=6, height=8.4) # 6.97 x 8.4
### histograms
colorManual.two = c("#1F77B4", "#FF7F0E", "#E7298A", "#1B9E77") # orange, blue, dark palette
df.TRP.kpreconstrain = as.data.frame(cbind(c("kinitFCless1", "kinitFCspan1"),
                   c(sum(newSummary.TRP.Kpre.constrained$kinitFCcategory == "kinitFC_less_1"),
                     sum(newSummary.TRP.Kpre.constrained$kinitFCcategory == "kinitFC_span_1"))))
                               #sum(newSummary.TRP$kinitFCcategory == "kinitFC_more_1"))))
colnames(df.TRP.kpreconstrain) = c("typeOfData", "countOfData")
df.TRP.kpreconstrain$typeOfData = factor(df.TRP.kpreconstrain$typeOfData, levels = c("kinitFCless1", "kinitFCspan1"))
df.TRP.kpreconstrain$countOfData = as.numeric(df.TRP.kpreconstrain$countOfData)
sum_all = sum(df.TRP.kpreconstrain$countOfData)
df.TRP.kpreconstrain$countOfData = df.TRP.kpreconstrain$countOfData / sum_all

ggplot(df.TRP.kpreconstrain, aes(x = "All Genes", y = countOfData, fill=typeOfData)) +
         #geom_bar(stat="identity", position=position_fill(reverse=T)) +
         geom_bar(stat="identity", position=position_stack(reverse=T)) +
         #scale_y_continuous(labels = scales::percent) +
         ylab("Fraction of Genes") +
         #labs(fill = "") +
         #xlab("") +
         #scale_fill_brewer(palette = "Dark2",
         scale_fill_manual(values = colorManual.two,
                           #labels= c("Decrease Initiation Rate", "Ambiguous Change", "Increase Initiation Rate")) +
                           labels= c("Decrease Initiation Rate", "Uncertain Change")) +
         # scale_color_manual(values = colorManual.one,
        #                  labels= c("Decrease Initiation Rate", "Uncertain Change")) +
                           #  labels= c("Decrease Initiation Rate", "Ambiguous Change", "Increase Initiation Rate")) +
         theme_bw() + fontTheme +
         theme(legend.justification = "center",
               axis.text.x = element_text(size=31),
              axis.text.y = element_text(size=27),
               axis.title.y = element_text(size=32),
              axis.title.x = element_blank(), legend.title = element_blank(), aspect.ratio = 9)

ggsave("CountFrequency_stacked_KinitFC_profile_TRP_treat_Kpre_constrained.pdf", height=8.4, width=8) # 6.53 x 8.4

# distribution of mean Kinit, Krel, kpre, kelong
# first - compare between krel_base/krel_expr and kpre
ggplot(newSummary.TRP.Kpre.constrained) +
	geom_violin(aes(x="release_base_mean",y= Release_base_mean), scale="width") +
        geom_violin(aes(x="release_expr_mean",y= Release_expr_mean), scale="width") +
        geom_violin(aes(x="kinit_base_mean",y= kinit_base_mean), scale="width") +
        geom_violin(aes(x="kinit_expr_mean",y= kinit_expr_mean), scale="width") +
        geom_violin(aes(x="kpre_mean",y= kpre_mean), scale="width") +
        geom_violin(aes(x="krel_base_mean",y= krel_base_mean), scale="width") +
        geom_violin(aes(x="krel_expr_mean",y= krel_expr_mean), scale="width") +
        geom_violin(aes(x="terminal_base_mean",y= Terminal_base_mean), scale="width") +
        geom_violin(aes(x="terminal_expr_mean",y= Terminal_expr_mean), scale="width") +

        geom_point(aes(x="release_base_mean",y= Release_base_mean), position="jitter", alpha=0.02) +
        geom_point(aes(x="release_expr_mean",y= Release_expr_mean), position="jitter", alpha=0.02) +
        geom_point(aes(x="kinit_base_mean",y= kinit_base_mean), position="jitter", alpha=0.02) +
        geom_point(aes(x="kinit_expr_mean",y= kinit_expr_mean), position="jitter", alpha=0.02) +
        geom_point(aes(x="kpre_mean",y= kpre_mean), position="jitter", alpha=0.02) +
        geom_point(aes(x="krel_base_mean",y= krel_base_mean), position="jitter", alpha=0.02) +
        geom_point(aes(x="krel_expr_mean",y= krel_expr_mean), position="jitter", alpha=0.02) +
        geom_point(aes(x="terminal_base_mean",y= Terminal_base_mean), position="jitter", alpha=0.02) +
        geom_point(aes(x="terminal_expr_mean",y= Terminal_expr_mean), position="jitter", alpha=0.02) +

        geom_boxplot(aes(x="release_base_mean",y= Release_base_mean), alpha=0.1, width=0.36, linewidth=0.3) +
        geom_boxplot(aes(x="release_expr_mean",y= Release_expr_mean), alpha=0.1, width=0.36, linewidth=0.3)+
        geom_boxplot(aes(x="kinit_base_mean",y= kinit_base_mean), alpha=0.1, width=0.36, linewidth=0.3) +
        geom_boxplot(aes(x="kinit_expr_mean",y= kinit_expr_mean), alpha=0.1, width=0.36, linewidth=0.3) +
        geom_boxplot(aes(x="kpre_mean",y= kpre_mean), alpha=0.1, width=0.36, linewidth=0.3) +
        geom_boxplot(aes(x="krel_base_mean",y= krel_base_mean), alpha=0.1, width=0.36, linewidth=0.3) +
        geom_boxplot(aes(x="krel_expr_mean",y= krel_expr_mean), alpha=0.1, width=0.36, linewidth=0.3) +
        geom_boxplot(aes(x="terminal_base_mean",y= Terminal_base_mean), alpha=0.1, width=0.36, linewidth=0.3) +
        geom_boxplot(aes(x="terminal_expr_mean",y= Terminal_expr_mean), alpha=0.1, width=0.36, linewidth=0.3)+
                  scale_y_continuous(trans=log10_trans(),
                breaks=trans_breaks('log10', function(x) 10^x),
                labels=trans_format('log10', math_format(.x))) +
         scale_x_discrete(labels = c(expression(atop(italic(K[init]), "(control)")), expression(atop(italic(K[init]), "(triptolide)")),
                                     expression(italic(K[pre])),expression(atop(italic(K[rel]), "(control)")),expression(atop(italic(K[rel]), "(triptolide)")),
                                     expression(atop("Effective" ~ italic(K[rel]), "(control)")),expression(atop("Effective" ~ italic(K[rel]), "(triptolide)")),
                                     expression(atop("Effective" ~ italic(K[pre]), "(control)")), expression(atop("Effective" ~ italic(K[pre]), "(triptolide)")))) + # kinit_base_mean, kinit_expr_mean, kpre_mean, krel_base_mean,krel_expr_mean, release_base_mean,release_expr_mean, terminal_base_mean, terminal_expr_mean
         ylab(expression(log[10] ~ "Mean Raw Values - All Genes")) +
         theme_bw() + fontTheme + theme(axis.title.x = element_blank(),
         axis.text.x = element_text(size=32, face=1, vjust=0.9))

ggsave("TRP_trends_Kpre_Krel_withconstrainedKpre.pdf", width = 26, height=8)
## plot P values trend
ggplot(newSummary.TRP.Kpre.constrained) + 
	geom_violin(aes(x="P_base_mean",y= P_base_mean)) +
	geom_violin(aes(x="P_expr_mean",y= P_expr_mean)) +
	geom_violin(aes(x="B_base_mean",y= B_base_mean)) +
	 geom_violin(aes(x="B_expr_mean",y= B_expr_mean)) +

	geom_point(aes(x="P_base_mean",y= P_base_mean), position="jitter", alpha=0.05) +
	geom_point(aes(x="P_expr_mean",y= P_expr_mean), position="jitter", alpha=0.05) + 
	 geom_point(aes(x="B_base_mean",y= B_base_mean), position="jitter", alpha=0.05) +
	 geom_point(aes(x="B_expr_mean",y= B_expr_mean), position="jitter", alpha=0.05) +

	geom_boxplot(aes(x="P_base_mean",y= P_base_mean), alpha=0.1, width=0.36, linewidth=0.3) +
	geom_boxplot(aes(x="P_expr_mean",y= P_expr_mean), alpha=0.1, width=0.36, linewidth=0.3) +
	geom_boxplot(aes(x="B_base_mean",y= B_base_mean), alpha=0.1, width=0.36, linewidth=0.3) +
	geom_boxplot(aes(x="B_expr_mean",y= B_expr_mean), alpha=0.1, width=0.36, linewidth=0.3) +	
		 scale_y_continuous(trans=log10_trans(), 
		breaks=trans_breaks('log10', function(x) 10^x), 
		labels=trans_format('log10', math_format(.x)), limits = c(NA,1)) +
         ylab(expression(log[10] ~ "raw values - all genes")) + 
	 theme_bw() + fontTheme + theme(axis.title.x = element_blank(), 
					axis.text.x = element_text(size=26, angle=30, vjust=1, hjust=0.9))

ggsave("TRP_trends_P_B_Vals_withconstrainedKpre.pdf")

### plot Krel trends
newSummary.TRP.Kpre.constrained$krelFCcategory = "Undefined"
newSummary.TRP.Kpre.constrained[newSummary.TRP.Kpre.constrained$krel_FC_max < 1, ]$krelFCcategory = "krelFC_less_1"
newSummary.TRP.Kpre.constrained[newSummary.TRP.Kpre.constrained$krel_FC_min >= 1, ]$krelFCcategory = "krelFC_more_1"

subset.rateData.TRP.kpreConstrain = subset(rateData.Kpre.constrained, gene %in% gene.list.Init)
subset.rateData.TRP.kpreConstrain$krelFCcategory = NA

for (gene.name in gene.list.Init) {
	subset.rateData.TRP.kpreConstrain[subset.rateData.TRP.kpreConstrain$gene == gene.name, ]$krelFCcategory = subset(newSummary.TRP.Kpre.constrained, gene == gene.name)$krelFCcategory 
}

subset.rateData.TRP.kpreConstrain$krelFCcategory = factor(subset.rateData.TRP.kpreConstrain$krelFCcategory, levels = c("krelFC_less_1", "krelFC_more_1"))

subset.rateData.TRP.kpreConstrain$gene = factor(subset.rateData.TRP.kpreConstrain$gene, levels = gene.list.Init)

# krel FC changes

ggplot(subset.rateData.TRP.kpreConstrain) + 
	geom_violin(aes(x=gene, y=krel_FC, fill=krelFCcategory), scale="width") +
        geom_point(aes(x=gene,y=krel_FC, color=krelFCcategory), position="jitter", alpha=0.3 , size=0.16) +
	scale_y_continuous(trans=log2_trans(), breaks=trans_breaks('log2', function(x) 2^x), labels=trans_format('log2', math_format(.x))) +
        ylab(expression(log[2] ~ "Foldchange in" ~ italic(K[rel])))  +
        scale_fill_manual(values = c("#E7298A", "#1B9E77"),
                          labels = c("decrease Krel", "increase Krel")) +
        scale_color_manual(values = c("#E7298A", "#1B9E77")) +
	theme_bw() + fontTheme +
        theme(legend.position="none",axis.text.x = element_text(size=32, face=1.5, angle=30, vjust=0.9, hjust=1),
               axis.title.x = element_blank(), axis.title.y = element_text(size=32), legend.title = element_blank())

ggsave("log2_TRP_krelFC_distribution_Kpre_constrained.pdf", height=8.4, width = 6)

## krel FC index of dispersion

ggplot() +
        geom_violin(data=newSummary.TRP.Kpre.constrained, aes(x=krelFCcategory, y = IndexDisp_Krel, fill = krelFCcategory), scale="width") +
        geom_boxplot(data=newSummary.TRP.Kpre.constrained, aes(x=krelFCcategory, y = IndexDisp_Krel),alpha = 0.1, width=0.36, linewidth=0.3, outlier.shape=NA) +

        geom_point(data=subset(newSummary.TRP.Kpre.constrained, gene %in% gene.list.Init), aes(x=krelFCcategory, y = IndexDisp_Krel, fill = krelFCcategory), size=2) +
       geom_text(data=subset(newSummary.TRP.Kpre.constrained, gene %in% gene.list.Init), aes(x=krelFCcategory, y = IndexDisp_Krel, label=gene, fontface=2), nudge_x=-0.12, nudge_y=0.3, size=7) +
       #scale_y_continuous(trans=log2_trans(), breaks=trans_breaks('log2', function(x) 2^x), labels=trans_format('log2', math_format(.x))) +
       scale_y_continuous(trans=log10_trans(), breaks=trans_breaks('log10', function(x) 10^x), labels=trans_format('log10', math_format(.x)), limits = c(NA,1)) +
       scale_fill_manual(values = c("#E7298A", "#1B9E77"),
                          labels = c("decrease Krel", "increase Krel")) +
       scale_color_manual(values = c("#E7298A", "#1B9E77")) +
        theme_bw() + fontTheme +
        #ggtitle("Gene Categories - Triptolide 12.5 min") +
        theme(plot.title = element_text(size=23),
              axis.text.x = element_text(size=31),
              axis.text.y = element_text(size=26),
               axis.title.y = element_text(size=29),
              axis.title.x = element_blank(),
             # axis.text.x = element_text(size=28, angle=30, hjust=0.1, vjust=0.22),
                legend.position = "None") +
                  labs(fill = "Gene Categories") +
        #ylab("Index of Dispersion - TRP treatment")
        #ylab("log2 Index of Dispersion - TRP treatment")
        #ylab(expression(atop(log[10] ~ frac("Variance", "Mean") ~ "Foldchange in Initiation, TRP treatment", "(valid parameter sets, TRP treatment)"))) +
        ylab(expression(log[10] ~ frac("Variance (Foldchange" ~ italic(K[rel]) ~")", "Mean (Foldchange" ~ italic(K[rel]) ~")"))) +
        scale_x_discrete(labels = c(expression(atop("Decrease", italic(K[rel]))), expression(atop("Increase", italic(K[rel])))))

ggsave("log10_TRP_krelFC_indexOfDispersion_Kpre_constrained.pdf", width=6, height=8.4) # 6.97 x 8.4

# krel for kinitFC_less_1 and kinitFC_span_1 
TRP.kpreConstrain.kinitFCless1.krelFCless1 = subset(newSummary.TRP.Kpre.constrained, kinit_FC_max < 1 & krel_FC_max < 1)
TRP.kpreConstrain.kinitFCless1.krelFCmore1 = subset(newSummary.TRP.Kpre.constrained, kinit_FC_max < 1 & krel_FC_min >= 1)

df.TRP.kpreConstrain.krel.kinit = as.data.frame(cbind(c("krelFCless1", "krelFCmore1"),
                             c(nrow(TRP.kpreConstrain.kinitFCless1.krelFCless1) ,
			       nrow(TRP.kpreConstrain.kinitFCless1.krelFCmore1))))
colnames(df.TRP.kpreConstrain.krel.kinit) = c("typeOfData", "countOfData")
df.TRP.kpreConstrain.krel.kinit$typeOfData = factor(df.TRP.kpreConstrain.krel.kinit$typeOfData, levels = c("krelFCless1", "krelFCmore1"))
df.TRP.kpreConstrain.krel.kinit$countOfData = as.numeric(df.TRP.kpreConstrain.krel.kinit$countOfData)
sum_all = sum(df.TRP.kpreConstrain.krel.kinit$countOfData)
df.TRP.kpreConstrain.krel.kinit$countOfData = df.TRP.kpreConstrain.krel.kinit$countOfData / sum_all
df.TRP.kpreConstrain.krel.kinit$class = "Decreasing Kinit"
TRP.kpreConstrain.kinitFCspan1.krelFCless1 = subset(newSummary.TRP.Kpre.constrained, kinit_FC_max >= 1 & kinit_FC_min < 1 & krel_FC_max < 1)
TRP.kpreConstrain.kinitFCspan1.krelFCmore1 = subset(newSummary.TRP.Kpre.constrained, kinit_FC_max >= 1 & kinit_FC_min < 1 & krel_FC_min >= 1)

tmp.df = as.data.frame(cbind(c("krelFCless1", "krelFCmore1"),
                             c(nrow(TRP.kpreConstrain.kinitFCspan1.krelFCless1) , 
			       nrow(TRP.kpreConstrain.kinitFCspan1.krelFCmore1))))

colnames(tmp.df) = c("typeOfData", "countOfData")
tmp.df$countOfData = as.numeric(tmp.df$countOfData)
sum_all = sum(tmp.df$countOfData)
tmp.df$countOfData = tmp.df$countOfData / sum_all
tmp.df$class = "Uncertain\nChange"
df.TRP.kpreConstrain.krel.kinit = rbind(df.TRP.kpreConstrain.krel.kinit, tmp.df)
colors.New = c("#E7298A", "#1B9E77")
ggplot(df.TRP.kpreConstrain.krel.kinit, aes(x=class, y = countOfData, fill=typeOfData)) +
         #geom_bar(stat="identity", position=position_fill(reverse=T)) +
         geom_bar(stat="identity", position=position_stack(reverse=F), width=0.6) +
         #scale_y_continuous(labels = scales::percent) +
         ylab("Fraction of Genes") +
         #ylab(expression("Fraction of Genes (with" ~ italic(K[init]) ~ " decrease)")) +
         #xlab(expression("Genes with decreasing" ~ italic(K[init]) )) +
         #labsForLegend +
         #scale_fill_brewer(palette = "Dark2",
        scale_fill_manual(values = colors.New ,
                          labels= c(expression("Decrease" ~ italic(K[rel])), expression("Increase" ~ italic(K[rel])))) +
        #scale_fill_manual(values = colorManual.two,
        #              labels= c("Decrease Pause Release", "Increase Pause Release")) +
         #scale_fill_frontiers(labels= c("Decrease Pause Release", "Increase Pause Release")) +
        #labs(fill = "Gene Categories") +
        scale_x_discrete(labels = c(expression(atop("Decrease","in" ~ italic(K[init]))),
                                expression(atop("Uncertain",italic(K[init]) ~ "Change")))) +
        theme_bw() + fontTheme +
         theme(legend.justification = "center",
               axis.ticks.length.x.bottom = unit(4, "mm"),
	       axis.title.x= element_blank(),
               axis.text.x = element_text(size=30, angle=90, face=1.3, vjust=0.65, hjust=0.5),
              axis.text.y = element_text(size=26),
               axis.title.y = element_text(size=32), legend.title = element_blank())
ggsave("CountFrequency_stacked_KrelFC_TRP_treat_Kpre_constrained.pdf", height=10.4)

# plot absolute raets for gene Nup214 (selected). Other genes are Fgfr2 Clhc1 Dmrta2.
subset.rateData.TRP.kpreConstrain$P_base = subset.rateData.TRP.kpreConstrain$kinit_base / (subset.rateData.TRP.kpreConstrain$kpre + subset.rateData.TRP.kpreConstrain$krel_base)
subset.rateData.TRP.kpreConstrain$P_expr = subset.rateData.TRP.kpreConstrain$kinit_expr / (subset.rateData.TRP.kpreConstrain$kpre + subset.rateData.TRP.kpreConstrain$krel_expr)

ggplot(subset(subset.rateData.TRP.kpreConstrain, gene %in% c("Nup214"))) + 
#ggplot(subset(subset.rateData.TRP.kpreConstrain, gene %in% c("Fgfr2"))) + 
	geom_violin(aes(x="Premature\nTermination",y= kpre)) +
	#geom_violin(aes(x="Elongation",y= kelong)) + 
	geom_violin(aes(x="Initiation\n(control)",y= kinit_base)) +
	geom_violin(aes(x="Initiation\n(triptolide)",y= kinit_expr)) +
	geom_violin(aes(x="Pause Release\n(control)",y= krel_base)) +
	geom_violin(aes(x="Pause Release\n(triptolide)",y= krel_expr)) +
	geom_violin(aes(x="release_base",y= Release_base), scale="width") +
        geom_violin(aes(x="release_expr",y= Release_expr), scale="width") +
	geom_violin(aes(x="terminal_base",y= Terminal_base), scale="width") +
        geom_violin(aes(x="terminal_expr",y= Terminal_expr), scale="width") +


        geom_point(aes(x="release_base",y= Release_base), position="jitter", alpha=0.1) +
        geom_point(aes(x="release_expr",y= Release_expr), position="jitter", alpha=0.1) +
	geom_point(aes(x="terminal_base",y= Terminal_base), position="jitter", alpha=0.1) +
        geom_point(aes(x="terminal_expr",y= Terminal_expr), position="jitter", alpha=0.1) +
	geom_point(aes(x="Premature\nTermination",y= kpre), position="jitter", alpha=0.1) +
	geom_point(aes(x="Initiation\n(control)",y= kinit_base), position="jitter", alpha=0.1) + 
	geom_point(aes(x="Initiation\n(triptolide)",y= kinit_expr), position="jitter", alpha=0.1) +
	geom_point(aes(x="Pause Release\n(control)",y= krel_base), position="jitter", alpha=0.1) +
	geom_point(aes(x="Pause Release\n(triptolide)",y= krel_expr), position="jitter", alpha=0.1) +

	geom_boxplot(aes(x="terminal_base",y= Terminal_base), alpha=0.1, width=0.36, linewidth=0.3) +
        geom_boxplot(aes(x="terminal_expr",y= Terminal_expr), alpha=0.1, width=0.36, linewidth=0.3)+
	geom_boxplot(aes(x="release_base",y= Release_base), alpha=0.1, width=0.36, linewidth=0.3) +
        geom_boxplot(aes(x="release_expr",y= Release_expr), alpha=0.1, width=0.36, linewidth=0.3)+
	geom_boxplot(aes(x="Premature\nTermination",y= kpre), alpha=0.1, width=0.36, linewidth=0.3) +
	geom_boxplot(aes(x="Initiation\n(control)",y= kinit_base), alpha=0.1, width=0.36, linewidth=0.3) + 
	geom_boxplot(aes(x="Initiation\n(triptolide)",y= kinit_expr), alpha=0.1, width=0.36, linewidth=0.3) + 
	geom_boxplot(aes(x="Pause Release\n(control)",y= krel_base), alpha=0.1, width=0.36, linewidth=0.3) +
	geom_boxplot(aes(x="Pause Release\n(triptolide)",y= krel_expr), alpha=0.1, width=0.36, linewidth=0.3) +
		 scale_y_continuous(trans=log10_trans(), 
		breaks=trans_breaks('log10', function(x) 10^x), 
		labels=trans_format('log10', math_format(.x)), limits=c(NA,1)) +
         ylab(expression(log[10] ~ "Raw values for gene" ~ italic("Nup214"))) + 
         #ylab(expression(log[10] ~ "Raw values for gene" ~ italic("Fgrf2"))) + 
	 scale_x_discrete(labels = c(expression(atop(italic(K[init]), "(control)")),
				     expression(atop(italic(K[init]), "(triptolide)")),
				     expression(atop(italic(K[rel]), "(control)")),
				     expression(atop(italic(K[rel]), "(triptolide)")),
				     expression(italic(K[pre])),
                                     expression(atop("Effective" ~ italic(K[rel]), "(control)")),expression(atop("Effective" ~ italic(K[rel]), "(triptolide)")),
                                     expression(atop("Effective" ~ italic(K[pre]), "(control)")), expression(atop("Effective" ~ italic(K[pre]), "(triptolide)")))) + 
	 theme_bw() + fontTheme + 
	 theme(axis.title.x = element_blank(),
	 axis.text.x = element_text(size=30, face=1.3),
         axis.text.y = element_text(size=26))
ggsave("TRP_treatment_values_gene_Nup214_kpre_constrained.pdf", height=8.4, width=26)
#ggsave("TRP_treatment_values_gene_Fgrf2_kpre_constrained.pdf", height=8.4)

# plot pause sums
ggplot(subset(subset.rateData.TRP.kpreConstrain, gene %in% c("Nup214"))) + 
	geom_violin(aes(x="Initiation\n(control)",y= kinit_base)) +
	geom_violin(aes(x="Elongation",y= kelong)) +
	geom_violin(aes(x="Pause Sum\n(control)",y= P_base)) +
	geom_violin(aes(x="Pause Sum\n(triptolide)",y= P_expr)) +
	
	geom_point(aes(x="Initiation\n(control)",y= kinit_base), position="jitter", alpha=0.1) + 
	geom_point(aes(x="Elongation",y= kelong), position="jitter", alpha=0.1) +
	geom_point(aes(x="Pause Sum\n(control)",y= P_base), position="jitter", alpha=0.1) +
	geom_point(aes(x="Pause Sum\n(triptolide)",y= P_expr), position="jitter", alpha=0.1) +

	geom_boxplot(aes(x="Initiation\n(control)",y= kinit_base), alpha=0.1, width=0.36, linewidth=0.3) + 
	geom_boxplot(aes(x="Elongation",y= kelong), alpha=0.1, width=0.36, linewidth=0.3) + 
	geom_boxplot(aes(x="Pause Sum\n(control)",y= P_base), alpha=0.1, width=0.36, linewidth=0.3) +
	geom_boxplot(aes(x="Pause Sum\n(triptolide)",y= P_expr), alpha=0.1, width=0.36, linewidth=0.3) +
		 scale_y_continuous(trans=log2_trans(), 
		breaks=trans_breaks('log2', function(x) 2^x), 
		labels=trans_format('log2', math_format(.x))) +
         ylab(expression(log[2] ~ "Values for gene" ~ italic("Nup214"))) + 
	 scale_x_discrete(labels = c(expression(italic(K[elong])),
				     expression(atop(italic(K[init]), "(control)")),
				     expression(atop("Pause Sum", "(control)")),
				     expression(atop("Pause Sum", "(triptolide)"))
				     )) + 
	 theme_bw() + fontTheme + 
	 theme(axis.title.x = element_blank(),
	 axis.text.x = element_text(size=30, face=1.3),
         axis.text.y = element_text(size=26))


ggsave("TRP_treatment_P_values_gene_Nup214_kpre_constrained.pdf", height=8.4)

