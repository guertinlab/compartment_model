library(ggplot2)
library(dplyr)
library(ggsci)
library(RColorBrewer)
library(scales)
# https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html
source("plotting_composites_lattice.R")

percentTotal_trans <- function(x){
    TOTAL = 7527 # change it to the total you are using
    paste(round((x / TOTAL)*100, digits=2), "%")
}
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

density.TRP = readRDS("Jonkers_2014_Trp_FP_paper/density.TRP_12.5min_all.v.control_all.rds")
rateData.TRP = readRDS("InternalControl_TRP12min_repressed/repeatedParamScanData_TRP_12min_repressed_allParamsConstrained.rds") # kinit, kpre, krel between 1e-8 and 1.
## plot pause sum and gene body avg
rateData.TRP$P_base = rateData.TRP$kinit_base / (rateData.TRP$kpre + rateData.TRP$krel_base)
rateData.TRP$P_expr = rateData.TRP$kinit_expr / (rateData.TRP$kpre + rateData.TRP$krel_expr)
rateData.TRP$B_base = (rateData.TRP$krel_base * rateData.TRP$P_base) / rateData.TRP$kelong
rateData.TRP$B_expr = (rateData.TRP$krel_expr * rateData.TRP$P_expr) / rateData.TRP$kelong
rateData.TRP$Release_base = rateData.TRP$P_base * rateData.TRP$krel_base
rateData.TRP$Release_expr = rateData.TRP$P_expr * rateData.TRP$krel_expr
rateData.TRP$Terminal_base = rateData.TRP$kpre * rateData.TRP$P_base
rateData.TRP$Terminal_expr = rateData.TRP$kpre * rateData.TRP$P_expr

newSummary.TRP = group_by(rateData.TRP, gene) %>%  summarise_at(c("Release_base", "Release_expr", "P_base", "P_expr", "B_base", "B_expr", "kinit_base", "kinit_expr", "krel_base", "krel_expr", "kpre", "kelong", "kinit_FC", "krel_FC"), c("min", "max", "var", "mean"))

# update category of kinit and krel foldchange
newSummary.TRP$kinitFCcategory = "Undefined"
newSummary.TRP[newSummary.TRP$kinit_FC_max < 1, ]$kinitFCcategory = "kinitFC_less_1"
newSummary.TRP[newSummary.TRP$kinit_FC_min <= 1 & newSummary.TRP$kinit_FC_max > 1, ]$kinitFCcategory = "kinitFC_span_1"
#newSummary.TRP[newSummary.TRP$kinit_FC_min > 1, ]$kinitFCcategory = "kinitFC_more_1"
newSummary.TRP$kinitFCcategory = factor(newSummary.TRP$kinitFCcategory, levels = c("kinitFC_less_1", "kinitFC_span_1"))
# krel FC category
newSummary.TRP$krelFCcategory = NA
newSummary.TRP[newSummary.TRP$krel_FC_min >= 1, ]$krelFCcategory = "krelFCmore1"
newSummary.TRP[newSummary.TRP$krel_FC_max <= 1, ]$krelFCcategory = "krelFCless1"
subset(newSummary.TRP, is.na(krelFCcategory))[,c("gene","krel_FC_min", "krel_FC_max")] # 5 genes with krel_FC_max hovering around 1
newSummary.TRP[is.na(newSummary.TRP$krelFCcategory), ]$krelFCcategory = "krelFCless1"

newSummary.TRP$krelFCcategory = factor(newSummary.TRP$krelFCcategory,levels = c("krelFCless1", "krelFCmore1") )


# bar plot for initiation foldchange categories
df.TRP = as.data.frame(cbind(c("kinitFCless1", "kinitFCspan1"),
			     c(sum(newSummary.TRP$kinitFCcategory == "kinitFC_less_1"), 
			       sum(newSummary.TRP$kinitFCcategory == "kinitFC_span_1"))))
			       #sum(newSummary.TRP$kinitFCcategory == "kinitFC_more_1"))))
colnames(df.TRP) = c("typeOfData", "countOfData")
df.TRP$typeOfData = factor(df.TRP$typeOfData, levels = c("kinitFCless1", "kinitFCspan1"))
df.TRP$countOfData = as.numeric(df.TRP$countOfData)
sum_all = sum(df.TRP$countOfData)
df.TRP$countOfData = df.TRP$countOfData / sum_all 

colorManual.one <- c("#1F77B4", "#FF7F0E", "#2CA02C")

ggplot(df.TRP, aes(x = "All Genes", y = countOfData, fill=typeOfData)) +
         #geom_bar(stat="identity", position=position_fill(reverse=T)) +
         geom_bar(stat="identity", position=position_stack(reverse=T)) +
	 #scale_y_continuous(labels = scales::percent) + 
	 ylab("Fraction of Genes") +
	 #labs(fill = "") + 
	 #xlab("") + 
	 #scale_fill_brewer(palette = "Dark2",
	 scale_fill_manual(values = colorManual.one,
			   #labels= c("Decrease Initiation Rate", "Ambiguous Change", "Increase Initiation Rate")) +
			   labels= c("Decrease Initiation Rate", "Uncertain Change")) +
	 # scale_color_manual(values = colorManual.one,
	#		   labels= c("Decrease Initiation Rate", "Uncertain Change")) +
			   #  labels= c("Decrease Initiation Rate", "Ambiguous Change", "Increase Initiation Rate")) +
	 theme_bw() + fontTheme +
	 theme(legend.justification = "center",
	       axis.text.x = element_text(size=31, face = 1.3),
              axis.text.y = element_text(size=27),
               axis.title.y = element_text(size=32),
	      axis.title.x = element_blank(), legend.title = element_blank(), aspect.ratio = 9)
#ggsave("CountFrequency_KinitFC_profile_TRP_treat.pdf")
ggsave("CountFrequency_stacked_KinitFC_profile_TRP_treat.pdf", height=8.4) # 6.67 x 8

## selection of genes which have wide and narrow distribution of Kinit FC
newSummary.TRP$IndexDisp = newSummary.TRP$kinit_FC_var / newSummary.TRP$kinit_FC_mean
newSummary.TRP = newSummary.TRP[order(newSummary.TRP$IndexDisp, decreasing=T),]
TRP.kinitFC.less1 = subset(newSummary.TRP, kinit_FC_max < 1)
TRP.kinitFC.span1 = subset(newSummary.TRP, kinit_FC_max >= 1 & kinit_FC_min <= 1)
colNames = c("gene", "kinit_FC_max", "kinit_FC_min", "IndexDisp", "strand")
TRP.kinitFC.less1 = TRP.kinitFC.less1[order(TRP.kinitFC.less1$IndexDisp, decreasing = T),]
TRP.kinitFC.span1 = TRP.kinitFC.span1[order(TRP.kinitFC.span1$IndexDisp, decreasing = T),]
plus.TRP.kinitFC.less1 = subset(TRP.kinitFC.less1, strand == "+")
head(plus.TRP.kinitFC.less1[,colNames], n=12) # Nup214 (selected)
tail(TRP.kinitFC.less1[,colNames], n=12) #gene , Fgfr2 (selected), Top2a, Cox11, Cdkl3
TRP.kinitFC.span1 = TRP.kinitFC.span1[order(TRP.kinitFC.span1$IndexDisp, decreasing=T),]
head(TRP.kinitFC.span1[, colNames]) # Clhc1, Klhl35, Car2
tail(TRP.kinitFC.span1[,colNames], n=30) # Dmrta2, Gal
gene.list.Init = c("Nup214", "Fgfr2", "Clhc1", "Dmrta2")

# plot log2 kinit foldchange
colorManual.two = c("#1F77B4", "#FF7F0E", "#E7298A", "#1B9E77") # blue, orange, dark palette
subset.rateData.TRP = subset(rateData.TRP, gene %in% gene.list.Init)
#subset.rateData.TRP$kinitFCcategory = "undefined"
subset.rateData.TRP$kinitFCcategory = NA 
#subset.rateData.TRP$kinitFCcategory = factor(subset.rateData.TRP$kinitFCcategory, levels= c("kinitFC_less_1", "kinitFC_span_1","kinitFC_more_1"))
for (gene.name in gene.list.Init) {
 subset.rateData.TRP[subset.rateData.TRP$gene == gene.name,]$kinitFCcategory = subset(newSummary.TRP,gene == gene.name)$kinitFCcategory
}
subset.rateData.TRP$kinitFCcategory = factor(subset.rateData.TRP$kinitFCcategory, levels= c("kinitFC_less_1", "kinitFC_span_1"))
subset.rateData.TRP$gene = factor(subset.rateData.TRP$gene, levels = gene.list.Init)

ggplot(subset.rateData.TRP) + geom_violin(aes(x=gene, y=kinit_FC, fill=kinitFCcategory), scale="width") + 
	geom_point(aes(x=gene,y=kinit_FC), position="jitter", alpha=0.3 , size=0.16) + 
	scale_y_continuous(trans=log2_trans(), breaks=trans_breaks('log2', function(x) 2^x), labels=trans_format('log2', math_format(.x))) + 
	#ylab(expression(atop(log[2] ~ "Foldchange in" ~ K[init] ~ ", TRP treatment")))  + 
ylab(expression(log[2] ~ "Foldchange in" ~ italic(K[init])))  + 
	scale_fill_manual(values = colorManual.two,
			   labels= c("Decrease Initiation Rate", "Uncertain Change")) +
	theme_bw() + fontTheme + 
	theme(legend.position="none",axis.text.x = element_text(size=32, face=1.5, angle=30, vjust=0.9, hjust=1),
	       axis.title.x = element_blank(), axis.title.y = element_text(size=32), legend.title = element_blank())

ggsave("log2_TRP_kinitFC_distribution.pdf", height=8.4, width = 6)

# plot log10 kinit foldchnage index of dispersion
ggplot() + 
	geom_violin(data=newSummary.TRP, aes(x=kinitFCcategory, y = IndexDisp, fill = kinitFCcategory), scale="width") + 
	geom_boxplot(data=newSummary.TRP, aes(x=kinitFCcategory, y = IndexDisp),alpha = 0.1, width=0.36, linewidth=0.3, outlier.shape=NA) + 
	
	geom_point(data=subset(newSummary.TRP, gene %in% gene.list.Init), aes(x=kinitFCcategory, y = IndexDisp, fill = kinitFCcategory), size=2) + 
       geom_text(data=subset(newSummary.TRP, gene %in% gene.list.Init), aes(x=kinitFCcategory, y = IndexDisp, label=gene, fontface=2), nudge_x=-0.12, nudge_y=0.3, size=7) + 
       #scale_y_continuous(trans=log2_trans(), breaks=trans_breaks('log2', function(x) 2^x), labels=trans_format('log2', math_format(.x))) +     
       scale_y_continuous(trans=log10_trans(), breaks=trans_breaks('log10', function(x) 10^x), labels=trans_format('log10', math_format(.x))) +
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
	#scale_x_discrete(labels = c("Decrease\nInitiation", "Ambiguous\nChange", "Increase\nInitiation"))  
	#theme(aspect.ratio = 7 / 4)
#ggsave("TRP_kinitFC_indexOfDispersion.pdf")
#ggsave("log2_TRP_kinitFC_indexOfDispersion.pdf")
ggsave("log10_TRP_kinitFC_indexOfDispersion.pdf", width=6, height=8.4) # 6.97 x 8.4

## plot for krel foldchange
df.TRP.krel = as.data.frame(cbind(c("krelFCless1", "krelFCmore1"),
				  c(sum(newSummary.TRP$krelFCcategory == "krelFCless1"), 
			       sum(newSummary.TRP$krelFCcategory == "krelFCmore1"))))	
colnames(df.TRP.krel) = c("typeOfData", "countOfData")
df.TRP.krel$typeOfData = factor(df.TRP.krel$typeOfData, levels = c("krelFCless1", "krelFCmore1"))
df.TRP.krel$countOfData = as.numeric(df.TRP.krel$countOfData)
sum_all = sum(df.TRP.krel$countOfData)
df.TRP.krel$countOfData = df.TRP.krel$countOfData / sum_all 
colorManual.two = c("#C7A76C", "#86B875")
ggplot(df.TRP.krel, aes(x = "genes", y = countOfData, fill=typeOfData)) +
         #geom_bar(stat="identity", position=position_fill(reverse=T)) +
         geom_bar(stat="identity", position=position_stack(reverse=F)) +
	 #scale_y_continuous(labels = scales::percent) + 
	 ylab("Fraction of Genes") + 
	 #xlab("") + 
	 #labsForLegend + 
	 #scale_fill_brewer(palette = "Dark2",
	scale_fill_manual(values = colorManual.two, 
		       labels= c("Decrease Pause Release", "Increase Pause Release")) +
	 #scale_fill_frontiers(labels= c("Decrease Pause Release", "Increase Pause Release")) + 
	#labs(fill = "Gene Categories") + 
	theme_bw() + fontTheme + 
	 theme(legend.justification = "center", axis.text.x = element_text(size=31),
              axis.text.y = element_text(size=26),
               axis.title.y = element_text(size=29),
              axis.title.x = element_blank(), legend.title = element_blank(), aspect.ratio = 9)

#ggsave("CountFrequency_KinitFC_profile_TRP_treat.pdf")
ggsave("CountFrequency_stacked_KrelFC_TRP_treat.pdf", height = 8.4) # 

# compute how many genes have krel fC > 1 , krelFC < 1 for kinitFC < 1
TRP.kinitFCless1.krelFCless1 = subset(newSummary.TRP, kinitFCcategory == "kinitFC_less_1" & krelFCcategory ==  "krelFCless1")
TRP.kinitFCless1.krelFCmore1 = subset(newSummary.TRP, kinitFCcategory == "kinitFC_less_1" & krelFCcategory ==  "krelFCmore1")
df.TRP.krel.kinit = as.data.frame(cbind(c("krelFCless1", "krelFCmore1"),
                             c(nrow(TRP.kinitFCless1.krelFCless1) , nrow(TRP.kinitFCless1.krelFCmore1))))
colnames(df.TRP.krel.kinit) = c("typeOfData", "countOfData")
df.TRP.krel.kinit$typeOfData = factor(df.TRP.krel.kinit$typeOfData, levels = c("krelFCless1", "krelFCmore1"))
df.TRP.krel.kinit$countOfData = as.numeric(df.TRP.krel.kinit$countOfData)
sum_all = sum(df.TRP.krel.kinit$countOfData)
df.TRP.krel.kinit$countOfData = df.TRP.krel.kinit$countOfData / sum_all
df.TRP.krel.kinit$class = "Decreasing Kinit"
TRP.kinitFCspan1.krelFCless1 = subset(newSummary.TRP, kinitFCcategory == "kinitFC_span_1" & krelFCcategory ==  "krelFCless1")
TRP.kinitFCspan1.krelFCmore1 = subset(newSummary.TRP, kinitFCcategory == "kinitFC_span_1" & krelFCcategory ==  "krelFCmore1")
tmp.df = as.data.frame(cbind(c("krelFCless1", "krelFCmore1"),
                             c(nrow(TRP.kinitFCspan1.krelFCless1) , nrow(TRP.kinitFCspan1.krelFCmore1))))
colnames(tmp.df) = c("typeOfData", "countOfData")
tmp.df$countOfData = as.numeric(tmp.df$countOfData)
sum_all = sum(tmp.df$countOfData)
tmp.df$countOfData = tmp.df$countOfData / sum_all
tmp.df$class = "Uncertain\nChange"
df.TRP.krel.kinit = rbind(df.TRP.krel.kinit, tmp.df)
# for all genes
#df.TRP$class = c("All Genes", "All Genes")
#df.TRP.krel.kinit = rbind(df.TRP.krel.kinit, df.TRP)
#df.TRP.krel.kinit$class = factor(df.TRP.krel.kinit$class, levels = c("All Genes", "Decreasing Kinit", "Kinit span 1"))
#df.TRP.krel.kinit$typeOfData = factor(df.TRP.krel.kinit$typeOfData, levels = c("kinitFCspan1", "kinitFCless1", "krelFCless1", "krelFCmore1"))
colors.New = c("#E7298A", "#1B9E77")
ggplot(df.TRP.krel.kinit, aes(x=class, y = countOfData, fill=typeOfData)) +
         #geom_bar(stat="identity", position=position_fill(reverse=T)) +
         geom_bar(stat="identity", position=position_stack(reverse=F), width=0.6) +
	 #scale_y_continuous(labels = scales::percent) + 
	 ylab("Fraction of Genes") + 
	 #ylab(expression("Fraction of Genes (with" ~ italic(K[init]) ~ " decrease)")) + 
	 #xlab(expression("Genes with decreasing" ~ italic(K[init]) )) + 
	 #labsForLegend + 
	 #scale_fill_brewer(palette = "Dark2",
	 scale_fill_manual(values = colors.New,	       
			  labels= c(expression("Decrease" ~ italic(K[rel])), expression("Increase" ~ italic(K[rel])))) +
	#scale_fill_manual(values = colorManual.two, 
	#	       labels= c("Decrease Pause Release", "Increase Pause Release")) +
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
       #	aspect.ratio = 4)
ggsave("CountFrequency_stacked_KinitFC_KrelFC_TRP_treat.pdf", height = 10.4) 

# plot krel FC for the genes selected
gene.list = gene.list.Init
subset.rateData.TRP = subset(rateData.TRP, gene %in% gene.list)
subset.rateData.TRP$gene = factor(subset.rateData.TRP$gene, levels = gene.list)
subset.rateData.TRP$krelFCcategory = NA
for (gene.name in gene.list) {
 subset.rateData.TRP[subset.rateData.TRP$gene == gene.name,]$krelFCcategory = newSummary.TRP[newSummary.TRP$gene == gene.name,]$krelFCcategory
}
subset.rateData.TRP$krelFCcategory = factor(subset.rateData.TRP$krelFCcategory, levels=c("krelFCless1", "krelFCmore1"))

ggplot(subset.rateData.TRP) + geom_violin(aes(x=gene, y=krel_FC, fill=krelFCcategory), scale="width") +
        geom_point(aes(x=gene,y=krel_FC,color=krelFCcategory), position="jitter", alpha=0.3 , size=0.16) +
        scale_y_continuous(trans=log2_trans(),
			   breaks=trans_breaks('log2', function(x) 2^x), 
			   labels=trans_format('log2', math_format(.x))) + 
	ylab(expression(log[2] ~ "Foldchange in" ~ italic(K[rel])))  + 
	xlab("Genes") + 
	scale_fill_manual(values = c("#E7298A", "#1B9E77"),
			  labels = c("decrease Krel", "increase Krel")) +   
	scale_color_manual(values = c("#E7298A", "#1B9E77")) + 
	#scale_color_brewer(palette = "Dark2") + 
	#scale_fill_brewer(palette = "Dark2") +
	#labs(fill = "Pause Release Rate", color = "Pause Release Rate") +
	#scale_color_startrek() + scale_fill_startrek()  + 
	#theme_light() + 
	theme_bw() + 
	fontTheme + 
	theme(legend.position="none",axis.text.x = element_text(size=32, face=1.5, angle=30, vjust=0.9, hjust=1),
               axis.title.x = element_blank(), axis.title.y = element_text(size=32), legend.title = element_blank())
ggsave("log2_krel_FC_TRP_distribution.pdf", height=8.4, width=6)
# index of dispersion - krel - always tight because B_fold_change = Krel_foldchange * P_foldchange
newSummary.TRP$IndexDisp_krel = newSummary.TRP$krel_FC_var / newSummary.TRP$krel_FC_mean

ggplot() + 
	geom_violin(data=subset(newSummary.TRP, IndexDisp_krel > 0), aes(x=krelFCcategory, y = IndexDisp_krel, fill=krelFCcategory), scale="width") +
        geom_boxplot(data=subset(newSummary.TRP, IndexDisp_krel > 0), aes(x=krelFCcategory, y = IndexDisp_krel),alpha = 0.1, width=0.36, linewidth=0.3, outlier.shape=NA) +
	geom_point(data=subset(newSummary.TRP, IndexDisp_krel > 0 & gene %in% gene.list), aes(x=krelFCcategory,y=IndexDisp_krel), size=2) +
        geom_text(data=subset(newSummary.TRP, gene %in% gene.list), aes(x=krelFCcategory,y=IndexDisp_krel,label=gene, fontface=2), nudge_x=-0.12, nudge_y=0.3, size=7) +
	scale_y_continuous(trans=log10_trans(), breaks=trans_breaks('log10', function(x) 10^x), labels=trans_format('log10', math_format(.x)), limits = c(NA,1)) +
	ylab(expression(log[10] ~ frac("Variance (Foldchange" ~ italic(K[rel]) ~")", "Mean (Foldchange" ~ italic(K[rel]) ~")"))) +
	#ylab(expression(atop(log[10] ~ "Index of Dispersion", "(valid parameter sets, TRP treatment)")))  + 
	#ylab("log10 index of dispersion - krel") + 
	xlab("Gene Categories") +
	#labs(fill = "Pause Release Rate") + 
	scale_fill_manual(values = c("#E7298A", "#1B9E77")) + 
	#scale_fill_brewer(palette = "Dark2") +
	#scale_fill_startrek() + 
       	theme_bw() + fontTheme + 
       theme(legend.position = "none", 
			     axis.title.y = element_text(size=29),
              axis.title.x = element_blank(),
	     axis.text.x = element_text(size=31)) +
       scale_x_discrete(labels = c(expression(atop("Decrease", "in" ~ italic(K[rel]))), 
				   expression(atop("Increase", "in" ~ italic(K[rel])))))  
       #theme(aspect.ratio = 6/4)

ggsave("log10_krel_FC_TRP_IndexOfDispersion.pdf", width=6, height=8.4)


