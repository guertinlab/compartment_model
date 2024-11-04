library(bigWig)
library(zoo)
library(lattice)
library(stringr)

############################################################################################
## find.pause.regions
############################################################################################

#' Get start, end coordinates of pause window and gene body for each gene using the bigWig files representing all conditions x replicates. The start and end coordinates are treated as transcription start site (TSS) and transcription termination site (TTS) respectively. The pause window has a 50bp size around the peak signal. The gene body starts from pause window end and ends at TTS. 
#' @param bed.input A bed6 file containing genes of interest, where start and end are transcription start site and transciption termination site respectively.
#' @param combined.plus.bw bigWig file for the plus strand data
#' @param combined.minus.bw bigWig file for the minus strand data
#' @param upwardThreshold
#' @return A bed6 file containing start, end coordinates of pause window (50bp in size), gene body () starting from seventh column.
#' @export
#' @examples
#' master.pTA.coords.file = 'mastercoords.bed' # bed file created from primary Transcript Annotation.
#' pTA <- read.table(master.pTA.coords.file)
#' df.pause.body <- find.pause.regions(pTA, bw.plus, bw.minus) # bw.* files are corresponding bigWig files.
find.pause.regions <- function(bed.input,
                               combined.plus.bw,
                               combined.minus.bw, searchWindow = 200) {
                                        #find TSS of each gene
    bed.input[,7] <- bed.input[,2] # set pause.start to start
    bed.input[bed.input[,6] == '-',7] <- bed.input[bed.input[,6] == '-',3] # set pause.start of '-' strand to end(-) of mastercoords.bed
    bed.input[,8] <- bed.input[,7]  # pause.end
    #bed.input[,8] <- bed.input[,7] + 1   # OLD pause.end

                                        #take searchWindow sized window centered around TSS
    tss.df <- bed.input[,c(1,7,8,4,7,6)]  # TSS in tss.df is now pause.start from bed.input
    # PLUS STRAND
    tss.df[tss.df[,6] == '+', 3] = tss.df[tss.df[,6] == '+',2] + searchWindow # 200bp window (300bp for GRO-seq)
    # MINUS STRAND
    tss.df[tss.df[,6] == '-', 2] = tss.df[tss.df[,6] == '-', 2] - searchWindow # 200bp window (300bp for GRO-seq) 
    # OLD CODE, DON'T USE - thousand base pair window here
    #tss.df[,3] = tss.df[,2] + 500
    #tss.df[,2] = tss.df[,2] - 500
    colnames(tss.df) <- c('chr','start','end','gene','TSS','strand')
    
                                        #find maximum signal within window - this is the pause summit
                                        #read through bigWig vignette to understand probeQuery
    dat.x = bed6.step.probeQuery.bigWig(combined.plus.bw, combined.minus.bw, tss.df, step = 1, abs.value=T,as.matrix = TRUE, op = "sum", follow.strand = FALSE)
    dat.x[is.na(dat.x)] <- 0
    # 50bp window rolling mean,
    dat.x.win = t(apply(dat.x, 1, function(x){rollapply(x, width = 50, FUN = mean, by = 1, by.column = FALSE,align = "left")}))
    index.max = apply(dat.x.win, 1, which.max)
    the.max = apply(dat.x.win, 1, max)

                                        #take 50 bp window around pause summit - this is the pause region
    pauseregion.df = tss.df # this is already 1000bp window
    # use maximum index of plus and minus strands separately
    # set pause start column of + and - strands to existing starts + index.max of bed6.stepQuery
    pauseregion.df[2][pauseregion.df$strand == '-',] = pauseregion.df[2][pauseregion.df$strand == '-',] + index.max[pauseregion.df$strand == '-']
    pauseregion.df[2][pauseregion.df$strand == '+',] = pauseregion.df[2][pauseregion.df$strand == '+',] + index.max[pauseregion.df$strand == '+']
    pauseregion.df[3] = pauseregion.df[2] + 50

                                        #define 'body region' as end of pause region to the end of the gene (strand specific)
    body.df = bed.input[,c(1:4,7,6)]
    colnames(body.df) <- c('chr','start','end','gene','TSS','strand') 
    body.df[body.df$strand == '+',2] = pauseregion.df[pauseregion.df$strand == '+',3]+1 # pause region end + 1
    body.df[body.df$strand == '-',3] = pauseregion.df[pauseregion.df$strand == '-',2]-1 # pause region end - 1

    final.df <- cbind(bed.input[,c(1:4,7,6)],pauseregion.df[,2:3],body.df[,2:3])
    colnames(final.df) <- c('chr','start','end','gene','TSS','strand','pause.start','pause.end','body.start','body.end')
    
    # fix edge cases
    for (i in seq(1:nrow(final.df))) {
         if (final.df[i,]$strand == "-" && final.df[i,]$body.start > final.df[i,]$body.end) {
            final.df[i,]$body.start = final.df[i,]$body.end - 500
         }
        else if (final.df[i,]$strand == "+" && final.df[i,]$body.start > final.df[i,]$body.end) {
            final.df[i,]$body.end = final.df[i,]$body.start + 500 
     }
    }
    return(final.df)
}


############################################################################################
## total.region.density
############################################################################################

#' appends "pause sum" and "body average" columns to a bed6 object created by [find.pause.regions()]. It also appends "pause index" column which is the ratio of "pause sum" and "body average".
#' @param df.input bed6 file created by [find.pause.regions()] containing the start and end of pause window and gene body for genes of interest.
#' @param cond1.plus.bw bigWig file containing data from plus strand in cond1 ("cond1" is baseline.)
#' @param cond1.minus.bw bigWig file containing data from minus strand in cond1 ("cond1" is baseline.)
#' @param cond1 name of the condition (baseline)
#' @param cond2.plus.bw bigWig file containing data from plus strand in cond2 ("cond1" is the condition of interest being compared against the baseline.)
#' @param cond2.minus.bw bigWig file containing data from minus strand in cond2 ("cond1" is the condition of interest being compared against the baseline.)
#' @param cond2 name of the condition of interest
#' @return the input bed6 object appended with columns containing pause sum, body average, pause index for each gene present in the input file.
#' @export
#' @examples
#' see documentation of [find.pause.regions()] to see how bed.object was created.
#' df.input = find.pause.regions(bed.object, bw.plus, bw.minus)
#' df.region.density <- total.region.density(df.input, cond1.plus.bw, cond1.minus.bw, "cond1",
#'                                            cond2.plus.bw, cond2.minus.bw, "cond2")
total.region.density <- function(df.input,
                                 cond1.plus.bw,
                                 cond1.minus.bw,
                                 cond1 = 'cond1',
                                 cond2.plus.bw,
                                 cond2.minus.bw,
                                 cond2 = 'cond2', useColnames=F) {

    #isolate pause and body regions into their own objects
    if (useColnames) {
      pauseregion.df <- df.input[,c("chr","pause.start","pause.end","gene","score","strand")]
      body.df <- df.input[,c("chr", "body.start", "body.end", "gene", "score", "strand")] 
    } else {
      pauseregion.df <- df.input[,c(1,7,8,4,5,6)] 
      body.df <- df.input[,c(1,9,10,4,5,6)]
    }
   
    #calculate polymerase density for pause region (sum) and body region (avg) for each condition
    cond1.df <- df.input
    
    cond1.df$pause.sum = bed6.region.bpQuery.bigWig(cond1.plus.bw,cond1.minus.bw, pauseregion.df, op = "sum", abs.value=T)
    cond1.df$body.avg = bed6.region.bpQuery.bigWig(cond1.plus.bw,cond1.minus.bw, body.df, op = "avg", abs.value=T)
    cond1.df$cond = cond1

    cond2.df <- df.input
    
    cond2.df$pause.sum = bed6.region.bpQuery.bigWig(cond2.plus.bw,cond2.minus.bw, pauseregion.df, op = "sum", abs.value=T)
    cond2.df$body.avg = bed6.region.bpQuery.bigWig(cond2.plus.bw,cond2.minus.bw, body.df, op = "avg", abs.value=T)
    cond2.df$cond = cond2

    #combine objects for each condition
    final.df <- rbind(cond1.df,cond2.df)
    final.df$pause.index = final.df[,'pause.sum'] / final.df[,'body.avg']

    return(final.df)        
}

############################################################################################
## generate.composite.df
############################################################################################

#' creates composite signal starting from -upstream to +downstream based upon bigWig files passed per condtion.
#' @param df.input a bed6 input file - simular to output of [find.pause.regions()] - containing pause.start and pause.end coordinates at seventh and eighth columns. 
#' @param cond1.plus.bw bigWig file for plus strand data in cond1 (baseline).
#' @param cond1.minus.bw bigWig file for plus strand data in cond1 (baseline).
#' @param cond1 name of cond1 (baseline) to be used as to fill cell values under `cond` column in the returned object,
#' @param cond2.plus.bw bigWig file for plus strand data in cond2 (cond2 being the condition of interest being compared against baseline)
#' @param cond2.minus.bw bigWig file for minus strand data in cond2 (cond2 being the condition of interest being compared against baseline)
#' @param cond2 name of the condition which is being compared against baseline
#' @param upstream the number of bases upstream from pause summit (default 350)
#' @param downstream the number of bases downstream from pause summit (default 1400)
#' @param roll,avg the rolling size used in rollmean() to calculate the composites (default 10)
#' @param step a factor used to define the ends of genes by (step * roll.avg)/2 
#input: output dataframe from find.pause.regions and plus / minus .bigWigs for each condition
#output: dataframe with composite polymerase density at each position within window and for each condition
generate.composite.df <- function(df.input,
                                 cond1.plus.bw,
                                 cond1.minus.bw,
                                 cond1 = 'cond1',
                                 cond2.plus.bw,
                                 cond2.minus.bw,
                                 cond2 = 'cond2',
                                 useColnames = F,
				 upstream = 350,
                                 downstream = 1400,
                                 roll.avg = 10, step = 5) {

    #define new dataframe defined by upstream and downstream windows around pause summit
    if (useColnames) {
      print("using column names inside generate.composite.df()")
      window.df <- df.input[,c("chr","pause.start","pause.end","gene","score","strand")]
    } else {
      window.df <- df.input[,c(1,7,8,4,5,6)] 
    }
    #window.df = df.input[,c(1,7,8,4,5,6)]
    window.df[,2] = rowMeans(window.df[,2:3]) # center point between 2 and 3.
    window.df[,3] = window.df[,2] + 1 # end is start + 1
    window.df = fiveprime.bed(window.df, upstreamWindow = upstream, # set start to (start - upstream - 1)
                              downstreamWindow = downstream)  # set end to (end + downstream + 1)

    coordin.start = (-upstream - 1)  + (step * roll.avg)/2
    coordin.end = downstream - (step * roll.avg)/2

    #calculate polymerase density throughout window for each condition
    tss.matrix.cond1 = bed6.step.bpQuery.bigWig(cond1.plus.bw, cond1.minus.bw, window.df,
                                          step = step, abs.value=T, as.matrix=TRUE, follow.strand=TRUE)

    composite.lattice.cond1 = data.frame(seq(coordin.start, coordin.end, by = step),
                                   rollmean(colMeans(tss.matrix.cond1), roll.avg),
                                           cond = cond1,
                                   stringsAsFactors = FALSE)
    colnames(composite.lattice.cond1) = c('x', 'est', 'cond')


    tss.matrix.cond2 = bed6.step.bpQuery.bigWig(cond2.plus.bw, cond2.minus.bw, window.df,
                                                step = step, abs.value=T, as.matrix=TRUE, follow.strand=TRUE)

    composite.lattice.cond2 = data.frame(seq(coordin.start, coordin.end, by = step),
                                   rollmean(colMeans(tss.matrix.cond2), roll.avg),
                                   cond = cond2,
                                   stringsAsFactors = FALSE)
    colnames(composite.lattice.cond2) = c('x', 'est', 'cond')

    composite.df <- rbind(composite.lattice.cond1,composite.lattice.cond2)
    composite.df$x <- as.numeric(composite.df$x)
    return(composite.df)

}

############################################################################################
## plot.composites
############################################################################################

#' uses output of [generate.composite.df()] to generate plot of composite signal.
#' @param dat output of [generate.composite.df()]
#' @param fact name of the factor in the condition being compared against
#' @param comp a small string to describe the condition being made
#' @param summit a string which will be used as a xlabel - distance from <summit>
#' @param y.low the lower limit of the y axis
#' @param y.high the max limit of the y axis. If not provided, it estimates it using the maximum estimate present in input data
#' @param col.lines colors to be used for drawing the composites 
#' @return NULl
plot.composites <- function(dat, fact = 'GR', comp='20.v.0',
                            summit = 'TSS', y.low =0, y.high = NULL, 
                            col.lines = c(rgb(0,0,1,1/2), rgb(1,0,0,1/2),  # rgb takes r,g,b and transparency values
                                          rgb(0.1,0.5,0.05,1/2), rgb(0,0,0,1/2),  
                                          rgb(1/2,0,1/2,1/2), rgb(0,1/2,1/2,1/2), 
                                          rgb(1/2,1/2,0,1/2))) {
    
  if (length(y.high)==0) {
    y.high <- max(dat$est) * 1.05
  }
  print(sprintf("will use %f as ylimit", y.high))
  filename = paste('composite_PolII_density_around_',fact,'_genes_',comp,'.pdf', sep='')
  pdf(filename, width=3.4, 
      height=3.4) 
  print(xyplot(est ~ x, group = cond, data = dat, # plot estimate against x coordinate
               type = 'l', #lin	e
               scales=list(x=list(cex=1,relation = "free"), y =list(cex=1, relation="free",axs = 'i')), # axis scales are free, axs="i" for exact data range
               col = col.lines, # color
               auto.key = list(points=F, lines=T, cex=1.0,font=1), # suitable legend for the grouping variable (legend will be a line) 
               par.settings = list(strip.background=list(col="grey85"), # add settings while plotting. background strip is of a color.
                                   superpose.symbol = list(pch = c(16), # symbols to be used while plotting, pch => plot character
                                                           col=col.lines, cex =0.5), 
                                   superpose.line = list(col = col.lines, lwd=c(3), # lwd is line width, lty is line type
                                                         lty = c(1))),
               par.strip.text=list(cex=0.9, font=1, col='black'), # param to define strp.text appearance
               aspect=1.0, # aspect ratio
               between=list(y=0.5, x=0.5), # space between panels
               lwd=3,
               ylab = list(label = "RNA Polymerase Density", cex =1.0,font=1),
               xlab = list(label = paste("Distance from ", summit, sep=''), cex =1.0,font=1),
               xlim = c(-80,950),
               ylim = c(y.low, y.high),
               panel = function(x, y, ...) {
           panel.xyplot(x, y, ...)
#           panel.abline(h = 0, lty =1, lwd = 1.5, col = 'grey80')
 #          panel.abline(v = -25, lty =2, lwd = 0.5, col = 'red')
 #          panel.abline(v = 25, lty =2, lwd = 0.5, col = 'red')
       }
       )
  )
  dev.off()
  print(sprintf("saved %s", filename))
}

############################################################################################
## plot.bw.plots
############################################################################################

#' creates a plot containing box and whisker plot, violin plot and scatter plot of the change in pause indices between a condition of interest and a baseline.    
#' @param input.df an object containing pause indices for each gene of interest (oer condition,) Ideally, an output of [total.region.density()]
#' @param cond1 name of the baseline condition
#' @param cond2 nqme of the condition which is being compared against baseline
#' @param color color to be used in the violin plot
#' @return NULL 
plot.bw.plots <- function(input.df,
                          cond1 = 'cond1',
                          cond2 = 'cond2',
                          color = 'lightblue') {
    plot.df <- data.frame(genes = unique(input.df$gene),
                          fc = 0)

    for(i in 1:nrow(plot.df)) {
        gene <- plot.df$gene[i]
        plot.df$fc[i] <- input.df[input.df$gene == gene & input.df$cond == cond2,'pause.index'] /
            input.df[input.df$gene == gene & input.df$cond == cond1,'pause.index']
    }
   
    # remove undefined values - Inf, -Inf, NaN 
    #plot.df <- plot.df[!plot.df$fc %in% c(Inf,NaN,-Inf),]
    plot.df <- plot.df[is.finite(plot.df$fc),] 
    plot.df$x <- 1
    
    # calculate ylimits
    log2.fc = log2(plot.df$fc)
    log2.fc = log2.fc[is.finite(log2.fc)]
    log2.ylim.min = min(log2.fc) - 0.03
    log2.ylim.max = max(log2.fc) + 0.03
    print(sprintf("box-whisker-violin: setting min(ylim): %f, max(ylim): %f", log2.ylim.min, log2.ylim.max))
   
    filename = paste0('pause.index.fc.',cond2,'.v.',cond1,'.pdf')
    pdf(filename, width=2.7)

    trellis.par.set(box.umbrella = list(lty = 1, col="#93939380", lwd=2), #whisker plot attributes
                    box.rectangle = list(col = '#93939380', lwd=1.6), # boxplpt attributes
                    plot.symbol = list(col='#93939380', lwd=1.6, pch ='.')) # symbol type

    print(bwplot(log2(fc) ~ x, data = plot.df,
                 between=list(y=1.0, x = 1.0),
                 scales=list(x=list(draw=FALSE),relation="free",rot = 45, alternating=c(1,1,1,1),cex=1,font=1), # rot is rotation, 
                      # alternating is whether labels should alternate between one panel to another
                 xlab = '',
                                        #main = "Pause Index (PI) Ratio",
                 ylab = substitute(paste('log'[2]*'',nn), list(nn=paste0('(',cond2,' PI / ',cond1,' PI)'))), # create an expression using nn
                 horizontal = FALSE, 
                 col= 'black',
                 aspect = 2,
                 #ylim = c(-1.52, 1),
                 ylim = c(log2.ylim.min, log2.ylim.max),
                 par.settings=list(par.xlab.text=list(cex=1.2,font=1),
                                   par.ylab.text=list(cex=1.2,font=1),
                                   par.main.text=list(cex=1.2, font=1),
                                   plot.symbol = list(col='black', lwd=1.6, pch =19, cex = 0.25)),
                 strip = function(..., which.panel, bg) {
                                        #bg.col = c("#ce228e" ,"grey60", "#2290cf","grey90")
                     strip.default(..., which.panel = which.panel,
                                   bg = rep(bg.col, length = which.panel)[which.panel])
                 },
                 panel = function(..., box.ratio, col) {
                     panel.abline(h = 0, col = 'grey45', lty = 2)
                     panel.violin(..., col = color,
                                  varwidth = FALSE, box.ratio = box.ratio, outer = FALSE)
                     panel.stripplot(..., col='#54545380', do.out=FALSE, jitter.data=TRUE, amount = 0.2, pch = 16)
                     panel.bwplot(..., pch = '|', do.out = FALSE)    
                 }))
    dev.off()
    
    print(sprintf("saved %s", filename))
}

############################################################################################
## merge.pbody.deseq2
############################################################################################

#' a helper function to merge two data frames on gene names. ASSUMES that both data frames have rownames set to gene names.  
#' @param df.pause.body a data frame containing start and end coordinate for pause windows and gene body. Ideally, an output of [find.pause.regions()]
#' @param deseq2.df a data object created from running DESeq2 workflow and containing response status of genes.
#' @return an object merged on gene names having columns in order of df.pause.body followed by deseq2.df
merge.pbody.deseq2 <- function(df.pause.body, deseq2.df) {
 # rownames are assumed to be set to gene names
 merged.df = merge(df.pause.body, deseq2.df) 
 merged.df = merged.df[,c(colnames(df.pause.body),colnames(deseq2.df))]
 return(merged.df) 
}

############################################################################################
## merge.pbody.bedfile
############################################################################################

#' a helper function to merge two data frames on gene names. ASSUMES that both data frames have rownames set to gene names.  
#' @param df.pause.body a data frame containing start and end coordinate for pause windows and gene body. Ideally, an output of [find.pause.regions()]
#' @param bed.df a data object containing coordinates of genes in a particular condition (activated/repressed). The column names are "chr", "bed.body.start", "bed.body.end", "gene", "type", "strand". The type can be - activated or repressed or any other values
#' @return an object merged on gene names having columns in order of df.pause.body followed by that of b.df
merge.pbody.bedfile <- function(df.pause.body, bed.df) {
 # rownames are assumed to be set to gene names
 merged.df = merge(df.pause.body, bed.df, by="row.names", all=T) 
 
 merged.df = merged.df[,c(colnames(df.pause.body),colnames(bed.df))]
 return(merged.df) 
}


############################################################################################
## run.plotting.steps
############################################################################################

#' run the plotting steps given a pause.body object (with desired genes) and pair of bigWigs from two conditions being compared
#' @param df.pausebody an object containing start, end coordinates of pause window and gene body.
#' @param cond1.name name of baseline condition
#' @param bw.cond1.plus bigWig file containing plus strand data from baseline condition
#' @param bw,cond1.minus bigWig file containing minus strand data from baseline condition
#' @param cond2.name name of the condition which is being compared against the baseline condition
#' @param bw.cond2.plus bigWig file containing plus strand data from the condition which is being compared against baseline condition
#' @param bw.cond2.minus bigWig file containing minus strand data from the condition which is being compared against baseline condition'
#' @param factor.name descriptive name of the factor involved in the comparison of condition to baseline
#' @param color.names to be used as colors for plotting. The order of colors chosen will be same as the alphabetical order of names of two conditions
#' @param fileidentifier to be used as substring in filenames while creating various objects 
#' @return NULL
run.plotting.steps <- function(df.pausebody, cond1.name, bw.cond1.plus, bw.cond1.minus, cond2.name, bw.cond2.plus, bw.cond2.minus, 
			       factor.name="efferocytosis",
			       fileidentifier = NULL, useColnames = F,
			       color.names=c(rgb(0,0,1,1/2), rgb(1,0,0,1/2))) {
 df.composite = generate.composite.df(df.pausebody,
                                   bw.cond1.plus,
                                   bw.cond1.minus,
                                   cond1.name,
                                   bw.cond2.plus,
                                   bw.cond2.minus,
                                   cond2.name, useColnames)
 # plot composites
 trim.cond1.name = str_replace_all(str_trim(cond1.name), " ", "")
 trim.cond2.name = str_replace_all(str_trim(cond2.name), " ", "") 
 if (is.null(fileidentifier)) {
   fileidentifier = sprintf("%s.v.%s", trim.cond2.name, trim.cond1.name)
 }
#df.composite <- readRDS(sprintf("composite.%s.rds", fileidentifier))
 saveRDS(df.composite, file=sprintf("composite.%s.rds", fileidentifier))
 print(sprintf("saved composite.%s.rds", fileidentifier))
 plot.composites(dat = df.composite,
                fact = factor.name,
                comp = fileidentifier,
                summit = 'Pause Peak', col.lines=color.names)
 
 # plot box and whisker for pause indices
 df.density <- total.region.density(df.pausebody,
                                   bw.cond1.plus,
                                   bw.cond1.minus,
                                   cond1.name,
                                   bw.cond2.plus,
                                   bw.cond2.minus,
                                   cond2.name, useColnames)
 saveRDS(df.density, file=sprintf("density.%s.rds", fileidentifier))
 print(sprintf("saved density.%s.rds", fileidentifier))
 plot.bw.plots(input.df = df.density,
              cond1 = cond1.name,
              cond2 = cond2.name,
              color = 'orange')
}

