# follows param sensitivity analyses - page 16 onwards of compartmental modeling vignettte
library(deSolve)
#library(dplyr)
library(pksensi)
options(error = function() traceback(4))
#setwd("/home/FCAM/rmukherjee/Summer2022Rotation/R_analysis/paramScanCompartment")
# these data objects were created using paramScannerCompartmentModel.R by setting samplesize = 45000, repl=100, fold=15
out <- readRDS("../data/paramScanOut.n45000.repl100.rds")
x <- readRDS("../data/x.rfast99.n45000.repl100.rds")

pdf('parameter_space_distribution.n45000.repl100.pdf', width=8, height=8, useDingbats=FALSE)
par(mfrow=c(4,4),mar=c(0.8,0.8,0.8,0),oma=c(4,4,2,1), pch =16)
for (j in c("kinit", "kpre", "kelong", "krel")) {
  if ( j == "krel") {
    plot(x$a[,1,j], ylab = "krel", cex = 0.3)
  } else plot(x$a[,1,j], xaxt="n", cex = 0.3, ylab = "")
  for (i in 2:3) {
  if ( j == "krel") {
    plot(x$a[,i,j], ylab = "", yaxt="n", cex = 0.3)
    } else plot(x$a[,i,j], xaxt="n", yaxt="n", cex = 0.3, ylab = "")
  }
  hist <- hist(x$a[,,j], plot=FALSE,
               breaks=seq(from=min(x$a[,,j]), to=max(x$a[,,j]), length.out=20))
  barplot(hist$density, axes=FALSE, space=0, horiz = T, main = j)
}
mtext("Model evaluation", SOUTH<-1, line=2, outer=TRUE)
dev.off()
print('saved parameter_space_distribution.n45000.repl100.pdf')

pdf('pause_density_parameter_sensitivity_index.n45000.repl100.pdf', width=8, height=8, useDingbats=FALSE)
 par(pty="s")
 par(mfrow=c(4,1))
 plot(out, vars = "P")
dev.off()
pdf('body_density_parameter_sensitivity_index.n45000.repl100.pdf', width=8, height=8, useDingbats=FALSE)
par(pty="s")
plot(out, vars = "B1")
dev.off()
#body
pdf('body_density_parameter_sensitivity_steady_state.n45000.repl100.pdf', width=8, height=3, useDingbats=FALSE)
	par(mfcol=c(1,4),mar=c(0.8,2.2,0,0),oma=c(4,4,2,1), pch =16, pty="s") #mfcol(r,c) - plot r x c images 
	for (j in c("kinit", "kpre", "kelong", "krel")){
	    plot(x$a[,1,j], out$y[,1, "100.01","B1"], cex = 0.3, col = '#8c8c8c33',
		 main = paste0("\n", j))
	}
	mtext("Parameter Value", SOUTH<-1, line=0, outer=TRUE)
	mtext("Gene Body Density", WEST<-2, line=0, outer=TRUE)
dev.off()
#pause
pdf('pause_density_parameter_sensitivity_steady_state.n45000.repl100.pdf', width=8, height=3, useDingbats=FALSE)
par(mfcol=c(1,4),mar=c(0.8,2.2,0,0),oma=c(4,4,2,1), pch =16, pty="s") # oma stands for outer.margin, pch is the character symbol index, pty is type of plot rgion - "s" being square, "m" being maximal
  for (j in c("kinit", "kpre", "kelong", "krel")){
    plot(x$a[,1,j], out$y[,1, "100.01","P"], cex = 0.3, col = '#8c8c8c33', #cex- magnification from fefault
    main = paste0("\n", j))
  }
mtext("Parameter Value", SOUTH<-1, line=0, outer=TRUE)
mtext("Pause Density", WEST<-2, line=0, outer=TRUE)
dev.off()

