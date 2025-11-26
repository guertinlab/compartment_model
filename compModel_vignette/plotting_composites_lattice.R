library(lattice)
#source("/Users/rudradeepmukherjee/Documents/UConn Health/comparment.model.project/datasets, genome/R_files/dynamic_traces.R")
source("/Users/rudradeepmukherjee/Documents/UConn Health/comparment.model.project/datasets, genome/dynamic_traces.R")

plot.pro.simulation.composites <- function(input,
                                           filename = 'dynamic_pro_model') {
    colors.ramp = colorRampPalette(c( "pink", "#ff0000", "#cd0000", "#8b0000"),
                                   bias=1, alpha = FALSE)(21)
    pdf(paste0(filename, '.pdf'), width=3.43, height = 9.43)
    print(
        xyplot(input[,2] ~ input[,3], groups = input[,1], data = input, type = 'l',
               scales=list(x=list(cex=0.8,relation = "free"),
                   y =list(cex=0.8,alternating=FALSE)),
               col = colors.ramp,
               par.settings = list(strip.background=list(col="grey85"), 
                                   superpose.symbol = list(pch = c(16),
                                                col = colors.ramp, cex =0.7),
                           superpose.line = list(lwd=c(2.5), col = colors.ramp,
                            lty = c(1,1))),
       aspect=1.0,
#    auto.key = list(points=F, lines=T, cex=0.8),
       lwd=2,
       ylab = list(label = "Simulated PRO-seq Signal", cex =1,font=1),
       xlab = list(label = 'Distance along gene body', cex =1,font=1),
       panel = function(x, y, ...) {
           panel.xyplot(x, y, ...)
           panel.abline(h = 0, lty =1, lwd = 1.5, col = 'grey80')
       }
       )
       )
    dev.off()
}

#combining the previous three functions for a dynamic PRO-seq profile
dynamic.pro.profile <- function(input, tau = 20,min.pk = 0, max.pk =1,
                                dpk =0.001, gene.len = 1000, filename = 'dynamic_pro_model_density') {
    bp = seq(0,gene.len-1)
    colnames(input) = c("time","Pause", "Body") #in context of composite_processes.R, time is 1, Pause is pause.peak and Body is body.peak
    body.parameters = mapply(find.body.param, # the funcrion where multiple param pairs are applied - one at a time i.e; params at row 1 will be applied together, followed by row 2 params and so on. 
                             bpeak = input$Body,# desired gene body level
                             tau = tau,
                             pausepeak = input$Pause,
                             min.pk = min.pk, # min and max level for gene body peak
                             max.pk = max.pk,
                             dpk = dpk)  # the original code say precision of implicit solution 
    input$bodyParam = body.parameters
    df.out = as.data.frame(matrix(NA, nrow=0, ncol = 3))
    bp = seq(0,gene.len-1)
    for (i in 1:length(input$Pause)) {#need to revisit and use apply func.
        x = get.pro.waveform(pausepeak=input$Pause[i],
                             bpeak=input$Body[i],
                             bodypeak=input$bodyParam[i],
                             bp.seq=bp,
                             tau=tau)$vec
        y = as.data.frame(cbind(i, x, bp))
        colnames(y) = colnames(df.out) 
        df.out = rbind(df.out, y)
    }
    plot.pro.simulation.composites(df.out, 
                               filename = filename) 
    return(df.out)
}

# plot gene ody and pause region density based on the input 
plot.changes.wrt <- function(input, filename = 'dynamic_pro_model') {
#    reformat to lattice
    lat.b = cbind(input[,c(1,3)], 'Body') # time and B1 (gene body)
    colnames(lat.b) = c('time', 'signal', 'region')
    lat.p = cbind(input[,c(1:2)], 'Pause') # time and P (promoter pause region)
    colnames(lat.p) = colnames(lat.b)
    lattice.result = rbind(lat.p, lat.b)
  colors.ramp = colorRampPalette(c( "pink", "#ff0000", "#cd0000", "#8b0000"),
                                   bias=1, alpha = FALSE)(21)
#plot
    pdf(paste0(filename, '.pdf'), width=4.43, height = 9.43)
    print(
        xyplot(signal ~ time| region, data = lattice.result, type = c('l', 'p'),
               scales=list(x=list(cex=0.8, relation = "free"),
                   y =list(cex=0.8, relation = "free", alternating=FALSE)),
               col = 'black',
               par.settings = list(strip.background=list(col="grey85"),
                                   superpose.symbol = list(pch = c(16),
                                                           col='black', cex =0.7),
                                   superpose.line = list(col = colors.ramp, lwd=c(2.5),
                                                         lty = c(1,1,1,1,1,1,1,1,1))),
       aspect=1.0,
       between=list(y=0.5, x=0.5),
       lwd=2,
       ylab = list(label = "Density", cex =1,font=1),
       xlab = list(label = 'Time', cex =1,font=1)
       )
)
    dev.off()
}

plot.pro.simulation.composites.2 <- function(input, pause.offset = 0,
					     filename = 'dynamic_pro_model', ylim = c(0, 0.0082),
					     trace.col = c('grey50', '#6aa3ce')) {
    colors.ramp = trace.col
    pdf(paste0(filename, '.pdf'), width=3.43, height = 3.43)
        print(xyplot(input[,2] ~ input[,3], groups = input[,1], data = input, type = 'l',
	            scales=list(x=list(cex=0.8, relation = "free"), 
	                        y=list(cex=0.8, alternating=FALSE, axs = 'i')),
	           xlim = c(0,400),
                   ylim = ylim,
	           col = colors.ramp, 
	           par.settings = list(strip.background = list(col="grey85"),
			               superpose.symbol = list(pch = c(16), col = colors.ramp, cex =0.7),
			               superpose.line = list(lwd=c(2.5), col = colors.ramp, lty = c(1,1))),
			       aspect=1.0,
		               auto.key = list(points=F, lines=T, cex=0.8),
			       lwd=2, ylab = list(label = "Simulated PRO-seq Signal", cex =1,font=1),
			       xlab = list(label = 'Distance from Pause Peak', cex =1,font=1),
			       panel = function(x, y, ...) {
					panel.xyplot(x, y, ...)
				        panel.abline(h = 0, lty =1, lwd = 1.5, col = 'grey80')
					panel.abline(v = 0 - pause.offset, lty =2, lwd = 0.5, col = 'red')
				        panel.abline(v = 50 - pause.offset, lty =2, lwd = 0.5, col = 'red')
					})) 
    dev.off()
   print(sprintf("saved %s/%s.pdf", getwd(),filename)) 
}

pro.integrated.peak <- function(input, tau = 20,min.pk = 0, max.pk =1,
			    dpk =0.001, gene.len = 2000, pause.height, time.char = '0 min') {
    bp = seq(0,gene.len-1)
    colnames(input) = c("iteration","Pause", "Body")
    body.parameters = mapply(find.body.param,
			 #print(body.parameters)
			 bpeak = input$Body,
			 tau = tau,
			 pausepeak = input$Pause,
			 min.pk = min.pk,
			 max.pk = max.pk,
			 dpk = dpk)
    input$bodyParam = body.parameters
    df.out = as.data.frame(matrix(NA, nrow=0, ncol = 4))
    bp = seq(0,gene.len-1)
    for (i in 1:length(input$Pause)) {#need to revisit and use apply func.
       x = get.pro.waveform(pausepeak=input$Pause[i],
				     bpeak=input$Body[i],
				     bodypeak=input$bodyParam[i],
	                             bp.seq=bp,tau=tau)$vec
       y = as.data.frame(cbind(i, x, bp))
	   integrated.intensity = sum(diff(y[1:120, 3])*rollmean(y[1:120, 2],2))
	   out.param = as.data.frame(cbind(input$Pause[i], integrated.intensity,
					    input$bodyParam[i], input$Body[i]))
	   colnames(out.param) = colnames(df.out)
	   df.out = rbind(df.out, out.param)
    }

    pause.height.body.param = df.out[which.min(abs(df.out[,2] - pause.height)),] #return(pause.height.body.param)
    waveform = get.pro.waveform(pausepeak=pause.height.body.param[1,1],
                                bpeak=pause.height.body.param[1,4],
	                        bodypeak=pause.height.body.param[1,3],
	                        bp.seq=bp, tau=tau)$vec
   df.wave = as.data.frame(cbind(time.char, waveform, bp))
   return(df.wave)
}

# use explicit values of kinit, krel, kelong to generate a simulated composite
# _base denotes control/baseline, _expr denotes treatment/experiment.
direct.simulated.plotter.noKpre <- function(kinit_base, krel_base, kinit_expr, krel_expr, kelong,
					    output.plot.name = "model_condnvsbaseline_composite",
					     baseline.condn.name = "control", condn.name = "condn",
					    seq.end = 0.1, by.step = 1e-4, tau = 10, dpk = 1e-5) {
	# plot the steady state profile
	# set dP, dB to zero and solve for peaks
	# dP = kinit - krel*P
	# dB = (krel)*P - kelong*B
	start.pause.peak = kinit_base/krel_base
	start.body.peak  = (krel_base*start.pause.peak)/kelong
	change.pause.peak = kinit_expr/krel_expr # substitute, the changed values
	change.body.peak  = (krel_expr*change.pause.peak)/kelong
	# make body objects
   input.body.baseline = data.frame(c(1:length(seq(0, seq.end, by = by.step))),
  				               seq(0, seq.end, by = by.step), start.body.peak)
   colnames(input.body.baseline) = c('iteration', 'Pause', 'Body')
   input.body.condition = data.frame(c(1:length(seq(0, seq.end, by = by.step))),
  				                seq(0, seq.end, by = by.step), change.body.peak)
   colnames(input.body.condition) = c('iteration', 'Pause', 'Body')
   # TBD: write on pro.integrated.peak
   # you get this error if your parameter space does not fit target data:
   # Error in bodypeak + exp(1) : non-numeric argument to binary operator
   peak.height.baseline = pro.integrated.peak(input.body.baseline,
					    pause.height = start.pause.peak, tau = tau,
					    time.char = baseline.condn.name,
					    dpk = 1e-5)
   peak.height.condition = pro.integrated.peak(input.body.condition,
					    pause.height = change.pause.peak, tau = tau,
					    time.char = condn.name, dpk = 1e-5)

   plot.lattice.model = rbind(peak.height.baseline, peak.height.condition)
   colnames(plot.lattice.model) = c('V1', 'V2', 'V3')
   plot.lattice.model[,2] = as.numeric(as.character(plot.lattice.model[,2]))
   plot.lattice.model[,3] = as.numeric(as.character(plot.lattice.model[,3]))
   # center on pause peak
   pause.off = plot.lattice.model[which.max(plot.lattice.model[,2]),3]
   plot.lattice.model[,3] = plot.lattice.model[,3] - pause.off
   y.end = max(plot.lattice.model[,2])+ max(plot.lattice.model[,2]) * 0.07 # 0.02172083
   #plot
   plot.pro.simulation.composites.2(plot.lattice.model, pause.offset = pause.off,
				   filename = output.plot.name, ylim = c(0, y.end))

}
# use explicit values of kinit, krel, kelong, kpre to generate a simulated composite
# _base denotes control/baseline, _expr denotes treatment/experiment.
direct.simulated.plotter <- function(kinit_base, krel_base, kinit_expr, krel_expr, kelong, kpre,
					    output.plot.name = "model_condnvsbaseline_composite",
					     baseline.condn.name = "control", condn.name = "condn",
					    seq.end = 0.1, by.step = 1e-4, tau = 10, dpk = 1e-5) {
	# plot the steady state profile
	# set dP, dB to zero and solve for peaks
	# dP = kinit - (krel + kpre)*P
	# dB = (krel)*P - kelong*B
	start.pause.peak = kinit_base/(krel_base + kpre)
	start.body.peak  = (krel_base*start.pause.peak)/kelong
	change.pause.peak = kinit_expr/(krel_expr + kpre) # substitute, the changed values
	change.body.peak  = (krel_expr*change.pause.peak)/kelong
	# make body objects
   input.body.baseline = data.frame(c(1:length(seq(0, seq.end, by = by.step))),
  				               seq(0, seq.end, by = by.step), start.body.peak)
   colnames(input.body.baseline) = c('iteration', 'Pause', 'Body')
   input.body.condition = data.frame(c(1:length(seq(0, seq.end, by = by.step))),
  				                seq(0, seq.end, by = by.step), change.body.peak)
   colnames(input.body.condition) = c('iteration', 'Pause', 'Body')
   # TBD: write on pro.integrated.peak
   # you get this error if your parameter space does not fit target data:
   # Error in bodypeak + exp(1) : non-numeric argument to binary operator
   peak.height.baseline = pro.integrated.peak(input.body.baseline,
					    pause.height = start.pause.peak, tau = tau,
					    time.char = baseline.condn.name,
					    dpk = 1e-5)
   peak.height.condition = pro.integrated.peak(input.body.condition,
					    pause.height = change.pause.peak, tau = tau,
					    time.char = condn.name, dpk = 1e-5)

   plot.lattice.model = rbind(peak.height.baseline, peak.height.condition)
   colnames(plot.lattice.model) = c('V1', 'V2', 'V3')
   plot.lattice.model[,2] = as.numeric(as.character(plot.lattice.model[,2]))
   plot.lattice.model[,3] = as.numeric(as.character(plot.lattice.model[,3]))
   # center on pause peak
   pause.off = plot.lattice.model[which.max(plot.lattice.model[,2]),3]
   plot.lattice.model[,3] = plot.lattice.model[,3] - pause.off
   y.end = max(plot.lattice.model[,2])+ max(plot.lattice.model[,2]) * 0.07 # 0.02172083
   #plot
   plot.pro.simulation.composites.2(plot.lattice.model, pause.offset = pause.off,
				   filename = output.plot.name, ylim = c(0, y.end))

}

# this helper function plots data for knit.krel changes and arranges code present from page 40 of compartment modeling vignette
plot.simulated.composites.helper <- function(input.data.identifier="activated.condnvsbaseline", 
					     output.plot.name = "model_condnvsbaseline_composite",
					     baseline.condn.name = "control", condn.name = "condn",
					     seq.end = 0.1, by.step = 1e-4,
					     tau = 10, dpk = 1e-5) {
   print(sprintf("working directory: %s", getwd()))
   kinit.krel.lists = readRDS(sprintf("kinit.krel.lists.%s.rds", input.data.identifier))
   kinit.vec <- readRDS(sprintf("kinit.vec.%s.rds",input.data.identifier))
   krel.vec <- readRDS(sprintf("krel.vec.%s.rds",input.data.identifier))
   kpre.vec <- readRDS(sprintf("kpre.vec.%s.rds",input.data.identifier))
   kelong.vec <- readRDS(sprintf("kelong.vec.%s.rds",input.data.identifier))
   out.pause <- readRDS(sprintf("out.pause.%s.rds",input.data.identifier))
   out.body <- readRDS(sprintf("out.body.%s.rds",input.data.identifier))
   index.kelong.41 = which.min(abs(kelong.vec - 41.66667)) # 2500 bases/per min ~ elongationr rate
   start.kelong = kelong.vec[index.kelong.41]
   start.kinit  = kinit.vec[index.kelong.41]
   start.krel   = krel.vec[index.kelong.41]
   start.kpre   = kpre.vec[index.kelong.41]
   start.body   = out.body[index.kelong.41]
   start.pause  = out.pause[index.kelong.41]
   # for "kinit krel" changes. constants - kelong, kpre
   change.kinit = mean(kinit.krel.lists[[1]][[index.kelong.41]])
   change.krel  = mean(kinit.krel.lists[[2]][[index.kelong.41]])
   print(sprintf('change in Kinit (%s): %f', input.data.identifier, change.kinit/start.kinit))
   print(sprintf('change in Krel (%s): %f', input.data.identifier, change.krel/start.krel))
   krel = start.krel
   kinit = start.kinit
   kelong = start.kelong
   kpre = start.kpre
   # plot the steady state profile
   # set dP, dB to zero and solve for peaks
   # dP = kinit - (kpre + krel)*P
   # dB = (krel)*P - kelong*B
   start.pause.peak = kinit/(kpre + krel)
   start.body.peak  = ((krel)*start.pause.peak)/kelong
   change.pause.peak = change.kinit/(kpre + change.krel) # substitute, the changed values
   change.body.peak  = ((change.krel)*change.pause.peak)/kelong
  
   # make body objects
   input.body.baseline = data.frame(c(1:length(seq(0, seq.end, by = by.step))),
  				               seq(0, seq.end, by = by.step), start.body.peak)
   colnames(input.body.baseline) = c('iteration', 'Pause', 'Body')
   input.body.condition = data.frame(c(1:length(seq(0, seq.end, by = by.step))),
  				                seq(0, seq.end, by = by.step), change.body.peak)
   colnames(input.body.condition) = c('iteration', 'Pause', 'Body')
   # TBD: write on pro.integrated.peak
   # you get this error if your parameter space does not fit target data:
   # Error in bodypeak + exp(1) : non-numeric argument to binary operator
   peak.height.baseline = pro.integrated.peak(input.body.baseline,
					    pause.height = start.pause.peak, tau = tau,
					    time.char = baseline.condn.name,
					    dpk = 1e-5)
   peak.height.condition = pro.integrated.peak(input.body.condition,
					    pause.height = change.pause.peak, tau = tau,
					    time.char = condn.name, dpk = 1e-5)

   plot.lattice.model = rbind(peak.height.baseline, peak.height.condition)
   colnames(plot.lattice.model) = c('V1', 'V2', 'V3')
   plot.lattice.model[,2] = as.numeric(as.character(plot.lattice.model[,2]))
   plot.lattice.model[,3] = as.numeric(as.character(plot.lattice.model[,3]))
   # center on pause peak
   pause.off = plot.lattice.model[which.max(plot.lattice.model[,2]),3]
   plot.lattice.model[,3] = plot.lattice.model[,3] - pause.off
   y.end = max(plot.lattice.model[,2])+ max(plot.lattice.model[,2]) * 0.07 # 0.02172083
   #plot
   plot.pro.simulation.composites.2(plot.lattice.model, pause.offset = pause.off,
				   filename = output.plot.name, ylim = c(0, y.end))
}

