# follows section 2.4 of vignette_compartment_modeling.pdf
library(dplyr)
library(parallel)
options(error = function() traceback(4))

source("../R/param_scanner_funcs.R")
#setwd("/home/FCAM/rmukherjee/Summer2022Rotation/R_analysis/paramScanCompartment")
# these data objects were created using paramScannerCompartmentModel.R by setting samplesize = 45000, repl=100, foldquery=15
out <- readRDS("../data/paramScanOut.n45000.repl100.rds")
x <- readRDS("../data/x.rfast99.n45000.repl100.rds")

# start working on pause.density objects (priginally created from callers/run_pauseworkflow.R)
#pdensity.obj.dir = "/home/FCAM/rmukherjee/Summer2022Rotation/R_analysis/pauseWindowsAndBWplots"
pdensity.obj.dir = "../data"
density.90minvs45min.active = readRDS(sprintf("%s/density.90min(activated).v.45min.rds", pdensity.obj.dir))
density.90minvscntrl.active = readRDS(sprintf("%s/density.90min(activated).v.control.rds", pdensity.obj.dir))
density.45minvscntrl.active = readRDS(sprintf("%s/density.45min(activated).v.control.rds", pdensity.obj.dir))


# initial pause index in 45 min case
pause.index.45min.in.90minv45min = mean(filter(density.90minvs45min.active, cond == "45min")$pause.sum, trim=0.05) /  # divide by
                                   mean(filter(density.90minvs45min.active, cond == "45min")$body.avg, trim=0.05)
# average fold change in pause.region
pause.sum.change.bw.90minv45min =  mean(filter(density.90minvs45min.active, cond == "90min (activated)")$pause.sum, trim=0.05) /  # divide by
                                   mean(filter(density.90minvs45min.active, cond == "45min")$pause.sum, trim=0.05)

body.change.bw.90minvs45min = mean(filter(density.90minvs45min.active, cond == "90min (activated)")$body.avg, trim=0.05) /  # divide by
                              mean(filter(density.90minvs45min.active, cond == "45min")$body.avg, trim=0.05)  # divide by

parallel.change.two.param.helper("45min", "90min", "activated", out, x, perc = 0.02,
                       pause.index.45min.in.90minv45min, pause.sum.change.bw.90minv45min, body.change.bw.90minvs45min,
                       foldquery=5.5, samplesize=20000)
# 45 minvscontrol case
pauseindex.cntrl.in.45minvscntrl = mean(filter(density.45minvscntrl.active, cond == "control")$pause.sum, trim=0.05) /  # divide by
                                   mean(filter(density.45minvscntrl.active, cond == "control")$body.avg, trim=0.05) 

pause.sum.change.bw.45minvscntrl =  mean(filter(density.45minvscntrl.active, cond == "45min (activated)")$pause.sum, trim=0.05) /  # divide by
                                   mean(filter(density.45minvscntrl.active, cond == "control")$pause.sum, trim=0.05)

body.change.bw.45minvscntrl = mean(filter(density.45minvscntrl.active, cond == "45min (activated)")$body.avg, trim=0.05) /  # divide by
                              mean(filter(density.45minvscntrl.active, cond == "control")$body.avg, trim=0.05)
parallel.change.two.param.helper("control", "45min", "activated", out, x, perc = 0.02,
                                pauseindex.cntrl.in.45minvscntrl, pause.sum.change.bw.45minvscntrl, body.change.bw.45minvscntrl,
                                foldquery=5.5, samplesize=20000)

# 90min vs control case
pauseindex.cntrl.in.90minvscntrl = mean(filter(density.90minvscntrl.active, cond == "control")$pause.sum, trim=0.05) /  # divide by
                                   mean(filter(density.90minvscntrl.active, cond == "control")$body.avg, trim=0.05) 

pause.sum.change.bw.90minvscntrl =  mean(filter(density.90minvscntrl.active, cond == "90min (activated)")$pause.sum, trim=0.05) /  # divide by
                                   mean(filter(density.90minvscntrl.active, cond == "control")$pause.sum, trim=0.05)

body.change.bw.90minvscntrl = mean(filter(density.90minvscntrl.active, cond == "90min (activated)")$body.avg, trim=0.05) /  # divide by
                              mean(filter(density.90minvscntrl.active, cond == "control")$body.avg, trim=0.05)

parallel.change.two.param.helper("control", "90min", "activated", out, x, perc = 0.02, 
                     pauseindex.cntrl.in.90minvscntrl, pause.sum.change.bw.90minvscntrl, body.change.bw.90minvscntrl,
                     foldquery=5.5, samplesize=20000)

