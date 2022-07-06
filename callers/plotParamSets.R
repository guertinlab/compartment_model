# the code deals with section 3 of the compartment modeling vignette. (page 32 onwards)
# this file has three sections for three comparisons between 90min, 45min and control.

source("../R/param_scanner_funcs.R")
#setwd("/home/FCAM/rmukherjee/Summer2022Rotation/R_analysis/paramScanCompartment")

# first handle 90min vs 45min case
kinit.krel.lists.45min = readRDS("kinit.krel.lists.activated.45min.v.control.rds")
kinit.kelong.lists.45min = readRDS("kinit.kelong.lists.activated.45min.v.control.rds")
kinit.kpre.lists.45min = readRDS("kinit.kpre.lists.activated.45min.v.control.rds")
kpre.krel.lists.45min = readRDS("kpre.krel.lists.activated.45min.v.control.rds")
kpre.kelong.lists.45min = readRDS("kpre.kelong.lists.activated.45min.v.control.rds")
krel.kelong.lists.45min = readRDS("krel.kelong.lists.activated.45min.v.control.rds")
kinit.vec.45min <- readRDS("kinit.vec.activated.45min.v.control.rds")
krel.vec.45min <- readRDS("krel.vec.activated.45min.v.control.rds")
kpre.vec.45min <- readRDS("kpre.vec.activated.45min.v.control.rds")
kelong.vec.45min <- readRDS("kelong.vec.activated.45min.v.control.rds")

kinit.krel.lattice = two.parameter.bw.plot.lattice(kinit.krel.lists.45min, kinit.vec.45min, krel.vec.45min,
                              params = c("Kinit", "Krel"),
                              constant.params = c('Kelong', 'Kpre'))
kinit.kelong.lattice = two.parameter.bw.plot.lattice(kinit.kelong.lists.45min, kinit.vec.45min, kelong.vec.45min,
                              params = c("Kinit", "Kelong"),
                              constant.params = c('Krel', 'Kpre'))
# for 45min vs control, this created empty results
#kinit.kpre.lattice = two.parameter.bw.plot.lattice(kinit.kpre.lists.45min, kinit.vec.45min, kpre.vec.45min,
#                              params = c("Kinit", "Kpre"),
#                              constant.params = c('Krel', 'Kelong'))
kpre.krel.lattice = two.parameter.bw.plot.lattice(kpre.krel.lists.45min, kpre.vec.45min, krel.vec.45min,
                              params = c("Kpre", "Krel"),
                              constant.params = c('Kinit', 'Kelong'))
kpre.kelong.lattice = two.parameter.bw.plot.lattice(kpre.kelong.lists.45min, kpre.vec.45min, kelong.vec.45min,
                              params = c("Kpre", "Kelong"),
                              constant.params = c('Kinit', 'Krel'))
# krel.kelong had empty results
#krel.kelong.lattice = two.parameter.bw.plot.lattice(krel.kelong.lists.45min, krel.vec.45min, kelong.vec.45min,
#                              params = c("Krel", "Kelong"),
#                              constant.params = c('Kinit', 'Kpre'))

factor.param1.param2.lattice = rbind(kinit.krel.lattice,
                                     kinit.kelong.lattice,
                                     #kinit.kpre.lattice,
                                     kpre.krel.lattice,
                                     kpre.kelong.lattice
                                     #krel.kelong.lattice
                                    )
plot.two.parameter.bw(factor.param1.param2.lattice, factor = 'efferocytosis.45minvcontrol', y.lab = 'Fold Change at activated genes\n45min vs control')

## 90min vs control case
kinit.krel.lists.90min = readRDS("kinit.krel.lists.activated.90min.v.control.rds")
kinit.kelong.lists.90min = readRDS("kinit.kelong.lists.activated.90min.v.control.rds")
kinit.kpre.lists.90min = readRDS("kinit.kpre.lists.activated.90min.v.control.rds")
kpre.krel.lists.90min = readRDS("kpre.krel.lists.activated.90min.v.control.rds")
kpre.kelong.lists.90min = readRDS("kpre.kelong.lists.activated.90min.v.control.rds")
krel.kelong.lists.90min = readRDS("krel.kelong.lists.activated.90min.v.control.rds")
kinit.vec.90min <- readRDS("kinit.vec.activated.90min.v.control.rds")
krel.vec.90min <- readRDS("krel.vec.activated.90min.v.control.rds")
kpre.vec.90min <- readRDS("kpre.vec.activated.90min.v.control.rds")
kelong.vec.90min <- readRDS("kelong.vec.activated.90min.v.control.rds")

kinit.krel.lattice = two.parameter.bw.plot.lattice(kinit.krel.lists.90min, kinit.vec.90min, krel.vec.90min,
                              params = c("Kinit", "Krel"),
                              constant.params = c('Kelong', 'Kpre'))
# not present for 90min vs control 
#kinit.kelong.lattice = two.parameter.bw.plot.lattice(kinit.kelong.lists.90min, kinit.vec.90min, kelong.vec.90min,
#                              params = c("Kinit", "Kelong"),
#                              constant.params = c('Krel', 'Kpre'))
# for 90min vs control, this created empty results
#kinit.kpre.lattice = two.parameter.bw.plot.lattice(kinit.kpre.lists.90min, kinit.vec.90min, kpre.vec.90min,
#                              params = c("Kinit", "Kpre"),
#                              constant.params = c('Krel', 'Kelong'))
kpre.krel.lattice = two.parameter.bw.plot.lattice(kpre.krel.lists.90min, kpre.vec.90min, krel.vec.90min,
                              params = c("Kpre", "Krel"),
                              constant.params = c('Kinit', 'Kelong'))
# not found in 90min vs control case 
#kpre.kelong.lattice = two.parameter.bw.plot.lattice(kpre.kelong.lists.90min, kpre.vec.90min, kelong.vec.90min,
#                              params = c("Kpre", "Kelong"),
#                              constant.params = c('Kinit', 'Krel'))
# krel.kelong had empty results
krel.kelong.lattice = two.parameter.bw.plot.lattice(krel.kelong.lists.90min, krel.vec.90min, kelong.vec.90min,
                              params = c("Krel", "Kelong"),
                              constant.params = c('Kinit', 'Kpre'))

factor.param1.param2.lattice = rbind(kinit.krel.lattice,
                                     #kinit.kelong.lattice,
                                     #kinit.kpre.lattice,
                                     kpre.krel.lattice,
                                     #kpre.kelong.lattice
                                     krel.kelong.lattice
                                    )
plot.two.parameter.bw(factor.param1.param2.lattice, factor = 'efferocytosis.90minvcontrol', y.lab = 'Fold Change at activated genes\n90min vs control')
######### --------- ###########
# 90min vs 45min case

kinit.krel.lists.90minvs45min = readRDS("kinit.krel.lists.activated.90min.v.45min.rds")
kinit.kelong.lists.90minvs45min = readRDS("kinit.kelong.lists.activated.90min.v.45min.rds")
kinit.kpre.lists.90minvs45min = readRDS("kinit.kpre.lists.activated.90min.v.45min.rds")
kpre.krel.lists.90minvs45min = readRDS("kpre.krel.lists.activated.90min.v.45min.rds")
kpre.kelong.lists.90minvs45min = readRDS("kpre.kelong.lists.activated.90min.v.45min.rds")
krel.kelong.lists.90minvs45min = readRDS("krel.kelong.lists.activated.90min.v.45min.rds")
kinit.vec.90minv45min <- readRDS("kinit.vec.activated.90min.v.45min.rds")
krel.vec.90minv45min <- readRDS("krel.vec.activated.90min.v.45min.rds")
kpre.vec.90minv45min <- readRDS("kpre.vec.activated.90min.v.45min.rds")
kelong.vec.90minv45min <- readRDS("kelong.vec.activated.90min.v.45min.rds")

kinit.krel.lattice = two.parameter.bw.plot.lattice(kinit.krel.lists.90minvs45min, kinit.vec.90minv45min, krel.vec.90minv45min,
                              params = c("Kinit", "Krel"),
                              constant.params = c('Kelong', 'Kpre'))
# not found in 90min vs 45min case
#kinit.kelong.lattice = two.parameter.bw.plot.lattice(kinit.kelong.lists.90minvs45min, kinit.vec.90minv45min, kelong.vec.90minv45min,
#                              params = c("Kinit", "Kelong"),
#                              constant.params = c('Krel', 'Kpre'))
# for 45min vs control, 90minvs45min -  this created empty results
#kinit.kpre.lattice = two.parameter.bw.plot.lattice(kinit.kpre.lists.90minvs45min, kinit.vec.90minv45min, kpre.vec.90minv45min,
#                              params = c("Kinit", "Kpre"),
#                              constant.params = c('Krel', 'Kelong'))
kpre.krel.lattice = two.parameter.bw.plot.lattice(kpre.krel.lists.90minvs45min, kpre.vec.90minv45min, krel.vec.90minv45min,
                              params = c("Kpre", "Krel"),
                              constant.params = c('Kinit', 'Kelong'))
# not found in 90min vs 45min case
#kpre.kelong.lattice = two.parameter.bw.plot.lattice(kpre.kelong.lists.90minvs45min, kpre.vec.90minv45min, kelong.vec.90minv45min,
#                              params = c("Kpre", "Kelong"),
#                              constant.params = c('Kinit', 'Krel'))
# krel.kelong had empty results
krel.kelong.lattice = two.parameter.bw.plot.lattice(krel.kelong.lists.90minvs45min, krel.vec.90minv45min, kelong.vec.90minv45min,
                              params = c("Krel", "Kelong"),
                              constant.params = c('Kinit', 'Kpre'))

factor.param1.param2.lattice.90minv45min = rbind(kinit.krel.lattice,
                                     #kinit.kelong.lattice,
                                     #kinit.kpre.lattice,
                                     kpre.krel.lattice,
                                     #kpre.kelong.lattice
                                     krel.kelong.lattice
                                    )

plot.two.parameter.bw(factor.param1.param2.lattice.90minv45min, factor = 'efferocytosis.90minvs45min', y.lab = 'Fold Change at activated genes\n90min vs 45min')


