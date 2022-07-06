# this tries to run through param scans to reproduce change in pause region and body desnity between conditions
# this tries to implement Section 2.3 of vignette_compartment_modeling.pdf on parallel cores
# follows "GNU MCSim model code" from this webpage: https://nanhung.github.io/pksensi/articles/pbtk1cpt.html 
library(pksensi)
library(dplyr)
library(deSolve)
library(parallel)

#this is a simplified version, since the gene body densities converge over time
# write this ODE function as compartment.model (following "pbtk1cpt.model" in the link above) 
# the file will be compiled and passed onto GNU MCSim code
# Prerequisites
# - run mcsim_install() from pksensi package to install GNU MCSim - https://www.gnu.org/software/mcsim/  - which has Mod program called internally by compile_model() function
# example - mcsim_install(library = "/home/FCAM/rmukherjee/GNUMCsim")
# as administrator priviledges won't be available on most servers, the software request for
# `GNU MCSim` might be separately raised with SysAdmins
## The compartment.model file will look like this (based on density.one.body) function. Please include the "#>" before every line in the .model file.
#> ## compartment.model (Based on Compartment Model, Dr. Michael Guertin Lab) ----
#>
#> States  = { P, B1 };
#> Outputs = { P, B1 };
#>
#> krel = 0.1;   # release from pause
#> kpre = 0.1;   # premature termination
#> kinit = 0.1;  # initialization rate
#> kelong = 50;  # elongation rate
#>
#> Dynamics {
#>   dt (P) = kinit  - ( kpre + krel ) * P;
#>   dt (B1) = krel * P - kelong * B1;
#> }
#> End.
# the init param values
kinit = 0.1
krel = 0.1
kpre = 0.1 #is it worth excluding this because we cannot distinguish between #contrasting effects of initiation and premature nonproductive release? kelong = 50
kelong = 50

mName <- "compartment" # as our model file is named `compartment.model`
compile_model(mName, application = "R") # needs `mod` from GNU MCSim package to create C file.
source(paste0(mName, "_inits.R")) # source the file to change initial conditions
parms <- initParms() # get params for setting them 
parms["kinit"] <- kinit
parms["krel"] <- krel
parms["kpre"] <- kpre
parms["kelong"] <- kelong
initState <- initStates(parms=parms)
initState["P"] <- 1
initState["B1"] <- 0.01
t <- seq(from = 0.01, to = 100.01, by = 5) 
compile_model(mName, application = "mcsim") # compile the model code into an format for execution

#y <- ode(y = initState, times = t, func = density.one.body, parms = parms)

params <- c("kinit", "krel", "kpre", "kelong")
# set up param vals and make a uninform distribution from [1/fold * param.val, fold * param.val] with the exception of kelong which stays between [0.5, 2] of its value
# similar to https://cran.r-project.org/web/packages/pksensi/ vignettes/pbtk1cpt.html
q <- c("qunif", "qunif", "qunif", "qunif")
fold = 20
q.arg <- list(list(min = kinit / fold, max = kinit * fold),
              list(min = krel / fold, max = krel * fold),
              list(min = kpre / fold, max = kpre * fold),
              list(min = kelong/2, max = kelong*2))
set.seed(1234) 
x <- rfast99(params, n = 80000, q = q, q.arg = q.arg, replicate = 100) # random sequences using random phase shift - see ?rfast99 for mathematical references
saveRDS(x, "x.rfast99.n80000.repl100.fold20.rds")
outputs <- c("P", "B1")
conditions <- c("P=1", "B1=0.01")
# TBD: explore parallel option in solve_mcsim() as not sure what it is being set to
out <- solve_mcsim(x, mName = mName, params = params, # params was used to create `x` from`rfast99` 
                               vars = outputs, time = t, 
                               condition = conditions, 
                               parallel = 5) # `parallel=val` is used to call makeCluster(val) from parallel package. Muliple copies of R process are created, so care should be taken in setting it to a reasonable number 
saveRDS(out, file="paramScanOut.n80000.repl100.fold20.rds")
print(sprintf("saved %s/paramScanOut.n80000.repl100.fold20.rds", getwd()))

