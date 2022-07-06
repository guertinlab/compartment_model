# this tries to run through param scans to reproduce change in pause region and body desnity between conditions
# reproduced from vignette_compartment_modeling.pdf
library(pksensi)
library(dplyr)
library(deSolve)
library(parallel)

#this is a simplified version, since the gene body densities converge over time
density.one.body <- function(t, initial.state, params = params)  {
    with(as.list(c(params, initial.state)), {
        dP = kinit - (kpre + krel)*P
        dB1 = (krel)*P - kelong*B1
        res = c(dP, dB1)
        list(res) 
    })
}

kinit = 0.1
krel = 0.1
kpre = 0.1 #is it worth excluding this because we cannot distinguish between #contrasting effects of initiation and premature nonproductive release? kelong = 50
kelong = 50

parms = c(kinit, krel, kpre, kelong)
initState = c(P = 1, B1 = 0.01)
t <- seq(from = 0.01, to = 100.01, by = 5) 
y <- ode(y = initState, times = t, func = density.one.body, parms = parms)

params <- c("kinit", "krel", "kpre", "kelong")
# set up param vals and make a uninform distribution from [1/10 * param.val, 10 * param.val] with the exception of kelong which stays between [0.5, 2] of its value
# similar to https://cran.r-project.org/web/packages/pksensi/ vignettes/pbtk1cpt.html
q <- c("qunif", "qunif", "qunif", "qunif")
fold = 15
q.arg <- list(list(min = kinit / fold, max = kinit * fold),
              list(min = krel / fold, max = krel * fold),
              list(min = kpre / fold, max = kpre * fold),
              list(min = kelong/2, max = kelong*2))
set.seed(1234) # seed is being set to tune the random number generator in order to reprroduce the value of x below 
x <- rfast99(params, n = 45000, q = q, q.arg = q.arg, replicate = 100) # random sequences using random phase shift - see ?rfast99 for mathematical references
saveRDS(x, "x.rfast99.n45000.repl100.rds")
print(sprintf("saved %s/x.rfast99.n45000.repl100.rds", getwd()))
outputs <- c("P", "B1")
out <- solve_fun(x, time = t, func = density.one.body,
                 initState = initState, outnames = outputs)
saveRDS(out, file="paramScanOut.n45000.repl100.rds")
print(sprintf("saved %s/paramScanOut.n45000.repl100.rds", getwd()))

