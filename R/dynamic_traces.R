library(zoo)
library(deSolve)
library(lattice)
source('R/plotting_composites_lattice.R')


############################################################################################
## density.prime
############################################################################################

#' Compartment model adapted from our G&D paper (doi:10.1101/gad.328237.119):
#' Declare differential equations as function
#' dP/dt = kinit - (kpre +krel)p
#' dB/dt = krel * p - kelong * b
#' P is the first dependent variable, promoter density; dP is its derivative wrt time
#' B is the second dependent variable, body density; dB is its derivative wrt time
#' @param t time points to run the ODE model
#' @param initial.state a list containing initial values of the variables being evolved.
#' @param params list of other variables (rates) being used in the model.
#' @return this function is called within ode function, and returns pause body with values in three body compartments.
#' @export
#' @examples
#' library(deSolve)
#' initial.state = c(P = 1, B1 = 0.01)
#' parms = c(kinit = 0.1, krel = 0.1, kpre = 0.1, kelong = 50)
#' t <- seq(from = 0.01, to = 100.01, by = 10)
#' y <- ode(y = initital.state, times = t, func = density.prime, parms = parms)
density.prime <- function(t, initial.state, params = params)  {
   # params contain kinit, kpre, kelong, krel
   # initial state contains c(P, B1, B2, B3)
   with(as.list(c(params, initial.state)), { # with executes code in a new environment with the data fields provided
        dP = kinit - (kpre + krel)*P
        dB1 = (krel)*P - kelong*B1
        dB2 = (kelong)*B1 - kelong*B2
        dB3 = (kelong)*B2 - kelong*B3
        res = c(dP, dB1, dB2, dB3)
        list(res)
    })
}

############################################################################################
## find.body.param
############################################################################################

#' code adapted from our G&D paper (doi:10.1101/gad.328237.119):
#' to visualize it in a composite profile form
#' function for finding the gene body parameter
#' @param bpeak desired gene body level
#' @param tau exponential decay constant
#' @param pausepeak desired paused region level
#' @param min.pk minimal level for gene body peak
#' @param max.pk maximal level for gene body peak
#' @param dpk resolution for the implicit solution
#' @return the desired body parameter
#' @export
#' @examples
#' body.param <- find.body.param(bpeak=0.08, tau=20, pausepeak=0.25, min.pk=0, max.pk=1, dpk=.001)
find.body.param <- function(bpeak=NULL,tau=NULL,
                            pausepeak=NULL,min.pk=0,max.pk=1,
                            dpk=0.001) {
    pk = seq(min.pk,max.pk,dpk) # sequence of peak parameter values
    dat = matrix(0,length(pk),2) # len(pk) x 2 matrix containing bodypeak and body values at each row
    colnames(dat) = c("body_param","body_assymp") 
    for (x in 1:length(pk)) {
        bodypeak = pk[x]
        #print(sprintf("bodypeak: %f", bodypeak))
	#print(bodypeak)
        root = tau*(bodypeak + exp(1)) * exp(-1)
        peak = (root/tau) * exp(-(root - tau)/tau) + bodypeak *
            (1 - exp(-root/tau))
        body = bodypeak * pausepeak / peak
        dat[x,1] = bodypeak
        dat[x,2] = body
    }
# look up the ratio of the body parameter over the assymptote
# use linear interpolation
    inter = findInterval(bpeak,dat[,2]) - 1 
    m = (dat[(inter+1),1] - dat[inter,1]) /
        (dat[(inter+1),2] - dat[inter,2])
    b = dat[inter,1] - m * dat[inter,2]
    est = m * bpeak + b
    return(est)
} 

############################################################################################
## get.pro.waveform
############################################################################################

#' function to get the PRO waveform
#' @param bpeak desired gene body level
#' @param pausepeak desired pause region level
#' @param bodypeak body peak value to obtain a level of bpeak
#' @param bp.seq base pair sequence
#' @param tau exponential decay constant
#' @return a list containing a vector with simulated PRO-seq singal and an object containing the simulated and requested signal 
#' @export
#' @examples
#' x <- get.pro.waveform(pausepeak=0.25, bpeak=0.02, bp.seq=seq(0,999),tau=20)
get.pro.waveform <- function(bpeak=NULL,pausepeak=NULL,bodypeak=NULL,
                             bp.seq=NULL,tau=NULL){
    root = tau*(bodypeak+exp(1))*exp(-1)
    peak = (root/tau) * exp(-(root - tau)/tau) + bodypeak * (1 - exp(-root/tau))
    vec = sapply(bp.seq,function(x){(pausepeak/peak)*((x/tau) * exp(-(x - tau)/tau) +
                                                      bodypeak * (1 - exp(-x/tau)))})
    vec = unname(vec)
    pars = as.data.frame(rbind(cbind(max(vec),
        vec[bp.seq[length(bp.seq)]]),c(pausepeak,bpeak)))
    names(pars) = c("Peak","Assymp")
    rownames(pars) = c("simulated","requested")
    out = list(vec=vec, pars=pars)
    return(out)
}
