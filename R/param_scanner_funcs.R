# functions reproduced from compartment_modeling_vignette.R
library(deSolve)
library(pksensi)
library(dplyr)
library(parallel)
library(lattice)

############################################################################################
## change.two.parameters
############################################################################################

#' Takes sets of two params to be changed and two param sets to be held constant, and selects all the param sets (being changed) which satisfy the target pause and body levels 
#' @param func.density a function which defines the compartment model and accomodates two of the param sets held constant and two params which are changed.
#' @param p.change change in pause density
#' @param b.change change in body density
#' @param baseline.pause a vector containing baseline pause densities
#' @param baseline.body a vector containing baseline body densities
#' @param param1 a vector containing individual realisations of first param to be changed
#' @param param2 a vector containing individual realisations of second param to be changed
#' @param params a string vector naming the two params being changed
#' @param constant.param1.vec a vector containing realisations of the first param to be kept constant
#' @param constant.param2.vec a vector containing realisations of the second param to be kept constant
#' @param perc percent threshold for matching targeted and realized values
#' @param fold.query the params to be changed will be sweeped from [1 / fold.query, fold.query] * original param vals
#' @param sample.size the number of samples randomly choosen from the sample space
#' @return a list containing two sets of modified parameters, i.e; param1 and param2 sets which satisify change in pause region (baseline.pause*p.change) and gene body (baseline.body*b.change) in the perc threshold provided
#' @export
#' @examples
#' #this function can be called separately. It is invoked within parallel.change.two.param.helper() where details are hidden from the user
#' change.two.parameters(func.density = density.prime.kpre.kelong,
#'                                  p.change = 1.15,
#'                                  b.change = 0.95,
#'                                  baseline.pause = out.pause,
#'                                  baseline.body = out.body,
#'                                  param1 = kpre.vec,
#'                                  param2 = kelong.vec,
#'                                  params = c("kpre", "kelong"),
#'                                  constant.param1.vec = kinit.vec,
#'                                  constant.param2.vec = krel.vec,
#'                                  perc = 0.05, fold.query = 3, sample.size = 10000)
change.two.parameters <- function(func.density = density.prime.kinit.krel, p.change = pause.change,
                         b.change = body.change, baseline.pause = out.pause, baseline.body = out.body,
                         param1 = kinit.vec, param2 = krel.vec, params = c("kinit", "krel"), constant.param1.vec = kpre.vec,
                         constant.param2.vec = kelong.vec, perc = 0.05, fold.query = 5, sample.size = 10000, ... ) {
	list.out.param1.values = list()
	list.out.param2.values = list()
        for (i in 1:length(param1)) { # start going through each param pair (param1[i] and param2[i]) one by one
                parms = eval(parse(text=sprintf('c("%s" = param1[i], "%s" = param2[i])',params[1], params[2])))    
		initState = c(P = 1, B1 = 0.01)
		t <- seq(from = 0.01, to = 100.01, by = 100)
                # ode is called with constant params as constant.param1[i], constant.param2[i]
		y <- ode(y = initState, times = t, func = func.density, parms = parms,
		     constant.param1 = constant.param1.vec[i], constant.param2 = constant.param2.vec[i])
                q <- c("qunif", "qunif")
		#allow rates to change by five fold in both directions (make fold change an option)
		q.arg <- list(list(min = parms[1] / fold.query, max = parms[1] * fold.query),
                              list(min = parms[2] / fold.query, max = parms[2] * fold.query))
		set.seed(1234)
		x <- rfast99(params, n = sample.size, q = q, q.arg = q.arg, replicate = 5)
		outputs <- c("P", "B1")
		out <- solve_fun(x, time = t, func = func.density,
		              constant.param1 = constant.param1.vec[i], constant.param2 = constant.param2.vec[i],
		              initState = initState, outnames = outputs) # what parameter values fit these densities?
		# allow the target values to differ by variable perc (default of 5%)
		indices.of.interest = out$y[,1, "100.01","B1"] > (baseline.body[i] * b.change * (1 - perc)) & 
                                      out$y[,1, "100.01","B1"] < (baseline.body[i] * b.change * (1 + perc)) &
		                      out$y[,1, "100.01","P"]  > (baseline.pause[i] * p.change * (1 - perc)) &
		                      out$y[,1, "100.01","P"]  < (baseline.pause[i] * p.change * (1 + perc))
                
                new.param1 = x$a[,1,params[1]][indices.of.interest]
		new.param2 = x$a[,1,params[2]][indices.of.interest]
		# in all cases the fold change in rate of pause release is # greater than the change in initiation
		list.out.param1.values[[i]] = new.param1
                list.out.param2.values[[i]] = new.param2
               # print(baseline.pause[i])
               # print(baseline.body[i])
	       # print(params[1])
               # print(param1[i]) 
               # print(new.param1) 
               # print(params[2]) 
               # print(param2[i])
               # print(new.param2)
	}
	return(list(list.out.param1.values, list.out.param2.values)) 

}

############################################################################################
## density.prime.kinit.krel
############################################################################################

#' an implementation of compartment model where kinit and krel are being modified, and kpre and kelong are held constant. (To be passed to ode() solver in package deSolve(). ) 
#' @param t series of time points where this function is used to evaluate the change in variables
#' @param initial.state the values of variables being used to start the model evolution
#' @param params containing values of rates being changed, i.e;  kinit and krel
#' @param constant.param1 the first param not being modified, here kpre
#' @param constant.param2 the second param not being modified, here kelong
#' @return two values corresponding to pause and body densities 
density.prime.kinit.krel <- function(t, initial.state, params = params,
                                     constant.param1, constant.param2)  {
    with(as.list(c(params, initial.state)), {
        dP = kinit - (constant.param1 + krel)*P
        dB1 = (krel)*P - constant.param2*B1
        res = c(dP, dB1)
        list(res) })}

############################################################################################
## density.prime.kinit.kelong
############################################################################################

#' an implementation of compartment model where kinit and kelong are being modified, and kpre and krel are held constant. (To be passed to ode() solver in package deSolve(). ) 
#' @param t series of time points where this function is used to evaluate the change in variables
#' @param initial.state the values of variables being used to start the model evolution
#' @param params containing values of rates being changed, i.e;  kinit and kelong
#' @param constant.param1 the first param not being modified, here kpre
#' @param constant.param2 the second param not being modified, here krel
#' @return two values corresponding to pause and body densities
density.prime.kinit.kelong <- function(t, initial.state, params = params,
                                       constant.param1, constant.param2)  {
    with(as.list(c(params, initial.state)), {
        dP = kinit - (constant.param1 + constant.param2)*P
        dB1 = (constant.param2)*P - kelong*B1
        res = c(dP, dB1)
        list(res)
    })}

############################################################################################
## density.prime.kinit.kpre
############################################################################################

#' an implementation of compartment model where kinit and kpre are being modified, and krel and kelong are held constant. (To be passed to ode() solver in package deSolve(). ) 
#' @param t series of time points where this function is used to evaluate the change in variables
#' @param initial.state the values of variables being used to start the model evolution
#' @param params containing values of rates being changed, i.e;  kinit and kpre
#' @param constant.param1 the first param not being modified, here krel
#' @param constant.param2 the second param not being modified, here kelong
#' @return two values corresponding to pause and body densities 
density.prime.kinit.kpre <- function(t, initial.state, params = params,
                                     constant.param1, constant.param2)  {
    with(as.list(c(params, initial.state)), {
        dP = kinit - (kpre + constant.param1)*P
        dB1 = (constant.param1)*P - constant.param2*B1
        res = c(dP, dB1)
        list(res)
}) }

############################################################################################
## density.prime.kpre.krel
############################################################################################

#' an implementation of compartment model where kpre and krel are being modified, and kinit and kelong are held constant. (To be passed to ode() solver in package deSolve(). ) 
#' @param t series of time points where this function is used to evaluate the change in variables
#' @param initial.state the values of variables being used to start the model evolution
#' @param params containing values of rates being changed, i.e;  kpre and krel
#' @param constant.param1 the first param not being modified, here kinit
#' @param constant.param2 the second param not being modified, here kelong
#' @return two values corresponding to pause and body densities 
density.prime.kpre.krel <- function(t, initial.state, params = params,
                                    constant.param1, constant.param2)  {
    with(as.list(c(params, initial.state)), {
        dP = constant.param1 - (kpre + krel)*P
        dB1 = (krel)*P - constant.param2*B1
        res = c(dP, dB1)
        list(res) })
}

############################################################################################
## density.prime.kpre.kelong
############################################################################################

#' an implementation of compartment model where kpre and kelong are being modified, and kinit and krel are held constant. (To be passed to ode() solver in package deSolve(). ) 
#' @param t series of time points where this function is used to evaluate the change in variables
#' @param initial.state the values of variables being used to start the model evolution
#' @param params containing values of rates being changed, i.e;  kinit and krel
#' @param constant.param1 the first param not being modified, here kinit
#' @param constant.param2 the second param not being modified, here krel
#' @return two values corresponding to pause and body densities
density.prime.kpre.kelong <- function(t, initial.state, params = params,
                                      constant.param1, constant.param2)  {
    with(as.list(c(params, initial.state)), {
        dP = constant.param1 - (kpre + constant.param2)*P
        dB1 = (constant.param2)*P - kelong*B1
        res = c(dP, dB1)
        list(res)
}) }

############################################################################################
## density.prime.krel.kelong
############################################################################################

#' an implementation of compartment model where krel and kelong are being modified, and kinit and kpre are held constant. (To be passed to ode() solver in package deSolve(). ) 
#' @param t series of time points where this function is used to evaluate the change in variables
#' @param initial.state the values of variables being used to start the model evolution
#' @param params containing values of rates being changed, i.e;  krel and kelong
#' @param constant.param1 the first param not being modified, here kinit
#' @param constant.param2 the second param not being modified, here kpre
#' @return two values corresponding to pause and body densities
density.prime.krel.kelong <- function(t, initial.state, params = params,
                                      constant.param1, constant.param2)  {
    with(as.list(c(params, initial.state)), {
        dP = constant.param1 - (constant.param2 + krel)*P
        dB1 = (krel)*P - kelong*B1
        res = c(dP, dB1)
        list(res)
}) }

############################################################################################
## parallel.change.two.param.helper
############################################################################################

#' a general function to call [change.two.parameters()] given the condition name and the objects associated with scanned parameters
#' @param cond1 name of baseline conditon
#' @param cond2 name of condition which is being compared with the baseline
#' @param direction nature of the change in the data being compared between condition of interest and the baseline. Say, "activated"
#' @param output.param.scan an object containing the model output from param scans.
#' @param x.fast99.obj an object created from using rfast99() function from pksensi package
#' @param perc percentage threshold used while matching simulated values with the target values
#' @param pause.index.baseline the pause index in the baseline condition
#' @param pause.sum.change the change in pause sum between conditon of interest and baseline
#' @param body.change the change in body signal between condition of interest and baseline
#' @param foldquery the param scan will be done with values ranging between [1/foldquery, foldquery] * original param values
#' @param samplesize the number of samples randomly selected from the param space defined by foldquery.
#' @export
#' @examples
#' # see paramFittingCompartmentModel.R for a description of calling this function
parallel.change.two.param.helper <- function(cond1="cond1", cond2="cond2", direction="activated",
                                 output.param.scan, x.fast99.obj, perc = 0.02, 
                                 pause.index.baseline, pause.sum.change, body.change,
                                 foldquery = 2, samplesize = 1000) {
 fileidentifier = sprintf("%s.%s.v.%s", direction, cond2, cond1)
 # just substituting variables for readable expressions
 out = output.param.scan
 pre.ratio = pause.index.baseline
 x = x.fast99.obj

 indices.of.interest = (( out$y[,1,10.1,"P"] / (out$y[,1,10.1,"B1"])) > (pre.ratio*(1 - perc))) &
                       ((out$y[,1,10.1,"P"] / (out$y[,1,10.1,"B1"])) < (pre.ratio * (1 + perc)))
 #--
 filename = sprintf("indices.of.interest.%s.rds", fileidentifier)
 saveRDS(indices.of.interest, file=filename)
 print(sprintf("saved %s/%s",getwd(), filename))
 #--
 kinit.vec = x$a[,1,"kinit"][indices.of.interest]
 krel.vec = x$a[,1,"krel"][indices.of.interest]
 kelong.vec = x$a[,1,"kelong"][indices.of.interest]
 kpre.vec = x$a[,1,"kpre"][indices.of.interest]
 out.pause = out$y[,1,10.1,"P"][indices.of.interest]
 out.body = out$y[,1,10.1,"B1"][indices.of.interest] 
 # indices.of.interest can be directly figured out by adding this condition there, 
 # but following the code to maintain sync with the vignette whch explores some issues on page 23 and 24
 pause.constraint.indices = out.pause <= 2
 kinit.vec = kinit.vec[pause.constraint.indices]
 krel.vec = krel.vec[pause.constraint.indices]
 kelong.vec = kelong.vec[pause.constraint.indices]
 kpre.vec = kpre.vec[pause.constraint.indices]
 out.pause = out.pause[pause.constraint.indices]
 out.body = out.body[pause.constraint.indices] 
 # save the vectors so that they can be retrieved later
 saveRDS(kinit.vec, file=sprintf("kinit.vec.%s.rds", fileidentifier))
 saveRDS(krel.vec, file=sprintf("krel.vec.%s.rds", fileidentifier))
 saveRDS(kelong.vec, file=sprintf("kelong.vec.%s.rds", fileidentifier))
 saveRDS(kpre.vec, file=sprintf("kpre.vec.%s.rds", fileidentifier))
 saveRDS(out.pause, file=sprintf("out.pause.%s.rds", fileidentifier))
 saveRDS(out.body,  file=sprintf("out.body.%s.rds", fileidentifier))


 kinit.krel.lists <- mcparallel(change.two.parameters(func.density = density.prime.kinit.krel,
                                  p.change = pause.sum.change,
                                  b.change = body.change,
                                  baseline.pause = out.pause,
                                  baseline.body = out.body,
                                  param1 = kinit.vec,
                                  param2 = krel.vec,
                                  params = c("kinit", "krel"),
                                  constant.param1.vec = kpre.vec,
                                  constant.param2.vec = kelong.vec,
                                  perc = 0.05, fold.query = foldquery, sample.size = samplesize),
                                name="kinit.krel") # name to be associated with parallel job
# filename = sprintf('kinit.krel.lists.%s.rds', fileidentifier)
# saveRDS(kinit.krel.lists, file = filename)
# print(sprintf("saved %s/%s", getwd(), filename))

 kinit.kelong.lists = mcparallel(change.two.parameters(func.density = density.prime.kinit.kelong,
                                  p.change = pause.sum.change,
                                  b.change = body.change,
                                  baseline.pause = out.pause,
                                  baseline.body = out.body,
                                  param1 = kinit.vec,
                                  param2 = kelong.vec,
                                  params = c("kinit", "kelong"),
                                  constant.param1.vec = kpre.vec,
                                  constant.param2.vec = krel.vec,
                                  perc = 0.05, fold.query = foldquery, sample.size = samplesize),
                                 name="kinit.kelong") # name to be associated with parallel job
 #filename = sprintf('kinit.kelong.lists.%s.rds', fileidentifier)
 #saveRDS(kinit.kelong.lists, file = filename)
 #print(sprintf("saved %s/%s", getwd(), filename))
 
 kinit.kpre.lists = mcparallel(change.two.parameters(func.density = density.prime.kinit.kpre,
                                  p.change = pause.sum.change,
                                  b.change = body.change,
                                  baseline.pause = out.pause,
                                  baseline.body = out.body,
                                  param1 = kinit.vec,
                                  param2 = kpre.vec,
                                  params = c("kinit", "kpre"),
                                  constant.param1.vec = krel.vec,
                                  constant.param2.vec = kelong.vec,
                                  perc = 0.05, fold.query = foldquery, sample.size = samplesize),
                               name="kinit.kpre") # name to be associated with parallel job
 #filename = sprintf('kinit.kpre.lists.%s.rds', fileidentifier)
 #saveRDS(kinit.kpre.lists, file = filename)
 #print(sprintf("saved %s/%s", getwd(), filename))
 kpre.krel.lists = mcparallel(change.two.parameters(func.density = density.prime.kpre.krel,
                                  p.change = pause.sum.change,
                                  b.change = body.change,
                                  baseline.pause = out.pause,
                                  baseline.body = out.body,
                                  param1 = kpre.vec,
                                  param2 = krel.vec,
                                  params = c("kpre", "krel"),
                                  constant.param1.vec = kinit.vec,
                                  constant.param2.vec = kelong.vec,
                                  perc = 0.05, fold.query = foldquery, sample.size = samplesize),
                               name="kpre.krel") # name to be associated with parallel job

 #filename = sprintf('kpre.krel.lists.%s.rds', fileidentifier)
 #saveRDS(kpre.krel.lists, file = filename)
 #print(sprintf("saved %s/%s", getwd(), filename))
 
 kpre.kelong.lists = mcparallel(change.two.parameters(func.density = density.prime.kpre.kelong,
                                  p.change = pause.sum.change,
                                  b.change = body.change,
                                  baseline.pause = out.pause,
                                  baseline.body = out.body,
                                  param1 = kpre.vec,
                                  param2 = kelong.vec,
                                  params = c("kpre", "kelong"),
                                  constant.param1.vec = kinit.vec,
                                  constant.param2.vec = krel.vec,
                                  perc = 0.05, fold.query = foldquery, sample.size = samplesize),
                               name="kpre.kelong") # name to be associated with parallel job

 #filename = sprintf('kpre.kelong.lists.%s.rds', fileidentifier)
 #saveRDS(kpre.kelong.lists, file = filename)
 #print(sprintf("saved %s/%s", getwd(), filename))
 krel.kelong.lists = mcparallel(change.two.parameters(func.density = density.prime.krel.kelong,
                                  p.change = pause.sum.change,
                                  b.change = body.change,
                                  baseline.pause = out.pause,
                                  baseline.body = out.body,
                                  param1 = krel.vec,
                                  param2 = kelong.vec,
                                  params = c("krel", "kelong"),
                                  constant.param1.vec = kinit.vec,
                                  constant.param2.vec = kpre.vec,
                                  perc = 0.05, fold.query = foldquery, sample.size = samplesize),
                               name="krel.kelong") # name to be associated with parallel job

 #filename = sprintf('krel.kelong.lists.%s.rds', fileidentifier)
 #saveRDS(krel.kelong.lists, file = filename)
 #print(sprintf("saved %s/%s", getwd(), filename))
 res <- mccollect(list(kinit.krel.lists, kinit.kelong.lists, kinit.kpre.lists,
                  kpre.krel.lists, kpre.kelong.lists, krel.kelong.lists))
 saveRDS(res, file=sprintf("mccollect.res.%s.rds", fileidentifier))
 print(sprintf("saved mccollect.res.%s.RDS", fileidentifier))
 listnames = c("kinit.krel", "kinit.kelong", "kinit.kpre", "kpre.krel", "kpre.kelong", "krel.kelong")
 # save results of parallel run - not necessary as mcccollect obejct is being saved, but maintaining sync
 # with the vignette. 
 for (i in 1:length(res))  {
  filename = sprintf('%s.lists.%s.rds', listnames[i], fileidentifier)
  saveRDS(eval(parse(text=sprintf("res[i]$%s",listnames[i]))), file = filename)
  print(sprintf("saved %s/%s", getwd(), filename))
 }

}

############################################################################################
## two.parameter.bw.plot.lattice
############################################################################################

#' calculates fold change between given param list and original values and returns a data frame
#' @param param1.param2.list a list containing the derived param values from running [change.two.parameters()] 
#' @param param1.vec original values of param1
#' @param param2.vec original values of param2
#' @param params names of params which were modified i.e; the identity of param1, param2
#' @param constant.params names of parameters which were held constant
#' @return data frame containing three columns - fold.change, (modified) rates and description of constant params. The description in 3rd column is constant_<param1>_<param2>
#' @export
#' @examples
#' #the function was taken from page 28 of compartment modeling vignette 
#' #see ./callers/plotParamSets.R for how these functions were called. The "two param" lists were created by parallel.change.two.param.helper()
two.parameter.bw.plot.lattice <- function(param1.param2.lists, param1.vec,
                                          param2.vec, params = c("Kinit", "Krel"),
                                          constant.params = c('Kelong', 'Kpre')) {
    vec.change.param1 = c()
    vec.change.param2 = c()
    count = 0
    for (i in 1:length(param1.param2.lists[[1]])) {
        if(length(param1.param2.lists[[1]][[i]]) != 0){
        count = count + 1
        vec.change.param1[count] = mean(param1.param2.lists[[1]][[i]])/param1.vec[i]
        vec.change.param2[count] = mean(param1.param2.lists[[2]][[i]])/param2.vec[i]
        }
    }
    print(sprintf("found %d sets of params.", count))
    factor.param1.param2 = data.frame(rbind(cbind(vec.change.param1, params[1]),
                                   cbind(vec.change.param2, params[2])), stringsAsFactors = FALSE)
    factor.param1.param2[,3] = paste0('constant_', constant.params[1], '_', constant.params[2])
    colnames(factor.param1.param2) = c('fold.change', 'rate', 'constants')
    factor.param1.param2[,1] = as.numeric(factor.param1.param2[,1])
    return(factor.param1.param2)
}

#####################################################
## plot.two.parameter.bw
#####################################################

#' plots fold change for new parameters (obtained from param scan) w.r.t to original parameters.
#' @param factor.param1.param2 an object created from [two.parameter.bw.plot.lattice()] containing three columns fold.change, (modified) rates and description of constant params.
#' @param factor name of the factor/experiment which will be used in the filnename as <factor>_change_two_parameters.pdf
#' @param y.lab label for y-axis
#' @param y.lim two values giving the lower and upper limit of y-axis. If not provided, the upper limit is calculated using the data and lower limit is set to zero.
#' @return NULL
plot.two.parameter.bw <- function(factor.param1.param2, factor = "GR", y.lab, y.lim = NULL) {
  filename = paste0(factor ,'_change_two_parameters.pdf')
  if (is.null(y.lim)) {
    y.lim = c(0, max(factor.param1.param2$fold.change)+1.2)
  }
  pdf(filename,
        useDingbats = FALSE, width=7,
        height=3.2* ceiling((length(unique(factor.param1.param2$constants))/3))) 
  trellis.par.set(box.umbrella = list(lty = 1, col="#93939380", lwd=1),
                    box.rectangle = list(col = '#93939380', lwd=1),
                    plot.symbol = list(col='#93939380', lwd=1, pch ='.'))
  print(bwplot(fold.change ~ rate | constants, data = factor.param1.param2,
        between=list(y=1.0, x = 1.0),
        scales=list(x=list(draw=TRUE, rot =45),
                    #y = list(axs = 'i'),
                    relation="free",
        alternating=c(1,1,1,1),cex=1,font=1),
        xlab = 'rate', ylab = y.lab, horizontal =FALSE, col= 'black',
        aspect = 1.35,
        #ylim = y.lim,
        box.width = 0.7,
        drop.unused.levels = TRUE,
        #yout = c(3,2),
        #index.cond = list(c(1,2,3,4,5)),
        par.settings=list(par.xlab.text=list(cex=1.2,font=1),
                  par.ylab.text=list(cex=1.2,font=1),
                  par.main.text=list(cex=0.7, font=1),
                  plot.symbol = list(col='black', lwd=1.6, pch =19, cex = 0.25)),
        strip = function(..., which.panel, bg) {
                  bg.col = c("grey85")
                  strip.default(..., which.panel = which.panel,
                  bg = rep(bg.col, length = which.panel)[which.panel])
        },
        panel = function(..., box.ratio, col) {
            #panel.abline(h = 0, col = 'grey45', lty = 2)
            #panel.violin(..., col = '#736fff',
            #varwidth = FALSE, box.ratio = box.ratio, outer = FALSE) 
            panel.stripplot(..., col='#54545380', do.out=FALSE, jitter.data=TRUE,
                    amount = 0.2, pch = 16)
            panel.bwplot(..., pch = '|', do.out = FALSE)
   }))
 dev.off()
 print(sprintf('saved %s', filename))
}

