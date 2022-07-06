# we will start with 45min vs control, and try to see if they can be turned into functions
#source("/home/FCAM/rmukherjee/Summer2022Rotation/R_files/param_scanner_funcs.R")
source("../R/plotting_composites_lattice.R")
source("../R/dynamic_traces.R")
#setwd("/home/FCAM/rmukherjee/Summer2022Rotation/R_analysis/paramScanCompartment")

# 45min vs control 
plot.simulated.composites.helper(input.data.identifier = "activated.45min.v.control",
				 output.plot.name = "model_45min.v.control.composites", 
				 baseline.condn.name = "control", condn.name = "45min" )

# 90min vs control
plot.simulated.composites.helper(input.data.identifier = "activated.90min.v.control",
				 output.plot.name = "model_90min.v.control.composites", 
				 baseline.condn.name = "control", condn.name = "90min" )

# 90min vs 45min
plot.simulated.composites.helper(input.data.identifier = "activated.90min.v.45min",
				 output.plot.name = "model_90min.v.45min.composites", 
				 baseline.condn.name = "45min", condn.name = "90min",
				 seq.end = 0.4, by.step=1e-5, tau = 5, dpk = 1e-5)

