###################################################################################
#
#  -file = "sealevel_rahm_model.R"   Code written March 2014
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This function sources the Rahmstorf(2007) model to be sourced into the uncertainty methods as described
#       in Ruckert et al. (GRL 2015). For further
#       description and references, please read the paper.
#
# To use this function, simply source this file:
#   source("sealevel_rahm_model.R")
#
###################################################################################

rahmfunction = function(parm, T){ #inputs are parameters and temperature
    model.p = length(parm) #number of parameters in the model
  a = parm[1] # sensitivity of sea-level to temperature changes
  Ti = parm[2] # temperature when sea-level is zero
  initialvalue = parm[3] # initial value of sea-level in 1880
  
  rate = a*(T - Ti) # find the rate of sea-level change each year
  values = rep(NA,to)
  values[1] = initialvalue
  
  #Run a forward euler to estimate sea-level over time
  for(i in from:to){
    values[i] = values[i-1]+rate[i-1]*timestep
  }
  #return sea-level, sea-level rates, and number of parameters
  return(list(sle = values, slrate = rate, model.p = model.p))
}

################################### END ##########################################