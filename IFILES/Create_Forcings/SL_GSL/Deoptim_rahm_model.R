###################################################################################
#
#  -file = "DEoptim_rahm_model.R"   Code written March 2014
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This function sources the Rahmstorf(2007) model to be sourced in the DEoptim
#       R function. Deoptim is used in the bootstrap and MCMC codes as described
#       in Ruckert et al. (GRL 2015) to find good initial values for each of the
#       parameters.
#
# To use this function, simply source this file:
#   source("DEoptim_rahm_model.R")
#
###################################################################################

model = function(p){ # p represents the parameters in a vector
  dH=p[1]*(T-p[2]) # dH represent the rate of the sea-level rise
  #p[1] = sensitivity of sea-level to temperature changes
  #p[2] = temperature when sea-level is zero
  
  # sea-level values
  H_1 = rep(NA, to)
  H_1[1]=p[3] # sea-level in 1880
  for (i in from:to){
    H_1[i]=H_1[i-1]+dH[i-1]*timestep
  }
  return(H_1)
}

#################################### END ##########################################