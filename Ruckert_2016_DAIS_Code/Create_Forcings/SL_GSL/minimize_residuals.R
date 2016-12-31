#################################################################################
#
#  -file = "minimize_residuals.R"   Code written March 2014
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This function find the sum of the absolute residuals estimated from the model
#       function in "DEoptim_rahm_model.R". Finding the sum is used with the DEoptim
#       R function to minimize the residuals and find the best values for the
#       parameters. The best parameter values are used as initial values for the
#       bootstrap and MCMC codes as described in Ruckert et al. (GRL 2015).
#
# To use this function, simply source this file:
#   source("minimize_residuals.R")
#
##################################################################################

min_res = function(p){
    sum(abs(slr - model(p))) # returns the sum of the absolute values
                             # from sea-level minus simulated model values
}

################################## END ###########################################