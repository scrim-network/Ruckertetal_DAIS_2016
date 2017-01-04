###################################################################################
#
#  -file = "sealevel_rahm_model.R"   Code written March 2014
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This function sources the Rahmstorf(2007) model.
#
# To use this function, simply source this file:
#   source("sealevel_rahm_model.R")
#
# INPUTS:
#   vector of parameter values:
#       parameters[1]: alpha; sensitivity of sea-level to temperature changes
#       parameters[2]: T_0; temperature when sea-level anomaly is zero
#       parameters[3]: H_0; initial sea-level anomaly
#
#   vector of temperatures
#==============================================================================
# Copyright 2014 Kelsey Ruckert
# This file is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This file is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this file.  If not, see <http://www.gnu.org/licenses/>.
#
###################################################################################

rahmfunction = function(parm, T){ #inputs are parameters and temperature
    
  # Determine number of model parameters
  model.p = length(parm)
  
  # Extract parameter values
  a = parm[1]            # sensitivity of sea-level to temperature changes
  Ti = parm[2]           # temperature when the sea-level anomaly is zero
  initialvalue = parm[3] # initial value of sea-level in 1880
  
  # Estimate the rate of sea-level change each year, Equation (S17) in Ruckert et al. (2016)
  rate = a*(T - Ti)
  
  # Set up empty vector for sea level anomalies.
  values = rep(NA,to)
  values[1] = initialvalue
  
  # Run a forward euler to estimate sea-level over time
  for(i in from:to){
    values[i] = values[i-1]+rate[i-1]*timestep
  }
  
  # Return sea-level anomalies, sea-level rates, and number of parameters
  return(list(sle = values, slrate = rate, model.p = model.p))
}

################################### END ##########################################