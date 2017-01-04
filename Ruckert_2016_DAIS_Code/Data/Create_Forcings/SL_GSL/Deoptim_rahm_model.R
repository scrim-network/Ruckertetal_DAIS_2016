###################################################################################
#
#  -file = "DEoptim_rahm_model.R"   Code written March 2014
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This function sources the Rahmstorf(2007) model to be sourced in the DEoptim
#       R function.
#
# To use this function, simply source this file:
#   source("DEoptim_rahm_model.R")
#
# INPUTS:
#   vector of parameter values:
#       p[1]: alpha; sensitivity of sea-level to temperature changes
#       p[2]: T_0; temperature when sea-level anomaly is zero
#       p[3]: H_0; initial sea-level anomaly
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

model = function(p){ # p represents the parameters in a vector
   
  # Estimate the rate, dH, of sea-level change each year, Equation (S17) in Ruckert et al. (2016)
  dH=p[1]*(T-p[2])
  #p[1] = sensitivity of sea-level to temperature changes
  #p[2] = temperature when sea-level is zero
  
  # Set up empty vector for sea level anomalies.
  H_1 = rep(NA, to)
  H_1[1]=p[3] # sea-level in 1880
  
  # Run a forward euler to estimate sea-level over time
  for (i in from:to){
    H_1[i]=H_1[i-1]+dH[i-1]*timestep
  }
  
  # Return sea-level anomalies
  return(H_1)
}

#################################### END ##########################################