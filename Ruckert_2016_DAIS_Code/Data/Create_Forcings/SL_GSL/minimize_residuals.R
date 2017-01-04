#################################################################################
#
#  -file = "minimize_residuals.R"   Code written March 2014
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This function finds the sum of the absolute residuals estimated from the model
#       function in "DEoptim_rahm_model.R". Finding the sum is used with the DEoptim
#       R function to minimize the residuals and find the best values for the
#       parameters.
#
# To use this function, simply source this file:
#   source("minimize_residuals.R")
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
##################################################################################

min_res = function(p){
    # Return the sum of the absolute values
    # from sea-level minus simulated model values
    sum(abs(slr - model(p)))
}

################################## END ###########################################