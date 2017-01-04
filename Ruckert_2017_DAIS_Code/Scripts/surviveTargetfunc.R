#############################################################
## -file = surviveTargetfunc.R
## -Input: The target range interested in [window] which includes:
##         The lower bound of the range [lbw] &
##         The upper bound of the range [upw]
##
##         The vector of data [sv]
##
## -Output: Which runs in a vector fall within a certain range [Target]
#############################################################
## Copyright 2016 Kelsey Ruckert
## This file is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This file is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this file.  If not, see <http://www.gnu.org/licenses/>.
##
## -Author: Kelsey Ruckert (klr324@psu.edu)
## -June 10 2015
###############################################################

surviveTarget = function(window, surtarget){
  lbw = window[1]
  upw = window[2]
  sv = surtarget
  
  which((sv >= lbw & sv <= upw))
}
