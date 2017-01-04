########################################################################################
#
#  -file = "future_sea-level_DAIS.R"   Code written April 2015
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This program loads in temperature and global sea-level data for use in model
#       described in Rahmstorf [2007]. For further description and references, 
#       please read the paper. This program estimates the best estimates for sea-level and 
#       rates of sea-level based on parameters estimated with Differential evolution optimization (DEoptim).
#
#   -NOTE: This file contains data that is sourced into the model. Information
#       regarding this data can be found below:
#
#       -RCP8.5 is used to create temperature simulations to 2300
#       -RCP8.5 simulates "Business as usual" and is similar to the A2 scenario
#       -The RCP8.5 temperatures are from the CNRM-CM5 model
#           (Centre National de Recherches Meteorologiques)
#       -These RCP8.5 temperatures closly resemble historical temperatures
#           from the Goddard Insitute for Space Studies (GISS)
#
#       -Annual global land and ocean temperature anomalies (C)
#       -Anomalies with respect to the 1961-1990 average
#       -http://data.giss.nasa.gov/gistemp/tabledata_v3/GLB.Ts+dSST.txt
#
#       -The SL data in Shaffer (2014) has already adjusted the Church and White 2011 data to
#       -the present-day (AD 1961-1990) value of zero
#==============================================================================
# Copyright 2015 Kelsey Ruckert
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
########################################################################################
#------------------------------ Clear the global environment -------------------------
rm(list =ls()) #Clear all previous variables
library(DEoptim)
library(compiler)
enableJIT(3)
set.seed(1234)

#---------------------------------- Read in the data ---------------------------------
#Read in temperature and sea level data
data = read.csv("giss&RCP85tempscenarios_23_04_2015.csv", skip=1)

#Historical time frame and temperatures from NOAA
T = data[1:130,7] #temperature data
alltime = data[,8] #1880-2300 (1 yr incriments)

#RCP 8.5 temperatures from 2015-2300 and NOAA historical temps from 1880-2014
rcp85 = data[,10]    #4.284 C in 2100 and 9.849284058 C in 2300

# Historical global mean sea-levels from tide gauges & estimated errors
# church = read.table("CSIRO_Recons_gmsl_yr_2011.txt")
# year = church[, 1] #timeframe 1880 to 2009 in 1yre incriments
# err.obs = church[,3]/1000 #read in the error from the estimates

# read in the Church and White 2011 data adjusted to the present-day mean from Shaffer (2014)
if (!exists('SL')) {
  SL = scan("../../future_SL.txt")     #Reconstructed sea-level
  GSL = scan("../../future_GSL.txt")   #Time rate of change of sea-level
}

slr <- SL[239880:240009]
# Set the observational errors by adding and subtracting the errors to the sea-level values
err_pos=slr+err.obs
err_neg=slr-err.obs


#------------------------ Source the SLR Model & Run ----------------------
hindcast_length = 130 # there are 130 years from 1880 to 2009
projection_length = 421 # from 1880 to 2300
timestep=1  # timesteps are 1 year
from=2 # start from the second year since the first year is an uncertain parameter
to=hindcast_length #130

#Run DEoptim in R to find good initial parameters
source("Deoptim_rahm_model.R") # source the model
source("minimize_residuals.R") # find the minimum residuals
lower=c(0,-4,err_neg[1])  
upper=c(0.02,2,err_pos[1])
iter=1000  # specify number of iterations
outDEoptim <- DEoptim(min_res, lower, upper, 
                      DEoptim.control(itermax=iter,
                                      trace=FALSE))
print(outDEoptim$optim$bestmem)# find best initial parameters

#Note parameter estimates should be close to the best estimated MCMC heteroskedastic parameters found in Ruckert et al. (2015)
parms = c(outDEoptim$optim$bestmem[1], outDEoptim$optim$bestmem[2], outDEoptim$optim$bestmem[3])
source("sealevel_rahm_model.R") #sealevel_rahm_model.R is the model equation
slr.est = rahmfunction(parms, T)
to=projection_length  #number of years in the projection
proj.est = rahmfunction(parms, rcp85)

#---------------------------- Test Results with Plot ---------------------
plot(year, slr, pch=20, xlab="Year", ylab="Sea-level Rise [m]", sub="Respect to 1961-1990 time period")
lines(year, slr.est$sle, lwd=2, col="blue")
abline(h=0, v=1990, lty=2, col="gray")

plot(year, slr, pch=20, xlab="Year", ylab="Sea-level Rise [m]", sub="Respect to 1961-1990 time period", xlim=c(1880,2300),
     ylim=c(-0.15,3.5))
lines(alltime, proj.est$sle, lwd=2, col="blue")

plot(alltime, proj.est$slrate, lwd=2, typ="l", col="blue", xlab="Year", ylab="Sea-level Rate [m/yr]", 
     sub="Respect to 1961-1990 time period",xlim=c(1880,2300),ylim=c(0,0.02))
points(year, GSL[239880:240009], pch=20)

plot(year, slr.est$slrate, lwd=2, typ="l", col="blue", xlab="Year", ylab="Sea-level Rate [m/yr]", 
     sub="Respect to 1961-1990 time period",ylim=c(-0.01,0.01))
lines(year, GSL[239880:240009], lwd=2)
#----------------------------- Create Matrix & Export the Data to a CSV File for DAIS model ---------------------------#
futureSeaLevelData = matrix(c(alltime, proj.est$sle, proj.est$slrate),nrow=421,ncol=3)
colnames(futureSeaLevelData, do.NULL = FALSE)
colnames(futureSeaLevelData) = c("Year","Sea Level Rise [m] (1961-1990)","Sea Level Rate [m/yr] (1961-1990)")
write.csv(futureSeaLevelData, file = "../Output_Files/FutureRCP85_SLR_for_DAIS_new.csv")

######################################## END ###############################################