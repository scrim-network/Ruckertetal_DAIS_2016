#############################################################
## -file = DAIS_data_C.R
## -Antarctic Ice Sheet (AIS) model
#############################################################
## -Model can be found in Shaffer_GMDD_2014
## -Matlab codes and forcings can be found at:
## -www.dcess.dk under "DCESS models"
## -Sea-level values are reconstructed from:
## -Waelbroeck et al 2002(21-240kyr BP), Clark et al 2012(7000-21000BP),
## -Lambeck et al 2010(6000BP-1869AD), and Church & White 2011(1870-2010AD)
## -Future Sea level values and and rates are calculated using the Rahmstorf (2007) model, DEoptim, and RCP8.5
## -Future air and ocean temperatures were calculated by Robert Nicholas and Varada Vaidya using the CNRM-CMIP5 model
## -Author: Kelsey Ruckert (klr324@psu.edu)
##############################################################
## -June 17, 2014 #Updates June 10 2015
##==============================================================================
## Copyright 2015 Kelsey Ruckert
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
##==============================================================================
##############################################################
# read in forcing data
if (!exists('SL')) {
  date = seq(-239999,300,1) #240 Kyr BP to 2100AD at one year intervals of the forcings
  GSL = scan("Data/future_GSL.txt")   #Time rate of change of sea-level
  Ta = scan("Data/future_TA.txt")     #Antarctic temperature reduced to sea-level
  Toc = scan("Data/future_TO.txt")     #High latitude subsurface ocean temperature
  SL = scan("Data/future_SL.txt")     #Reconstructed sea-level
}

#Set model parameters at their standard values
Tf = -1.8             #Freezing temperature of sea water
rho_w = 1030              #Density of sea water [g/cm^3]
rho_i = 917            #Density of ice water [g/cm^3]
rho_m = 4000             #Density of rock [g/cm^3]
b0 = 775                #Height of bed at the center of the continent [m]
slope = 0.0006              #Slope of the bed
f0 = 1.2                #Constant of proportionality for ice speed [m/yr]
h0 = 1471               #Initial value for runoff line calculation [m]
c = 95                  #Second value for runoff line calculation [m]
Rad0 = 1.8636e6          #Steady state AIS radius for present day Ta and SL [m]
Volo = 2.4789e16        #Steady state AIS volume for present day Ta and SL [m^3]
Toc_0 = 0.72              #Present day high latitude ocean subsurface temperature [K]

del = rho_w/rho_i
eps1 = rho_i/(rho_m - rho_i)
eps2 = rho_w/(rho_m - rho_i)

#Initial condition for integration
R = Rad0

#Setup AIS melt ranges and specific dates so there is no use of magic numbers
#in the code
last.interglacial = c(date[110000], date[120000], date[130000])
last.glacialmax = c(date[219000], date[220000], date[221000])
holocene = c(date[233800], date[234000], date[234200])
SL.1961_1990 = 239961:239990
kyrbp_25 = 220000
kyrbp_6 = 234000
AD_1880 = 239880
present = 240002
end = 240000
enddate = 240010
