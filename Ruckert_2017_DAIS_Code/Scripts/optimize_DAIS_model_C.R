###################################
## file: optimize_DAIS_model_C.R
###################################
## Author and copyright: Kelsey Ruckert
## Pennsylvania State University
# klr324@psu.edu
## Date: April 2015; updated July 2016
###################################
## Function of the Shaffer DAIS model 2014
## This function/DAIS model returns the root mean square error of the hindcast
## compared to observational constraints. The point of this function is to be
## used for optimization.
##==============================================================================
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
##==============================================================================
###################################

# DAIS Model
optimize_dais = function(parameters, forcings, standards, windows, obs.years){
  
  Gamma = parameters[1]         # sensitivity of ice flow to sea level
  alpha = parameters[2]         # sensitivity of ice flow to ocean subsurface temperature
  mu = parameters[3]            # Profile parameter related to ice stress [m^(1/2)]
  nu = parameters[4]            # Relates balance gradient to precipitation [m^(-1/2) yr^(-1/2)]
  P0 = parameters[5]            # Precipitation at 0C [m of ice/yr]
  kappa = parameters[6]         # Relates precipitation to temperature [K^-1]
  f0 = parameters[7]            # Constant of proportionality for ice speed [m/yr]
  h0 = parameters[8]            # Initial value for runoff line calculation [m]
  c = parameters[9]             # Second value for runoff line calculation [m]
  b0 =  parameters[10]          # Height of bed at the center of the continent [m]
  slope =  parameters[11]       # Slope of the bed

  
  Ta = forcings[,1]
  Toc = forcings[,2]
  GSL = forcings[,3]
  SL = forcings[,4]
  
  Tf    = standards[1]          # Freezing temperature of sea water
  rho_w = standards[2]          # Density of sea water [g/cm^3]
  rho_i = standards[3]          # Density of ice water [g/cm^3]
  rho_m = standards[4]          # Density of rock [g/cm^3]
  Toc_0 = standards[5]          # Present day high latitude ocean subsurface temperature [K]
  Rad0  = standards[6]          # Steady state AIS radius for present day Ta and SL [m]
  Volo  = standards[7]          
  
  # Initialize intermediate parameters.
  del  = rho_w/rho_i            # Ratio sea water and ice density [-]
  eps1 = rho_i/(rho_m - rho_i)  # Ratio ice density and density difference between rock and ice [-]
  eps2 = rho_w/(rho_m - rho_i)  # Ratio sea water density and density difference between rock and ice [-]
  R  = Rad0                     # gets updated at end of loop

  Rad = rep(NA, length(Ta))
  Vol = rep(NA, length(Ta))
  SLE = rep(NA, length(Ta))
  
  # Run model
  for(i in 1:length(Ta)){
    # Function that is used in calculating Ice speed and Ice Flux (modified from equation 11)
    f = f0*((1 - alpha) + alpha*((Toc[i] - Tf)/(Toc_0 - Tf))^2)/((slope*Rad0 - b0)^(Gamma - 1))
    hr= h0 + c * Ta[i] # equation 5
    rc = (b0 - SL[i])/slope # application of equation 1 (paragraph after eq3)
    P = P0 * exp(kappa * Ta[i]) # equation 6
    beta = nu * sqrt(P) # equation 7 (corrected with respect to text)
    rR = R - (abs(hr - b0 + slope * R)*(hr - b0 + slope*R))/mu # Distance from the continent center to where the runoff line intersects the ice sheet surface.
    
    if (R<=rR && R<=rc) {
      # Total mass accumulation on ice sheet (equation 8)
      Btot = pi * P * R * R
      # In case there is no marine ice sheet / grounding line
      F = 0     # no ice flux
      ISO = 0   # (third term equation 14) NAME?
      fac = pi * (1 + eps1) * (4/3 * mu^0.5 * R^1.5 - slope*R*R) # ratio dV/dR
      
    } else if (R>rR && R<=rc) {
      # Total mass accumulation on ice sheet minus runoff
      Btot = pi * P * R * R - pi * beta*(hr - b0 + slope * R) *
      (R*R - rR*rR) - (4 * pi * beta * mu^0.5)/5 * (R - rR)^2.5 +
      (4 * pi * beta * mu^0.5)/3 * R * (R - rR)^1.5
      # In case there is no marine ice sheet / grounding line
      F = 0     # no ice flux
      ISO = 0   # (third term equation 14) NAME?
      fac = pi * (1 + eps1) * (4/3 * mu^0.5 * R^1.5 - slope*R*R) # ratio dV/dR
      
    } else if (R<=rR && R>=rc) {
      # Total mass accumulation with marine ice sheet / grounding line
      Btot = pi * P * R * R
      Hw = slope * R - b0 + SL[i]  # (equation 10)
      F = 2 * pi * R * f * del * Hw^(Gamma+1)   # Ice flux (equation 9)
      ISO = 2 * pi * eps2 * (slope * rc * rc - b0/slope * rc) * GSL[i] # third term equation 14 !! NAME?
      fac = pi * (1 + eps1) * (4/3 * mu^0.5 * R^1.5 - slope*R*R) - 2*pi*eps2*(slope*R*R - b0*R)
      
    } else {
      # Total mass accumulation minus runoff with marine ice sheet / grounding line
      Btot = pi * P * R * R - pi * beta * (hr - b0 + slope*R) *
      (R*R - rR*rR) - (4 * pi * beta * mu^0.5)/5 *
      ((R-rR)^2.5) + (4 * pi * beta * mu^0.5)/3 * (R * (R-rR)^1.5)
      Hw = slope * R - b0 + SL[i]  # (equation 10)
      F = 2 * pi * R * f * del * Hw^(Gamma+1)   # Ice flux (equation 9)
      ISO = 2 * pi * eps2 * (slope * rc * rc - b0/slope * rc) * GSL[i] # third term equation 14 !! NAME?
      fac = pi * (1 + eps1) * (4/3 * mu^0.5 * R^1.5 - slope*R*R) - 2*pi*eps2*(slope*R*R - b0*R)
    }
    
    dR = (Btot-F+ISO)/fac
    R = R+dR # Estimate new radius
    V = 8/15 * pi * mu^0.5 * R^2.5 - 1/3 * pi * slope * R^3
    
    # Calculate sea volume
    Vsea = pi * (2/3 *slope * (R^3 - rc^3) - b0 * (R^2 - rc^2))
    
    # Calulate the volume change over time
    if(R<=rc){
      Volt = (1 + eps1) * V
    }
    else{
      Volt = (1 + eps1) * V - eps2 * Vsea
    }
    
    # Ice sheet volume (equation 13)
    Rad[i] = R
    Vol[i] = Volt
    SLE[i] = 57 * (1 - Vol[i]/Volo) # Takes steady state present day volume to correspond to 57m SLE
  }
  
  # Estimate the residuals
  #  wresid <- rep(NA,length(obs.years))
  #for(i in 1:length(obs.years)){
  #  wresid[i] <- (median(windows[i,]) - (SLE[obs.years[i]] - mean(SLE[SL.1961_1990])))
  #}
  
  wresid.1 <- (median(windows[1,]) - (SLE[obs.years[1]] - mean(SLE[SL.1961_1990])))
  wresid.2 <- (median(windows[2,]) - (SLE[obs.years[2]] - mean(SLE[SL.1961_1990])))
  wresid.3 <- (median(windows[3,]) - (SLE[obs.years[3]] - mean(SLE[SL.1961_1990])))
  wresid.4 <- (median(windows[4,]) - (SLE[obs.years[4]] - SLE[239992]))
  
  wresid <- c(wresid.1, wresid.2, wresid.3, wresid.4)
  
  # Estimate and return the root means square error
  wrmse = sqrt(mean(wresid^2))
  return(wrmse)
}
