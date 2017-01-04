###################################
## file: model.R  05 Feb. 2016
##
## Author: Robert Fuller   (rwf136@psu.edu)
## Edits by Kelsey Ruckert (klr324@psu.edu)
## Pennsylvania State University
##
## Note: This code calls the DAIS model written in different languages:
## R, FORTRAN, and C. The default calls the C model, which is the version
## used in Ruckert et al. (2017).
##
## Copyright 2016 Robert Fuller, Kelsey Ruckert
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
useCmodel <- T
useFmodel <- F

# Call the R DAIS function
if (!useCmodel && !useFmodel) {
    source("Scripts/DAIS_IceFlux_model.R")
}

# Call the FORTRAN DAIS function
if (useFmodel) {
    source("daisF.R")

    iceflux <- function(iceflux, forcings, standards)
    {
        Volume_F <- daisF(
          tstep = 1,
          b0    = iceflux[10],
          slope = iceflux[11],
          mu    = iceflux[3],
          h0    = iceflux[8],
          c     = iceflux[9],
          P0    = iceflux[5],
          kappa = iceflux[6],
          nu    = iceflux[4],
          f0    = iceflux[7],
          gamma = iceflux[1],
          alpha = iceflux[2],
          Tf    = -1.8,             #Freezing temperature of sea water
          rho_w = 1030,             #Density of sea water [g/cm^3]
          rho_i = 917,              #Density of ice water [g/cm^3]
          rho_m = 4000,             #Density of rock [g/cm^3]
          Toc_0 = 0.72,             #Present day high latitude ocean subsurface temperature [K]
          Rad0  = 1.8636e6,         #Steady state AIS radius for present day Ta and SL [m]
          Ta     = forcings[, 1], 
          SL     = forcings[, 4],
          Toc    = forcings[, 2],
          dSL    = forcings[, 3])

        return (Volume_F)
    }
}

# Call the C DAIS function
if (useCmodel) {
    source("roblib.R")
    dynReload("dais", srcname=c("dais.c", "r.c"), extrasrc="r.h")

    iceflux <- function(iceflux, forcings, standards)
    {
        mp <- c(
          b0    = iceflux[10],
          slope = iceflux[11],
          mu    = iceflux[3],
          h0    = iceflux[8],
          c     = iceflux[9],
          P0    = iceflux[5],
          kappa = iceflux[6],
          nu    = iceflux[4],
          f0    = iceflux[7],
          gamma = iceflux[1],
          alpha = iceflux[2],
          Tf    = standards[1],             #Freezing temperature of sea water
          rho_w = standards[2],             #Density of sea water [g/cm^3]
          rho_i = standards[3],             #Density of ice water [g/cm^3]
          rho_m = standards[4],             #Density of rock [g/cm^3]
          Toc_0 = standards[5],             #Present day high latitude ocean subsurface temperature [K]
          Rad0  = standards[6]          #Steady state AIS radius for present day Ta and SL [m]
        )

        np     <- nrow(forcings)
        Rad    <- numeric(length=np)               # Radius of ice sheet
        Vais   <- numeric(length=np)               # Ice volume
        SLE    <- numeric(length=np)               # Sea-level equivalent [m]

        .Call("daisOdeC", list(mp=mp, frc=forcings, out=list(SLE, Vais, Rad)), PACKAGE = "dais")

        return(SLE)
    }
}
