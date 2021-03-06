#--------- Markov Chain Monte Carlo DAIS Simulations
#  - file ~ DAIS_convergence_plotsR_C.R
#  - Code written: March 2016
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This program assesses the convergence of the DAIS model using the
#       Heidelberger and Welch's convergence diagnostic and the Potential
#       scale reduction factor. For this we assess the output from multiple
#       seeds (1234 & 1780).
#
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
########### Plot the results for further analysis on convergence #################################
#install.packages("coda")
library(coda)

################################## CONVERGENCE ####################################
# Test for MCMC chain (1234) convergence:
load("Scratch/DAIS_calib_MCMC_C1234_relative_8e5.RData")

NI = length(DAIS_chains[,1])
burnin.length = (NI*0.04)+1
#burnin.length = (NI*0.02)+1
results = DAIS_chains[burnin.length:NI,]
NI = length(results[,1])

# Heidelberger and Welch's convergence diagnostic
#conv = results[length(burnin):NI,]
heidel.diag(results, eps=0.1, pvalue=0.05)

# Test for MCMC chain (1780) convergence:
load("Scratch/DAIS_calib_MCMC_C1780_relative__8e5.RData") # seed 1780

NI = length(DAIS_chains1780[,1])
results_1780 = DAIS_chains1780[burnin.length:NI,]
NI = length(results_1780[,1])

# Heidelberger and Welch's convergence diagnostic
#conv = results[length(burnin):NI,]
heidel.diag(results_1780, eps=0.1, pvalue=0.05)

## ============ TEST 3rd seed (Not used in the analysis for the paper) ============
# Test for MCMC chain (1) convergence:
# load("Scratch/DAIS_calib_MCMC_C1.RData") # seed 1

# NI = length(DAIS_chains1780[,1])
# results_1 = DAIS_chains1780[burnin.length:NI,]
# NI = length(results_1[,1])

# # Heidelberger and Welch's convergence diagnostic
# conv = results[length(burnin):NI,]
# heidel.diag(results_1, eps=0.1, pvalue=0.05)
##==============================================================================

# Potential Scale Reduction Factor
# MCMC is converged when the potental scale reduction factor is less than 1.1
heter = as.mcmc(results)
heter2 = as.mcmc(results_1780)
# heter3 = as.mcmc(results_1)

heterlist = mcmc.list(list(heter, heter2))
# heterlist = mcmc.list(list(heter, heter2, heter3))
gelman.diag(heterlist)

################################## TRACE PLOTS ####################################
## 1) - 8) Plot the Trace and Histograms for each parameter in the chain to insure convergence:
jpeg(file="Converge/converge_gammainstpaleo.jpeg", family="Helvetica", width=1200, height=700, units="px", pointsize=20)
par(mfrow=c(2,1), mar=c(3, 7, 1, 7), mgp=c(1.5,.5,0))
plot(results[,1], type="l",
     ylab="Gamma [sensitivity of ice flow to sea level [dimensionless]]",
     xlab="Number of Runs", main="")
par(mar=c(3,7,1,7))
hist(results[,1], freq=FALSE, col="gray",
     xlab="Gamma [sensitivity of ice flow to sea level [dimensionless]]",
     ylab="Density [Dimensionless]", main="")
dev.off()

# ------------------------------------------------------------------------
jpeg(file="Converge/converge_alphainstpaleo.jpeg", family="Helvetica", width=1200, height=700, units="px", pointsize=20)
par(mfrow=c(2,1), mar=c(3, 7, 1, 7), mgp=c(1.5,.5,0))
plot(results[,2], type="l",
     ylab="Alpha [sensitivity of ice flow to ocean subsurface temperature [dimensionless]]",
     xlab="Number of Runs", main="")
par(mar=c(3,7,1,7))
hist(results[,2], freq=FALSE, col="gray",
     xlab="Alpha [sensitivity of ice flow to ocean subsurface temperature [dimensionless]]",
     ylab="Density [Dimensionless]", main="")
dev.off()

# ------------------------------------------------------------------------
jpeg(file="Converge/converge_muinstpaleo.jpeg", family="Helvetica", width=1200, height=700, units="px", pointsize=20)
par(mfrow=c(2,1), mar=c(3, 7, 1, 7), mgp=c(1.5,.5,0))
plot(results[,3], type="l",
     ylab="Mu [Profile parameter related to ice stress [m^(1/2)]]",
     xlab="Number of Runs", main="")
par(mar=c(3,7,1,7))
hist(results[,3], freq=FALSE, col="gray",
     xlab="Mu [Profile parameter related to ice stress [m^(1/2)]]",
     ylab="Density [Dimensionless]", main="")
dev.off()

# ------------------------------------------------------------------------
jpeg(file="Converge/converge_nuinstpaleo.jpeg", family="Helvetica", width=1200, height=700, units="px", pointsize=20)
par(mfrow=c(2,1), mar=c(3, 7, 1, 7), mgp=c(1.5,.5,0))
plot(results[,4], type="l",
     ylab="Nu [Relates balance gradient to precipitation [m^(-1/2) yr^(-1/2)]]",
     xlab="Number of Runs", main="")
par(mar=c(3,7,1,7))
hist(results[,4], freq=FALSE, col="gray",
     xlab="Nu [Relates balance gradient to precipitation [m^(-1/2) yr^(-1/2)]]",
     ylab="Density [Dimensionless]", main="")
dev.off()

# ------------------------------------------------------------------------
jpeg(file="Converge/converge_Poinstpaleo.jpeg", family="Helvetica", width=1200, height=700, units="px", pointsize=20)
par(mfrow=c(2,1), mar=c(3, 7, 1, 7), mgp=c(1.5,.5,0))
plot(results[,5], type="l",
     ylab="Po [Precipitation at 0C [m of ice/yr]]",
     xlab="Number of Runs", main="")
par(mar=c(3,7,1,7))
hist(results[,5], freq=FALSE, col="gray",
     xlab="Po [Precipitation at 0C [m of ice/yr]]",
     ylab="Density [Dimensionless]", main="")
dev.off()

# ------------------------------------------------------------------------
jpeg(file="Converge/converge_kappainstpaleo.jpeg", family="Helvetica", width=1200, height=700, units="px", pointsize=20)
par(mfrow=c(2,1), mar=c(3, 7, 1, 7), mgp=c(1.5,.5,0))
plot(results[,6], type="l",
     ylab="Kappa [Relates precipitation to temperature [K^-1]]",
     xlab="Number of Runs", main="")
par(mar=c(3,7,1,7))
hist(results[,6], freq=FALSE, col="gray",
     xlab="Kappa [Relates precipitation to temperature [K^-1]]",
     ylab="Density [Dimensionless]", main="")
dev.off()

# ------------------------------------------------------------------------
jpeg(file="Converge/converge_foinstpaleo.jpeg", family="Helvetica", width=1200, height=700, units="px", pointsize=20)
par(mfrow=c(2,1), mar=c(3, 7, 1, 7), mgp=c(1.5,.5,0))
plot(results[,7], type="l",
ylab="Fo [Constant of proportionality for ice speed [m/yr]]",
xlab="Number of Runs", main="")
par(mar=c(3,7,1,7))
hist(results[,7], freq=FALSE, col="gray",
xlab="Fo [Constant of proportionality for ice speed [m/yr]]",
ylab="Density [Dimensionless]", main="")
dev.off()

# ------------------------------------------------------------------------
jpeg(file="Converge/converge_hoinstpaleo.jpeg", family="Helvetica", width=1200, height=700, units="px", pointsize=20)
par(mfrow=c(2,1), mar=c(3, 7, 1, 7), mgp=c(1.5,.5,0))
plot(results[,8], type="l",
ylab="Ho [Initial value for runoff line calculation [m]]",
xlab="Number of Runs", main="")
par(mar=c(3,7,1,7))
hist(results[,8], freq=FALSE, col="gray",
xlab="Ho [Initial value for runoff line calculation [m]]",
ylab="Density [Dimensionless]", main="")
dev.off()

# ------------------------------------------------------------------------
jpeg(file="Converge/converge_coinstpaleo.jpeg", family="Helvetica", width=1200, height=700, units="px", pointsize=20)
par(mfrow=c(2,1), mar=c(3, 7, 1, 7), mgp=c(1.5,.5,0))
plot(results[,9], type="l",
ylab="Co [Second value for runoff line calculation [m]]",
xlab="Number of Runs", main="")
par(mar=c(3,7,1,7))
hist(results[,9], freq=FALSE, col="gray",
xlab="Co [Second value for runoff line calculation [m]]",
ylab="Density [Dimensionless]", main="")
dev.off()

# ------------------------------------------------------------------------
jpeg(file="Converge/converge_boinstpaleo.jpeg", family="Helvetica", width=1200, height=700, units="px", pointsize=20)
par(mfrow=c(2,1), mar=c(3, 7, 1, 7), mgp=c(1.5,.5,0))
plot(results[,10], type="l",
     ylab="Bo [Height of bed at the center of the continent [m]]",
     xlab="Number of Runs", main="")
par(mar=c(3,7,1,7))
hist(results[,10], freq=FALSE, col="gray",
     xlab="Bo [Height of bed at the center of the continent [m]]",
     ylab="Density [Dimensionless]", main="")
dev.off()

# ------------------------------------------------------------------------
jpeg(file="Converge/converge_sinstpaleo.jpeg", family="Helvetica", width=1200, height=700, units="px", pointsize=20)
par(mfrow=c(2,1), mar=c(3, 7, 1, 7), mgp=c(1.5,.5,0))
plot(results[,11], type="l",
     ylab="s [Slope of the bed [%]]",
     xlab="Number of Runs", main="")
par(mar=c(3,7,1,7))
hist(results[,11], freq=FALSE, col="gray",
     xlab="s [Slope of the bed [%]]",
     ylab="Density [Dimensionless]", main="")
dev.off()

# ------------------------------------------------------------------------
jpeg(file="Converge/converge_paleosigmainstpaleo.jpeg", family="Helvetica", width=1200, height=700, units="px", pointsize=20)
par(mfrow=c(2,1), mar=c(3, 7, 1, 7), mgp=c(1.5,.5,0))
plot(results[,12], type="l",
     ylab="Paleo Sigma^2 [Dimensionless]",
     xlab="Number of Runs", main="")
par(mar=c(3,7,1,7))
hist(results[,12], freq=FALSE, col="gray",
     xlab="Paleo Sigma^2 [Dimensionless]",
     ylab="Density [Dimensionless]", main="")
dev.off()

# ------------------------------------------------------------------------
jpeg(file="Converge/converge_instsigmainstpaleo.jpeg", family="Helvetica", width=1200, height=700, units="px", pointsize=20)
par(mfrow=c(2,1), mar=c(3, 7, 1, 7), mgp=c(1.5,.5,0))
plot(results[,13], type="l",
ylab="Inst. Sigma^2 [Dimensionless]",
xlab="Number of Runs", main="")
par(mar=c(3,7,1,7))
hist(results[,13], freq=FALSE, col="gray",
xlab="Inst. Sigma^2 [Dimensionless]",
ylab="Density [Dimensionless]", main="")
dev.off()

# ------------------------------------------------------------------------
## To check if subset is sufficient:
#They should be roughly similiar
#subset_N = NI/2500 #calculate the every nth number to get a subset of ~2500
#R_subset = round(subset_N,0)
#sschain = results[seq(1, length(results[,1]), R_subset),]
#subset_length = length(sschain[,1])

subset_length = 3500
sschain = results[sample(nrow(results), size=subset_length, replace=FALSE), ]

pdf(file="Converge/subset_checkinstpaleo.pdf")
par(mfrow=c(4,4))
plot(density(results[,1]), main="gamma", xlab="")
lines(density(sschain[,1]), col="red")
plot(density(results[,2]), main="alpha", xlab="")
lines(density(sschain[,2]), col="red")
plot(density(results[,3]), main="mu", xlab="")
lines(density(sschain[,3]), col="red")
plot(density(results[,4]), main="nu", xlab="")
lines(density(sschain[,4]), col="red")
plot(density(results[,5]), main="Po", xlab="")
lines(density(sschain[,5]), col="red")
plot(density(results[,6]), main="kappa", xlab="")
lines(density(sschain[,6]), col="red")
plot(density(results[,7]), main="fo", xlab="")
lines(density(sschain[,7]), col="red")
plot(density(results[,8]), main="ho", xlab="")
lines(density(sschain[,8]), col="red")
plot(density(results[,9]), main="c", xlab="")
lines(density(sschain[,9]), col="red")
plot(density(results[,10]), main="bo", xlab="")
lines(density(sschain[,10]), col="red")
plot(density(results[,11]), main="slope", xlab="")
lines(density(sschain[,11]), col="red")
plot(density(results[,12]), main="var.paleo", xlab="")
lines(density(sschain[,12]), col="red")
plot(density(results[,13]), main="var.inst", xlab="")
lines(density(sschain[,13]), col="red")
dev.off()

#------------------------------ END ---------------------------------------------------------------------


