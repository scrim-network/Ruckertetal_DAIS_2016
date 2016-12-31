#################################################################################
#
#  - file = "DAIScali_hetero_model_iid_mcmc_R_C1780.R"
#  - Code written: August 2015, updated March 2016
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This program runs a Markov Chain Monte Carlo analysis of the DAIS model
#       assuming heteroskedastic errors and IID residuals as described in Ruckert et al. (2016).
#       For further description and references, please read the paper.
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOT IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
###################################################################################
rm(list =ls()) # Clear global environment
library(compiler)
enableJIT(3)
enableJIT(3)
library(mcmc)

#Set the seed.
set.seed(1780)

# Read in hindcast/ forcing data, read in standard values, and read in AIS dates both specific and ranges so there is no use of magic numbers.
source("Data/DAIS_data_C.R")

############################## List of physcial model parameters ##############################
# We will set the initial parameters to specifications from Shaffer [2014]
#Define the parameters:
# [1] gamma = 2 				#sensitivity of ice flow to sea level
# [2] alpha = 0.35 			#sensitivity of ice flow to ocean subsurface temperature
# [3] mu = 8.7    			#Profile parameter related to ice stress [m^(1/2)]
# [4] nu = 0.012   			#Relates balance gradient to precipitation [m^(-1/2) yr^(-1/2)]
# [5] P0 = 0.35    			#Precipitation at 0C [m of ice/yr]
# [6] kappa = 0.04 			#Relates precipitation to temperature [K^-1]
# [7] f0 = 1.2          #Constant of proportionality for ice speed [m/yr]
# [8] h0 = 1471         #Initial value for runoff line calculation [m]
# [9] c = 95            #Second value for runoff line calculation [m]
# [10] b0 = 775         #Height of bed at the center of the continent [m]
# [11] slope = 0.0006   #Slope of the bed

# Create matrices for projections and hindcasts.
project.forcings = matrix(c(Ta, Toc, GSL, SL), ncol=4, nrow=240300)
hindcast.forcings = matrix(c(Ta[1:240010], Toc[1:240010], GSL[1:240010], SL[1:240010]), ncol=4, nrow=240010)

# SEt initial parameters to the best case (Case #4) from Shaffer (2014)
IP = c(2, 0.35, 8.7, 0.012, 0.35, 0.04, 1.2, 1471, 95, 775, 0.0006)

#Source the function with the standards and the initial parameters (IP) to
#get the best estimated AIS volume loss 

# Load in the physical model. Calls C model by default and set vector of standard values:
source("models.R") 
standards = c(Tf, rho_w, rho_i, rho_m, Toc_0, Rad0)

# Estimate the AIS volume loss hindcast for Case #4 with respect to the present day in sea level equivalence (SLE):
AIS_melt = iceflux(IP, hindcast.forcings, standards)

# Estimate the AIS volume loss projection for Case #4 with respect to the present day in sea level equivalence (SLE):
Project_melt = iceflux(IP, project.forcings, standards)

############################## Setup observational constraint info ##############################
# For this model the residuals are based off of the windowing approach.
# These windows are presented in Shaffer (2014) and calculated from figure 5 in Shepherd et al. (2012).

# Accummulate the sea-level equivalent in meters from 1992 to the year 2002
# using the 1992 to 2011 trend from Shepherd et al. 2012; -71 +/- 53 Gt per yr.
# Conversion: 360 Gt = 1 mm SLE
estimate.SLE.rate = abs(-71/360)/1000
time.years = 2002 - 1992
mid.cum.SLE_2002 = estimate.SLE.rate * time.years

# Accumulate the 1 sigma error for the year 2002 and estimate the 2-sigma error:
estimate.SLE.error = abs(-53/360)/1000 #1- sigma error
estimate.SLE.error = sqrt(time.years)*abs(-53/360)/1000 # 1-sigma error
# (*sqrt(10) because 10 years of potentially accumulated error:
#  total error^2 = year 1 error^2 + year 2 error^2 + ... year 10 error^2
#                = 10*year X error^2)
SE2_2002 = estimate.SLE.error*2 # 2-sigma error

# Add and subtract the 2 standard error to the mean value
positive_2SE = mid.cum.SLE_2002 + SE2_2002 
negative_2SE = mid.cum.SLE_2002 - SE2_2002

# Create observational constraint windows.
upper.wind = c(6.0, -6.9, -1.25, positive_2SE )
lower.wind = c(1.8, -15.8, -4.0, negative_2SE)
windows = matrix(c(lower.wind, upper.wind), nrow = 4, ncol=2)

# Determine observational error from windows: half-width of window = uncertainty; assume all windows are 2*stdErr (last one actually is)
obs.errs = (windows[,2]-windows[,1])*.5
 
# Create a vector with each observation year.
#            120 kyr, 20 Kyr,  6 kyr,  2002
obs.years = c(120000, 220000, 234000, 240002)

############################## CALCULATE RESIDUALS and INITIAL SIGMA VALUE ##############################
#resid <- rep(NA,length(obs.years))

# Estimate the residuals: modification from equation (1)
#for(i in 1:length(obs.years)){
#    resid[i] <- (median(windows[i,]) - (AIS_melt[obs.years[i]] - mean(AIS_melt[SL.1961_1990])))
#}

resid.1 <- (median(windows[1,]) - (AIS_melt[obs.years[1]] - mean(AIS_melt[SL.1961_1990])))
resid.2 <- (median(windows[2,]) - (AIS_melt[obs.years[2]] - mean(AIS_melt[SL.1961_1990])))
resid.3 <- (median(windows[3,]) - (AIS_melt[obs.years[3]] - mean(AIS_melt[SL.1961_1990])))
resid.4 <- (median(windows[4,]) - (AIS_melt[obs.years[4]] - AIS_melt[239992]))

resid <- c(resid.1, resid.2, resid.3)

# Calculate the variance, sigma^2
paleo_variance = sd(resid)^2
inst_variance = resid.4^2

############################## RUN MCMC CALIBRATION #######################################
# Set up priors.
bound.lower = IP - (IP*0.5)    ; bound.upper = IP + (IP*0.5)
print(bound.lower)
print(bound.upper)

# var.y has inverse gamma prior, so there is a lower bound at 0 but no upper bound
parnames    = c('gamma','alpha','mu'  ,'nu'  ,'P0' ,'kappa','f0' ,'h0'  ,'c'  , 'b0','slope' ,'var.paleo', 'var.inst')
bound.upper = c( 4.25 ,  1     , 13.05, 0.018,0.525,  0.06 , 1.8 ,2206.5, 142.5, 825 , 0.00075,   Inf,       0.0004)
bound.lower = c( 0.5  ,  0     , 4.35 , 0.006,0.175,  0.02 , 0.6 , 735.5,  47.5, 725 , 0.00045 ,    0,       0)

# Specify the number of model parameters.
# Variance is a statistical parameter and is not counted in the number.
model.p=11

# Load the likelihood model assuming heteroskedastic observation errors and non-correlated residuals.
source("Scripts/DAISobs_likelihood_iid.R")

# Optimize the likelihood function to estimate initial starting values.
p = c(IP, paleo_variance, inst_variance) # Shaffer [2014] Case #4 parameters
p0 = c(2.1, 0.29, 8, 0.015, 0.4, 0.04, 1.0, 1450, 90, 770, 0.0005, 0.6, (5.5e-5)^2) # Random guesses
p0 = optim(p0, function(p) - log.post(p))$par
print(round(p0,4))

# Set the step size and number of iterations.
#step = c(0.1, 0.015, 0.2, 0.035, 0.1, 0.01, 0.1, 50, 10, 25, 0.0005, 0.1)/5
#step = c(0.05, 0.01, 0.15, 0.035, 0.1, 0.01, 0.1, 50, 10, 30, 0.0005, 0.1)/5
#step = c(0.1, 0.015, 0.2, 0.025, 0.1, 0.01, 0.1, 50, 10, 20, 0.0005, 0.15)/5
NI = 8E5

# Run MCMC calibration.
#DAIS_mcmc_output1780 = metrop(log.post, p0, nbatch=NI, scale=step)
#DAIS_chains1780 = DAIS_mcmc_output1780$batch

# Print the acceptance rate as a percent. Should be ~ 25%
#acceptrate = DAIS_mcmc_output1780$accept * 100
#cat("Accept rate =", acceptrate, "%\n")

############################## RUN ADAPTIVE MCMC #######################################
library(adaptMCMC)

# Set optimal acceptance rate as # parameters->infinity (Gelman et al, 1996; Roberts et al, 1997)
accept.mcmc = 0.234

# Set number of iterations, burnin, and rate of adaptation (between 0.5 and 1, lower is faster adaptation)
gamma.mcmc = 0.5											#

# # Specify when to stop adapting (niter*1 => don't stop) and step size
# # stopadapt.mcmc = round(niter.mcmc*1.0)
step.mcmc = (bound.upper-bound.lower)*.05
step.mcmc = c(step.mcmc[1:11], 0.05, 1e-8)

## Actually run the calibration.
DAIS_mcmc_output1780 = MCMC(log.post, NI, p0, scale=step.mcmc, adapt=TRUE, acc.rate=accept.mcmc,gamma=gamma.mcmc, list=TRUE,
                            n.start=round(0.01*NI))
DAIS_chains1780 = DAIS_mcmc_output1780$samples

## Save workspace image - you do not want to re-simulate all those!
save.image(file = "DAIS_calib_MCMC_C1780_relative__8e5.RData")

########################################## END ########################################


