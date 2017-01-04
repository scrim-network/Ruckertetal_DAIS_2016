#####################################################################################
# R Function: DAISobs_likelihood_iid.R
# -Antarctic Ice Sheet (AIS) model
#
# compute (log) likelihood for observations
# The observations are indepent and identically distributed (IID)
#
# -Author: Yawen Guan (yig5031@psu.edu) and Kelsey Ruckert (klr324@psu.edu)
# -June 10 2015
##==============================================================================
## Copyright 2015 Yawen Guan, Kelsey Ruckert
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
#####################################################################################

log.lik = function(p) # model.p is the dimension of model parameters 
{ 
  par=p[1:model.p]
  
  var.paleo = p[model.p+1]
  var.inst = p[model.p+2]
  
  model.out = iceflux(par, hindcast.forcings, standards)
  y.mod = model.out
  
  llik.y  = 0
  
  # get the residuals
  # the paleo observational constraints are relative to 1961-1990 average
  r1 = median(windows[1,]) - (y.mod[120000] - mean(y.mod[SL.1961_1990]))
  r2 = median(windows[2,]) - (y.mod[220000] - mean(y.mod[SL.1961_1990]))
  r3 = median(windows[3,]) - (y.mod[234000] - mean(y.mod[SL.1961_1990]))
  #r4 = median(windows[4,]) - (y.mod[240002] - mean(y.mod[SL.1961_1990]))
  
  # the 2002 observational constraint is relative to the year 1992 and not 1961-1990
  r4 = median(windows[4,]) - (y.mod[240002] - y.mod[239992])

  #resid.y = c(r1, r2, r3, r4)
  resid.paleo = c(r1, r2, r3)
  resid.inst = r4
  sterr.y = obs.errs #This makes the model heteroskedastic
  
  #Calculate the likelihood. The observations are not correlated. They are independent
  llik.paleo = sum(dnorm(resid.paleo, mean=rep(0,3), sd = sqrt(var.paleo + sterr.y[1:3]^2), log=TRUE))
  llik.inst = sum(dnorm(resid.inst, mean=0, sd = sqrt(var.inst + sterr.y[4]^2), log=TRUE))
  
  llik = llik.paleo + llik.inst # assume residuals are independent
  llik  
}
log.pri = function(p)
{

  var.paleo = p[model.p+1]
  var.inst = p[model.p+2]

# var.paleo has inverse gamma prior, so there is a lower bound at 0 but no upper bound
  in.range = all(p > bound.lower) & all(p < bound.upper)
  
  alpha_var = 2
  beta_var = 1
  var_pri_paleo = 0
  var_pri_inst = 0
  
  if(in.range) {
    var_pri_paleo = (-alpha_var - 1)*log(var.paleo) + (-beta_var/var.paleo)
      
    lpri = 0 + var_pri_paleo + var_pri_inst
  } else {
    lpri = -Inf
  }

  lpri
}

# (log) posterior distribution:  posterior ~ likelihood * prior
log.post = function(p)
{  
  lpri = log.pri(p)
  if(is.finite(lpri)) { # evaluate likelihood if nonzero prior probability
    lpost = log.lik(p) + lpri
  } else {
    lpost = -Inf  
  }
  lpost
}
 

# calculate a scale matrix to transform a iid normal distribution,
# which defines a multivariate normal proposal distribution for MCMC
# (the proposal distribution is proportional to the covariance of
# a preliminary Markov chain of the posterior distribution to sample,
# tuned to be optimally scaled if the posterior is multivariate normal,
# plus a small nugget term to ensure ergodicity)
#
# from Gareth O. Roberts and Jeffrey S. Rosenthal,
# "Examples of adaptive MCMC", unpublished
proposal.matrix = function(prechain, mult=1, beta=0.05)
{
  # mult = overall scale factor to adjust all step sizes
  # beta = relative influence of nugget term
  
  p = ncol(prechain)
  precov = cov(prechain)
  
  propcov = (1-beta)*2.38^2*precov/p + beta*0.1^2*diag(p)/p
  propcov = mult*propcov
  
  mat = t(chol(propcov))
}

