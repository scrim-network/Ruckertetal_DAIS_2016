%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  - file = "DAIScali_hetero_model_iid_mcmc.m"
%  - Code written: August 2015, updated March 2016
%  - Author: Kelsey Ruckert (klr324@psu.edu)
%
%  -This program runs a Markov Chain Monte Carlo analysis of the DAIS model with a seed of 1234
%       assuming heteroskedastic errors and IID residuals as described in Ruckert et al. (2016).
%       For further description and references, please read the paper.
%
% THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
% NOT IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
% BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
% APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
% AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

rand('seed', 1234);  % For reproducibility
%rand('seed', 1780);  % For reproducibility

% step 1 define the boundary for parameters
run DAIS_data.m

%-------------------- SET INITIAL PARAMETERS --------------------%
% We will set the initial parameters to specifications from Shaffer [2014]
%Define the parameters:
% [1] gamma = 2 			  %sensitivity of ice flow to sea level
% [2] alpha = 0.35 			  %sensitivity of ice flow to ocean subsurface temperature
% [3] mu = 8.7    			  %Profile parameter related to ice stress [m^(1/2)]
% [4] eta = 0.012   		  %Relates balance gradient to precipitation [m^(-1/2) yr^(-1/2)]
% [5] Po = 0.35    			  %Precipitation at 0C [m of ice/yr]
% [6] kappa = 0.04 			  %Relates precipitation to temperature [K^-1]
% [7] fo = 1.2                %Constant of proportionality for ice speed [m/yr]
% [8] ho = 1471               %Initial value for runoff line calculation [m]
% [9] co = 95                 %Second value for runoff line calculation [m]
% [10] bo = 775               %Height of bed at the center of the continent [m]
% [11] s = 0.0006             %Slope of the bed

%Create a matrix that has 4 columns and 240300 rows
project_forcings(1:240300,1)=TA; project_forcings(1:240300,2)=TO;
project_forcings(1:240300,3)=GSL; project_forcings(1:240300,4)=SL;

% 240010x4 matrix
hindcast_forcings(1:240010,1)=TA(1:240010); hindcast_forcings(1:240010,2)=TO(1:240010);
hindcast_forcings(1:240010,3)=GSL(1:240010); hindcast_forcings(1:240010,4)=SL(1:240010);


% Best Case (Case #4) from Shaffer (2014)
IP = [2, 0.35, 8.7, 0.012, 0.35, 0.04, 1.2, 1471, 95, 775, 0.0006];

%Source the function with the standards and the initial parameters (IP) to
%get the best estimated AIS volume loss with respect to the present day in sea level equivalence (SLE):
standards = [Tice, eps1, del, eps2, TOo, Volo, Roa, R];
endinfo = [endsat, enddate];

%source("Scripts/DAIS_IceFlux_model.R")
AIS_melt = DAIS_IceFlux_model(IP, hindcast_forcings, standards, endinfo);

%set the end dates to the year 2300 to get future projections
endsat = 240298;
enddate = 240300;
endinfo = [endsat, enddate];
Project_melt = DAIS_IceFlux_model(IP, project_forcings, standards, endinfo);

%Set the end dates back to the hindcast period:
endsat = 240000;
enddate = 240010;
endinfo = [endsat, enddate];

%--------------------  CALCULATE RESIDUALS (PRIOR SIGMA) --------------------%
% These windows are presented in Shaffer (2014) and calculated from Shepherd et al. (2012)
% # Calculate the uncertainty with the +/- 2 standard error

% Create a vector with each observation year
% 120kyr, 20Kyr, 6kyr, 2002
obs_years = [120000, 220000, 234000, 240002];

% 1992 to 2011 trend from Shepherd et al. 2012 is -71 +/- 53 Gt per yr
% We want the cumulative sea-level equivalent in meters for the year 2002
% 360Gt = 1mm SLE
estimate_SLE_rate = abs(-71/360)/1000;
time_years = 2002-1992;
mid_cum_SLE_2002 = estimate_SLE_rate*time_years;

estimate_SLE_error = abs(-53/360)/1000; %1- sigma error
SE2_2002 = estimate_SLE_error*2; %2-sigma error

positive_2SE = mid_cum_SLE_2002 + SE2_2002; % Add the 2 standard error to the mean value
negative_2SE = mid_cum_SLE_2002 - SE2_2002; % Subtract the 2 standard error to the mean value

upper_wind = [6.0, -6.9, -1.25, positive_2SE];
lower_wind = [1.8, -15.8, -4.0, negative_2SE];
windows(1:4,1) = lower_wind;
windows(1:4,2) = upper_wind;

obs_errs = [abs(median(windows(1,:))-windows(1,1)); abs(median(windows(2,:))-windows(2,1));
            abs(median(windows(3,:))-windows(3,1)); SE2_2002];

%Set up equation to find the residuals and then the prior sigma
resid(1:length(obs_years)) = NaN;     %Create a vector of the residuals

for i = 1:length(obs_years)
resid(i) = (median(windows(i,:))-(AIS_melt(obs_years(i))-mean(AIS_melt(SL_1961_1990))));
end

sigma = std(resid)^2; %calculate the variance (sigma^2);

%%
%--------------------  SETUP MCMC --------------------%
%Set up the priors; the upper and lower bounds
bound_lower = IP - (IP*0.5)    ; bound_upper = IP + (IP*0.5);
bound_lower(1:2) = [1/2, 0]   ; bound_upper(1:2) = [17/4, 1]; %Set bounds for gamma and alpha
bound_lower(10:11) = [725, 0.00045]   ; bound_upper(10:11) = [825, 0.00075]; %Set bounds for bo and s

%bound_lower(12) = 0 ; bound_upper(12) = 1; % Prior uniform range for sigma (the variance)

%Shaffer [2014] best guess parameters
p = [IP, sigma];
p0 = [2.1, 0.29, 8, 0.015, 0.4, 0.04, 1.0, 1450, 90, 770, 0.0005, 0.8];

% Set up the structure for likelihood function
constraint.LIG = median(windows(1,:));
constraint.LGM = median(windows(2,:));
constraint.MH = median(windows(3,:));
constraint.IP = median(windows(4,:));
constraint

% SL_1961_1990 = 239961:239990;

%print out obs_errs
obs_errs

% Setup hindcast forcings to be called during MCMC calibration
global data; 
data = hindcast_forcings;

% Call physical and statistical model
logmodelprior=@log_pri_copy; % prior.
loglike=@log_lik_calibration_copy; % log likelihood.

%----------------- Estimate good initial parameters ------------%
%Initial parameters from optim in R
load OtimizedInitialParameters.txt;

minit(1,:) = OtimizedInitialParameters(:,2);
minit(12)=[1.16];
% minit=[2.2413, 0.3682, 8.0021, 0.0156, 0.4647, 0.0486, 1.2850, 1558.7492, 89.9974, 805.6928, 0.0006, 0.5782];
nsimu = 1000;  % number of simulations
step = minit/150;
thin=5;

%%
%-------------- Run the MCMC chain --------------%
rand('seed', 1234);  % For reproducibility
%rand('seed', 1780);  % For reproducibility

mmc=mcmc(minit,loglike,logmodelprior,step,nsimu,thin);
mmc1 = mmc;

%Estimate a more appropriate step size:
scale=proposal_matrix(mmc1,0.5);
step2(1:12) = diag(scale);

% New starting value:
minit2 = mmc1(length(mmc1),:);

rand('seed', 1234);  % For reproducibility
%rand('seed', 1780);  % For reproducibility

nsimu2 = 1.2e6;
mmc = mcmc(minit2, loglike, logmodelprior, step2, nsimu2, thin);
mmc2 = mmc;

save('DAIS_MCMCchain_1234', 'mmc1','mmc2')
%save('DAIS_MCMCchain_1780', 'mmc1','mmc2')

%--------------------  Analysis of the MCMC chain produced --------------------%
% m(1:100,:)=[]; %crop drift
% plotmatrix(m);
%
%% Plot some figures with the chain:
%figure(1); clf
%mcmcplot(mmc2,[],[],'chainpanel')

%figure(2); clf
%mcmcplot(mmc2,[],[],'dens')

%figure(3); clf
%mcmcplot(mmc2,[],[],'hist')

%save('DAIS_Matlab_MCMCcalibration','-v7.3') %Save the workspace into a R readable file

% The rest of the analysis will be run in R
%----------------------------------------------------------------------%

