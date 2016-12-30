#Sea-level contribution from the Antarctic Ice Sheet code for Ruckert et al. (Accepted)

README file last updated by Kelsey Ruckert, klr324-at-psu-dot-edu, Wed Dec 21 14:19:51 EST 2016

##Citation

This code is intended to accompany the results of

>Ruckert, KL, Shaffer, G, Pollard, D, Guan, Y, Wong, TE, Forest, CE and Keller, K. Assessing the impact of retreat mechanisms in a simple Antarctic ice sheet model using Bayesian calibration, (Accepted).

Please cite that paper and the Shaffer (2014) study when using any results generated with this code.

>Shaffer G (2014) Formulation, calibration and validation of the DAIS model (version 1), a simple Antarctic ice sheet model sensitive to variations of sea level and ocean subsurface temperature _Geoscientific Model Development_ **7**(4) 1803–1818, doi:10.5194/gmd-7-1803-2014.

##Overview
This study approximates *"the effects of Marine Ice Cliff Instability (MICI) by comparing our results to those from expert assessments with more realistic models and quantifies the bias during the last interglacial when MICI may have been triggered."* Accounting for retreat mechanisms is important because ignoring them can lead to not only a low-biased during warming period AIS melting, but potentially a low-bias in projected sea levels and flood risks.

This R code is intended to help users who wish to work with the sea-level projections or hindcasts from the Antarctic ice sheet in greater detail than provided in the text. Key functionality of these scripts include:

1. Antarctic ice sheet sea-level contribution from 240,000 years ago to the year 2300 with associated probabilities
2. How to calibrate the model to data constraints with heteroskedastic errors using Markov chain Monte Carlo
3. Produces plots from the paper

## Requirements
### Software
These scripts are written to run in R (tested under R v3.2.1; https://www.r-project.org/) with some help from C. They also require mulitple packages including:  
>adaptMCMC 
coda  
lhs  
pscl  
compiler  
mcmc  
ash  
fields  
RColorBrewer

You can install and open the packages in R as shown in the example below, which installs the lhs package.

```R
install.packages("adaptMCMC")
library(adaptMCMC)
```

## Instructions
* Download the IFILES folder.
* Open the `.R` files in R or Rstudio and edit the paths to the files according to your directory structure.

The IFILES directory contains all the scripts and data necessary to run the analysis along with a README file. _(Note that the user may have to edit the scripts according to their folder directory so that the scripts will locate the files/scripts needed to run. Additionally, the following empty folders need to be created before running the analysis: 'Workspace', 'Figures', and 'Converge'. Output will be saved to these folders.)_

Instructions on how to run the scripts can be found in the README file in the IFILES directory. The main files for running the analysis  are the **DAIScali_hetero_model_iid_mcmc_R_C** (MCMC calibration) and **DAIS_precali_LHS_C** (Latin hypercube sampling) files. Note, the MCMC calibration takes ~5.5 hours to complete.

##Contact
Kelsey Ruckert  
E-mail: <klr324@psu.edu>  

Corresponding author:  
Klaus Keller   
E-mail: <klaus@psu.edu>
