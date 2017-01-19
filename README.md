#Sea-level contribution from the Antarctic Ice Sheet code and pre-run files for Ruckert et al. (2017)

README file last updated by Kelsey Ruckert, klr324-at-psu-dot-edu, Wed Jan 4 2017.

##Citation

This code and pre-run data analysis are intended to accompany the results of

>Ruckert, KL, Shaffer, G, Pollard, D, Guan, Y, Wong, TE, Forest, CE and Keller, K (2017). Assessing the impact of retreat mechanisms in a simple Antarctic ice sheet model using Bayesian calibration,

The corrected code from Shaffer (2014) is provided as well. Please cite these papers when using any results generated with the Ruckert et al. (2017) or the Shaffer (2014) code. 

>Shaffer G (2014) Formulation, calibration and validation of the DAIS model (version 1), a simple Antarctic ice sheet model sensitive to variations of sea level and ocean subsurface temperature _Geoscientific Model Development_ **7**(4) 1803–1818, doi:10.5194/gmd-7-1803-2014.

##Overview
Ruckert et al. (2017) approximates *"the effects of Marine Ice Cliff Instability (MICI) by comparing our results to those from expert assessments with more realistic models and quantifies the bias during the last interglacial when MICI may have been triggered."* Accounting for retreat mechanisms is important because ignoring them can lead to not only a low-biased during warming period AIS melting, but potentially a low-bias in projected sea levels and flood risks.

This R code is intended to help users who wish to work with the sea-level projections or hindcasts from the Antarctic ice sheet in greater detail than provided in the text. Key functionality of these scripts include:

1. Antarctic ice sheet sea-level contribution from 240,000 years ago to the year 2300 with associated probabilities
2. How to calibrate the model to data constraints with heteroskedastic errors using Markov chain Monte Carlo
3. Produces plots from the paper

The pre-run analysis is also provided for those who do not want to run the codes (which take multiple hours to run) and wish to use the results from the paper.
**The saved workspace used in the paper analysis can be downloaded at this location [link coming soon].**

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
plotrix  

You can install and open the packages in R as shown in the example below, which installs the adaptMCMC package.

```R
install.packages("adaptMCMC")
library(adaptMCMC)
``` 

The original Shaffer (2014) scripts are written to run in matlab. For more details on these codes please see the "ReadMeDAIS" file within the "Shaffer_2014_DAIS_Code" folder.

## Instructions
* Download the Ruckert_2017_DAIS_Code folder.
* Open the `.R` files in R or Rstudio and edit the paths to the files according to your directory structure.
* Create `Converge` and 'Scratch' directories.

The Ruckert_2016_DAIS_Code directory contains all the scripts and data necessary to run the analysis along with a README file. _(Note that the user may have to edit the scripts according to their folder directory so that the scripts will locate the files/scripts needed to run. Additionally, the user will have to create directories named **Converge** and **Scratch**. The script will output files into these directories.)_

Instructions on how to run the scripts can be found in the README file in the Ruckert_2017_DAIS_Code directory. The main files for running the analysis  are the **DAIScali_hetero_model_iid_mcmc_R_C** (MCMC calibration) and **DAIS_precali_LHS_C** (Latin hypercube sampling) files. Note, the MCMC calibration and the optimization proccess in the LHS script take multiple hours to complete.

##Contact
Kelsey Ruckert  
E-mail: <klr324@psu.edu>  

Corresponding author:  
Klaus Keller   
E-mail: <klaus@psu.edu>
