#Sea-level contribution from the Antarctic Ice Sheet code for Ruckert et al. (in prep.)

README file last updated by Kelsey Ruckert, klr324-at-psu-dot-edu, Fri July 8 12:19:51 EST 2016

##Citation

This code is intended to accompany the results of

>Ruckert, KL, Shaffer, G, Pollard, D, Guan, Y, and Keller, K. Neglecting cliff instability in Antarctic ice sheet models can reduce melting during warming periods, (in prep.).

Please cite that paper and the Shaffer (2014) study when using any results generated with this code.

>Shaffer G (2014) Formulation, calibration and validation of the DAIS model (version 1), a simple Antarctic ice sheet model sensitive to variations of sea level and ocean subsurface temperature _Geoscientific Model Development_ **7**(4) 1803–1818, doi:10.5194/gmd-7-1803-2014.

##Overview

This code requires R and matlab with the following libraries:
- coda
- lhs
- mcmc (matlab script included)
- R.matlab
- RColorBrewer
- ash
- fields
- pscl
- DEoptim
- compiler

This R code is intended to help users who wish to work with the sea-level projections or hindcasts from the Antarctic ice sheet in greater detail than provided in the text. Key functionality of these scripts include:

1. Antarctic ice sheet sea-level contribution from 240,000 years ago to 2300 with associated probabilities
2. How to fit the model to data constraints with heteroskedastic errors using Markov Chain Monte Carlo
3. Produces plots from the paper

The IFILES directory contains all the scripts and data necessary to run the analysis along with a README file. The prerun analysis output used to generate the Ruckert et al. (in prep.) figures exceeds 100 MB. For access to the prerun analysis please contact the corresponding author. _(Note that the folder directory MUST be in the same format as when downloaded otherwise the scripts will not locate the files/scripts needed to run. Additionally, the following empty folders need to be created before running the analysis: 'Workspace' and 'Figures' with the subfolder, 'SuppFigures'. Output will be saved to these folders.)_

The most important functions are **DAISobs_likelihood_iid**, **DAIScali_hetero_model_iid_mcmc**, and **DAIScali_hetero_model_iid_mcmcmat** for MCMC calibration and projections, **DAIS_precali_LHS** for Latin hypercube pre-calibration,  **DAIS_convergence_plotsR** for testing MCMC convergence, and **DAIS_IceFlux_model**, which is the Shaffer (2014) DAIS model.

To fit the DAIS model using MCMC, simply open a terminal and run **DAIScali_hetero_model_iid_mcmc** (Note, this runs ~24hrs). Then source **DAIScali_hetero_model_iid_mcmcmat** and **MCMC_plots** to generate projections and plots. (Specific details can be found in the README.txt file.

To run the Latin hypercube pre-calibration, simply open a terminal and source **DAIS_precali_LHS** followed by **LHS_plots** to generate projections and plots. (Note, this runs ~2hrs; specific details can be found in the README.txt file.

##Contact
Kelsey Ruckert: <klr324@psu.edu>  
Corresponding author: Klaus Keller at <klaus@psu.edu>
