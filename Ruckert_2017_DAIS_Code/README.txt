"Assessing the impact of retreat mechanisms in a simple Antarctic ice sheet model using Bayesian calibration" Codes

Reference:
Ruckert, KL, Shaffer, G, Pollard, D, Guan, Y, Wong, TE, Forest, CE, and Keller, K (2017). Assessing the impact of retreat mechanisms in a simple Antarctic ice sheet model using Bayesian calibration, PLOS ONE 12(1): e0170052. doi: 10.1371/journal.pone.0170052.

Copyright 2016 Kelsey Ruckert (klr324@psu.edu)
This file is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This file is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>. 

Required software:
    R
    C

Required libraries:
    coda
    lhs
    mcmc
    RColorBrewer
    ash
    fields
    pscl
    compiler
    adaptMCMC
    plotrix

Please note that file locations may need to be edited according to the users structure of the directory.
Additionally, folders for output need to be created otherwise the R with produce errors.
The saved workspace/files used in the paper analysis can be downloaded at this location: 
https://download.scrim.psu.edu/Ruckert_etal_DAIS/.
 =============================================================================
| Create Output Directories:
|(NOTE: The directories created here will contain the files outputed from the scripts.)
|
|  1. Navigate to the directory with the codes
|  2. Create Converge directory
|  3. Create Scratch directory
 =============================================================================
 =============================================================================
| How to run LHS analysis:
|(NOTE:LHS analysis takes 2+ hrs to run. This is due to the optimization command and not LHS sampling)
|
|  1. Open R
|  2. Run DAIS_precali_LHS_C.R
|  
|Create Plots:
|  3. Run LHS_plots_C.R
 =============================================================================
| How to run MCMC analysis:
|(NOTE: MCMC analysis takes a ~5 to 6 hrs to run)
|(NOTE: The seed was not set for the MCMC analysis used in the paper. Even though the chains are converged, 
| estimates may vary ever so slightly because MCMC is a random number generator. This does not impact    
| the overall results and conclusions. The seeds 1234 and 1780 should yield the same distributions and figures.
|
| If you wish to use the exact workspace/files used in the paper analysis, then they can be downloaded at 
| this location: https://download.scrim.psu.edu/Ruckert_etal_DAIS/.)
|
|  1. Open R
|  2. Run DAIScali_hetero_model_iid_mcmc_R_C.R
|
|Check For Convergence:
|  3. Run DAIScali_hetero_model_iid_mcmc_R_C1780.R
|  4. Run DAIS_convergence_plotsR_C.R    
|  
|Create Plots:
|  8. Run MCMC_plots_C.R
 =============================================================================



 =============================================================================
| Short description of main scripts:
|
|  1. DAIS_data_C.R: Input forcing data, standards, and other important info.
|  2. dais.c: Antarctic ice sheet model.
|  3. DAIS_precali_LHS_C.R: Latin hypercube sampling of the AIS model.
|  4. DAIScali_hetero_model_iid_mcmcR_C.R: Heteroskedastic IID MCMC fitting of the AIS model.      
|  5. DAIS_convergence_plots_C.R: Tests MCMC convergence based on the Potential Scale Reduction 
|     factor, Trace plots, and Heidelberger and Welch's convergence diagnostic.
|  6. DAISobs_likelihood_iid.R: Estimates the likelihood functions for the MCMC methods.
|  7. LHS_plots_C.R & MCMC_plots_C.R: Generates figures shown in the paper.
 =============================================================================


Credits:
    DAIS Model: 
Shaffer, G.: Formulation, calibration and validation of the DAIS model (version 1), a simple Antarctic ice sheet model sensitive to variations of sea level and ocean subsurface temperature, Geoscientific Model Development, 7(4), 1803–1818, 2014.

Please enjoy the code and offer us any and all of your suggestions. Thanks!
Kelsey Ruckert (klr324@psu.edu)

