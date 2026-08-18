# Interference
Code and Data for Manuscript: Predator interference imposes similar foraging time costs compared to handling time across taxa.

This repository contains all code and data required to reproduce the analyses presented in Biagioli et al. (2026). 

The workflow combines Julia (for fitting functional response models) and R (for secondary analyses, abundance scaling, 
regressions, and plotting).

----------------------------------------------------------------------
OVERVIEW OF THE WORKFLOW
----------------------------------------------------------------------

1) Julia is used to fit functional response models to FoRAGE data and estimate posterior distributions of model 
parameters.

2) R is used to:
   - Merge fitted parameters with FoRAGE metadata
   - Estimate predator and prey abundances using mass–abundance scaling
   - Calculate derived quantities (e.g., w/ah, wC/ahR)
   - Perform Bayesian regressions
   - Generate figures for the main text

----------------------------------------------------------------------
JULIA ANALYSIS
----------------------------------------------------------------------

STEP 1: Fit functional responses

File: Fit_inter_Revision.jl

- Fits the Rogers Random Predator Equation (RPE) and related functional response models.
- Estimates posterior distributions of functional response parameters.

Required inputs:
- FoRAGE functional response curve data
- FoRAGE metadata file
- simulate_data_set_INT.jl
- inter_pkgs.jl (install required packages)

Outputs:
Running this script generates the following files:

1) Inter_fitted_parameters_Revision.csv
   - Median values of posterior parameter estimates for each study

2) posterior_h_matrix_Revision.csv
   - Posterior distributions of handling time (h)

3) posterior_scaled_a_matrix_Revision.csv
   - Posterior distributions of attack rate (a)

4) posterior_w_matrix_Revision.csv
   - Posterior distributions of interference (w)

----------------------------------------------------------------------
R ANALYSIS
----------------------------------------------------------------------

STEP 2: Merge fitted parameters with FoRAGE metadata

File: /Analysis/Interference_Dataframe_Revision.Rmd

- Combines Julia-derived parameter estimates with FoRAGE metadata.

Required inputs:
- Inter_fitted_parameters_Revision.csv
- FoRAGE metadata file

----------------------------------------------------------------------
STEP 3: Estimate mass–abundance scaling relationships

File: /Mass-Abundance/MassAbundance_BayesianReg.R

- Fits Bayesian log–log mass–abundance scaling relationships used in downstream analyses.

Required inputs:
- AbundanceScaling.csv from Hatton et al. (2019)

Outputs:
- Saves an R workspace image containing fitted scaling relationships (used in Step 4).

----------------------------------------------------------------------
STEP 4: Calculate expected predator and prey abundances

File: /Mass-Abundance Scaling/DataManipulation_Revision.R

- Uses mass–abundance scaling relationships to estimate expected predator and prey abundances for each FoRAGE study.

Outputs:
Running this script generates:

1) forage_modified_Revision.csv
   - FoRAGE metadata augmented with median, 10%, and 90% credible intervals for expected predator 
      and prey abundances

2) Abundance_prey_full.csv
   - Posterior distributions of expected prey abundance

3) Abundance_pred_full.csv
   - Posterior distributions of expected predator abundance

----------------------------------------------------------------------
STEP 5: Calculate posterior distributions of w/ah

File: /Analysis/Posterior_Prediction_Calculation_Revision.R

- Calculates posterior distributions of w/ah for each study using posterior samples of w, a, and h.

Required inputs:
- posterior_w_matrix.csv
- posterior_scaled_a_matrix.csv
- posterior_h_matrix.csv

Outputs:
  - w_ah_post.csv
  - Posterior distributions of w/ah for each study

- Also calculates posterior distributions of wC/ahR using:
  - Posterior distributions of parameter estimates (w, a, h) from Julia
  - Posterior distributions of expected predator and prey abundances (C and R)
    
Required inputs:
- Abundance_pred_full.csv
- Abundance_prey_full.csv
- forage_modified_bolker_Abund_CI.csv

Outputs:
1) wC_ahR_predict.csv
   - Posterior distributions of wC/ahR for each study

2)  forage_modified_Revision.csv
   - Updated to include median, 10%, and 90% credible intervals for wC/ahR

----------------------------------------------------------------------
STEP 6: Run Bayesian regressions

File: /Analysis/Inter_Regressions_Revision.Rmd

- Fits Bayesian multivariate linear regression models examining relationships between Interference (w) and:
  - Functional response parameters (a and h)
  - Predator and prey body masses

Required inputs:
- Original (unmodified) FoRAGE_Inter_Revision.csv
- Modified forage_modified_Revision.csv

Outputs:
- Saves an R workspace image containing regression fits used for plotting

----------------------------------------------------------------------
STEP 7: Generate regression plots

File: /Plotting/Inter_Plotting_Bolker_Abund_CI.Rmd

- Creates and saves regression plots presented in the main text.

Required inputs:
- Inter_Regression_fits.RData (from Step 7)

----------------------------------------------------------------------
STEP 8: Generate forest plots

File: /Plotting/Inter_Hand_Plotting_Revision.Rmd

- Creates forest plots of interference metrics.

Required inputs:
- Inter_Regression_fits.RData
- wC_ahR_predict.csv
- w_ah_post.csv

Outputs:
- Forest plots used in the main text
