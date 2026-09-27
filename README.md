# Spatial Dependence Monte Carlo

R scripts for AEM 6850 (Empirical Methods, Cornell, Fall 2025) that use a Monte Carlo simulation to compare fixed-effects OLS, Conley standard errors and a spatial error model on spatially correlated county panel data.

script_files:
- spatial_se_sim.R : builds a 5-year panel of 958 eastern US counties with a 10-nearest-neighbor weights matrix, simulates spatially correlated x2 and errors, runs 750 simulations of OLS (plm), Conley standard errors (0-500 mile cutoffs) and a spatial error model (splm), and plots the estimates and standard errors
- functions.R : mrnorm for correlated normal draws and conley for spatial HAC standard errors
- out.RDS : saved simulation results, loaded by spatial_se_sim.R if present
- Spatial-Dependence-Monte-Carlo.Rproj : RStudio project, open this first so the scripts find data/ and output_figure/

data:
- CountyDistrictList.csv : list matching each county to its agricultural statistics district

output_figure:
- OLS_verses_SEM_Spatial_Dependence.png : estimates and standard errors for OLS, SEM and Conley

Other files:
- readme.rtf : full readme with general, methodological and data-specific information
