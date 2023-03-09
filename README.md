# soil-carbon-statistics

This is a repository with code to run simulations and analyses to measure soil organic carbon (SOC). Note that a skeletal R package containing functions to run the nonparametric tests we developed to test for soil carbon sequestration can be installed easily by running `devtools::install_github("spertus/nptests")`. 

`functions.R` contains the key functions to simulate and estimate SOC and to compute optimal sampling and measurement designs given various parameters and costs. 

`composites_analysis.R` contains an analysis of optimal sampling and measurement schemes for the in progress paper "Optimal sampling and measurement for estimating soil organic carbon stocks" by Spertus. 

`tautges_replication.R` contains code to replicate the analysis and run permutation tests on data from the Russell Ranch experiment presented in Tautges et al "Deep soil inventories reveal that impacts of cover crops and compost on soil carbon sequestration differ in surface and subsurface soils."  

`soil_carbon_measurements.Rmd` and `soil_carbon_permutations.Rmd` are R markdown notebooks analyzing measurement errors in SOC assays and permutation analyses of SOC experiments, respectively.  

`hedgerow_analysis.R` is an analysis of Jessica Chiartas' hedgerow/row quasi-experiment in the Central Valley. The corresponding data files are `HR_Soil_Health_Updated.csv` and `HR_soil_carbon.csv`.

`california_SOC_analysis.R` is a script to analyze Paige Stanley's samples from Paicines Ranch, Jessica Chiartas' samples from Southern California cropland, and assay replicates done on their respective labs. The corresponding data file is `R_Heterogeneity_Master_PS_04132021.xlsx`. 

`saturation_hypothesis.Rmd` is an R markdown notebook evaluating the "saturation hypothesis" that SOC gains are moderated by baseline SOC levels. The hypothesis is articulated and critiqued in [Slessarev et al 2022](https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.16491).
