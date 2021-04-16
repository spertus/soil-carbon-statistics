# soil-carbon-simulations

This is a repository with code to run simulations and analyses to estimate soil organic carbon (SOC). 

`functions.R` contains the key functions to simulate and estimate SOC and to compute optimal sampling and measurement designs given various parameters and costs. 

`composites_analysis.R` contains an analysis of optimal sampling and measurement schemes for the in progress paper "Optimal sampling and measurement for estimating soil organic carbon stocks" by Spertus. 

`tautges_replication.R` contains code to replicate the analysis and run permutation tests on data from the Russell Ranch experiment presented in Tautges et al "Deep soil inventories reveal that impacts of cover crops and compost on soil carbon sequestration differ in surface and subsurface soils."  

`soil_carbon_measurements.Rmd` and `soil_carbon_permutations.Rmd` are R markdown notebooks analyzing measurement errors in SOC assays and permutation analyses of SOC experiments, respectively.  

`hedgerow_analysis.R` is an analysis of Jessica Chiartas' hedgerow/row quasi-experiment in the Central Valley.

`california_SOC_analysis.R` is a script to analyze Paige Stanley's samples from Paicines Ranch, Jessica Chiartas' samples from Southern California cropland, and assay replicates done on their respective labs. The corresponding data is `R_Heterogeneity_Master_PS_04132021.xlsx`. 
