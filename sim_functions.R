library(tidyverse)

##### function to collect samples from a simulated surface ####
#inputs:
#surface: a simulated surface, i.e. a matrix or array with SOC concentrations at each point
#outputs: samples collected along random transects, these are the true values (i.e. with no measurement error)  
collect_sample <- function(surface){
  #start randomly from lower left or lower right hand corner
  left <- sample(c(TRUE, FALSE), size = 1)
  if(left){
    vertical_bounds <- c(nrow(surface)-ceiling(nrow(surface)/10), nrow(surface))
    horizontal_bounds <- c(1,floor(ncol(surface)/10))
    
    initial_location <- c(sample(vertical_bounds[1]:vertical_bounds[2], size = 1), sample(horizontal_bounds[1]:horizontal_bounds[2], size = 1))
    vertical_steps <- floor(seq(0, initial_location[1], length.out = 9))
    horizontal_steps <- floor(seq(0, ncol(surface)-initial_location[2], length.out = 9))
    
    transect <- matrix(c(initial_location[1] - vertical_steps, initial_location[2] + horizontal_steps), nrow = length(vertical_steps))
  } else {
    vertical_bounds <- c(nrow(surface)-ceiling(nrow(surface)/10), nrow(surface))
    horizontal_bounds <- c(ncol(surface),ncol(surface) - floor(ncol(surface)/10))
    
    initial_location <- c(sample(vertical_bounds[1]:vertical_bounds[2], size = 1), sample(horizontal_bounds[1]:horizontal_bounds[2], size = 1))
    vertical_steps <- floor(seq(0, initial_location[1], length.out = 9))
    horizontal_steps <- floor(seq(0, initial_location[2], length.out = 9))
    
    transect <- matrix(c(initial_location[1] - vertical_steps, initial_location[2] - horizontal_steps), nrow = length(vertical_steps))
  }
  true_samples <- simulated_surface[transect]
  true_samples
}



########### function to perturb measurements of samples ##########
#inputs:
  #true_samples: true pct carbon (unknown) to be measured with error
  #pct_perturbations: a vector of any length that gives magnitude of pct perturbations (probably empirically determined) to be randomly sampled. Sign is randomly sampled, so perturbations are symmetric about 0 (leading to unbiased but more variable estimates).
perturbed_measurements <- function(true_samples, pct_perturbations){
  random_perturbation_magnitude <- sample(pct_perturbations, size = length(true_samples), replace= TRUE)
  perturbations <- true_samples * random_perturbation_magnitude * sample(c(-1,1), size = length(true_samples), replace = TRUE)
  measured_samples <- true_samples + perturbations
  measured_samples
}

############ functions to compute plot averages and ATEs ############
#inputs:
  #measurements: a matrix with samples in the rows and plots in the columns
  #bulk density: a vector of bulk densities (as total mass of plot) of length 1 (all BDs are the same) or length ncol(measurements) (each plot has its own BD)
  #treatement: a binary vector of length ncol(measurements) indicating whether a plot was treatment (1) or control (0)
#outputs:
#a vector with the following components
  #control_mean_SOC: intercept in ANOVA
  #ATE: estimate of average treatment effect, based on coefficient in ANOVA
  #ATE_se: the estimated standard error of the ATE estimate
  #pval: the p value of a test of whether the ATE is different from 0
estimate_ATE_anova <- function(measurements, bulk_density, treatment){
  mean_pct_SOC <- colMeans(measurements)
  total_SOC <- mean_pct_SOC * bulk_density
  ATE_anova <- lm(total_SOC ~ treatment)
  c(
    "control_mean_SOC" = coef(ATE_anova)[1],
    "ATE" = coef(ATE_anova)[2],
    "ATE_se" = summary(ATE_anova)$coefficients[2,2],
    "pval" = summary(ATE_anova)$coefficients[2,4]
  )
}



################## function to run simulations ################

run_sims <- function(n_sims = 1000, surfaces, pct_perturbations, treatment, bulk_density){
  results <- matrix(0, ncol = 4, nrow = n_sims)
  for(i in 1:n_sims){
    samples <- surfaces %>% 
      map(collect_sample)
    measured_samples <- samples %>%
      map(perturbed_measurements, pct_perturbations = pct_perturbations) %>%
      reduce(cbind)
    results[i,] <- estimate_ATE_anova(measurements = measured_samples, bulk_density = bulk_density, treatment = treatment)
  }
}


