library(tidyverse)
library(data.table)
library(gstat)
################## function to simulate a % SOC surface #################
simulate_truth <- function(i = 1, size = c(250,600), nugget = .01, sill = .05, range = 20, intercept = .01, y_trend = TRUE, max_mean = .2){
  x <- 1:size[1]
  y <- 1:size[2]
  xy <- expand.grid(x, y)
  colnames(xy) <- c("x", "y")
  
  #variogram (models covariance)
  vario <- vgm(nugget = nugget, psill = sill, range = range, model = "Exp")
  
  #should there be a trend in the mean concentration in the y-direction (e.g. a slope)?
  if(y_trend){
    #define intercept and slope on normal CDF scale
    beta <- c(qnorm(intercept), 0, (qnorm(max_mean) - qnorm(intercept)) / size[2])
    gstat_mod <- gstat(formula = z ~ 1 + x + y, locations = ~ x + y, dummy = TRUE, beta = beta, model = vario, nmax = 10)
  } else{
    #define intercept on normal CDF scale
    beta <- qnorm(intercept)
    gstat_mod <- gstat(formula = z ~ 1, locations = ~ x + y, dummy = TRUE, beta = beta, model = vario, nmax = 10)
  }
  
  
  simulation <- predict(gstat_mod, newdata = xy, nsim = 1) %>%
    mutate(z = pnorm(sim1)) %>%
    dplyr::select(x, y, z)
  simulation
}



##### function to collect samples from a simulated surface ####
#inputs:
#surface: a simulated surface, i.e. a matrix or array with SOC concentrations at each point
#outputs: samples collected along random transects, these are the true values (i.e. with no measurement error)  
collect_sample <- function(surface, transect = TRUE, n_samp = 9){
  #start randomly from lower left or lower right hand corner
  x_grid <- 1:max(surface$x)
  y_grid <- 1:max(surface$y)
  
  if(transect){
    #start randomly in lower left corner
    start_x <- sample(1:floor(max(surface$x)/n_samp), size = 1)
    start_y <- sample(1:floor(max(surface$y)/n_samp), size = 1)
    x_increment <- floor(max(surface$x)/n_samp)
    y_increment <- floor(max(surface$y)/n_samp)
    x_grid <- round(seq(start_x, x_increment*n_samp, length.out = n_samp))
    y_grid <- round(seq(start_y, y_increment*n_samp, length.out = n_samp))
    #also randomly start from lower left or right corner
    if(sample(c(TRUE, FALSE), size = 1)){
      x_grid <- rev(x_grid)
    }
    grid <- data.frame("x" = x_grid, "y" = y_grid)
  }
  
  true_samples <- surface %>% 
    inner_join(grid, by = c("x","y")) %>% 
    pull(z)
    
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
  #treatment: a binary vector of length ncol(measurements) indicating whether a plot was treatment (1) or control (0)
#outputs:
#a vector with the following components
  #control_mean_SOC: intercept in ANOVA
  #ATE: estimate of average treatment effect, based on coefficient in ANOVA
  #ATE_se: the estimated standard error of the ATE estimate
  #pval: the p value of a test of whether the ATE is different from 0
estimate_ATE_anova <- function(measurements, treatment){
  mean_pct_SOC <- colMeans(measurements)
  ATE_anova <- lm(mean_pct_SOC ~ treatment)
  c(
    "control_mean_SOC" = coef(ATE_anova)[1],
    "ATE" = coef(ATE_anova)[2],
    "ATE_se" = summary(ATE_anova)$coefficients[2,2],
    "pval" = summary(ATE_anova)$coefficients[2,4]
  )
}



################## function to run simulations ################
#inputs: 
  #n_sims: number of simulated samples to draw
  #surfaces: a collection of surfaces (i.e. plots) as a list
  #pct_perturbations: a vector of perturbations (i.e. measurement errors)
  #treatment: a binary vector indicating whether a plot is treated (1) or control (0)
#outputs:
  #a matrix of results with columns for control SOC estimate, treatment effect, estimated standard error, and a p value
run_sims <- function(n_sims = 1000, surfaces, pct_perturbations, treatment){
  results <- matrix(0, ncol = 4, nrow = n_sims)
  for(i in 1:n_sims){
    samples <- surfaces %>% 
      map(collect_sample) %>%
      reduce(cbind)
    #measured_samples <- samples %>%
    #  map(perturbed_measurements, pct_perturbations = pct_perturbations) %>%
    #  reduce(cbind)
    results[i,] <- estimate_ATE_anova(measurements = samples, treatment = treatment)
  }
  results
}


