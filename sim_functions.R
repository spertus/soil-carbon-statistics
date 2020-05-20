library(tidyverse)
library(data.table)
library(gstat)
################## function to simulate a % SOC surface #################
simulate_truth <- function(size = c(250,600), nugget = .01, sill = .05, range = 20, intercept = .01, y_trend = TRUE, max_mean = .2){
  #generates a 2-dimensional surface based on a variogram model, possibly with a trend
  #relies on the gstat::gstat and gstat::vgm to sample from a Gaussian random field, then transforms to a copula (bounded [0,1] variable) using pnorm()
  #inputs:
    #size: a length-2 vector with x and y dimensions. Defaults to c(250,600) which corresponds to a 25 x 60 meter plot with resolution at 10cm
    #nugget: (co)variance parameter, the nugget variance (distance 0 variance)
    #sill: (co)variance parameter, the sill variance (maximum variance)
    #range: (co)variance parameter, how far until sill variance is reached
    #intercept: the mean, if y_trend is FALSE. If not, the mean at y = 0
    #y_trend: boolean, should there be a trend in the y direction (e.g. to simulate a hill slope)?
    #max_mean: mean at y == max(size[2]) if y_trend == TRUE, if y_trend == FALSE then ignored
  #outputs:
    #a dataframe with 3 columns: x (the x coordinate), y (the y coordinate), and z (the value of % SOC)
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

############ function to plot a surface as a heatmap###############
plot_surface <- function(surface){
  ggplot(data = surface, aes(x = x, y = y, fill = z)) + geom_raster()
}

############### function to add samples at depth, given surface samples ##############
add_depth_samples <- function(surface, increments = 4, decrease = "exponential", gamma = -1, proportions){
  #SOC stocks are defined at depth, typically up to a meter
  #standard increments for depth for SOC sequestration measurement are: 0-10cm, 10-30cm, 30-50cm, and 50-100cm (see eg Ryals et al 2014)
  #some studies have found that soil carbon declines exponentially with depth (e.g Allen et al 2010 "a review of sampling designs...")
  #this function adds new surfaces at depth.
  #SOC at depth is completely determined by SOC at the surface (no randomness) with user-defined decreases
  #inputs:
  #surface: a simulated surface in (x,y) space, with %SOC at each location
  #increments: the number of depth increments, depths 1:increments will be defined and returned
  #decrease: how should the %SOC decrease with depth
  #"exponential": SOC at depth = (SOC at surface) * exp(gamma * depth) 
  #"manual": SOC at depth = (SOC at surface) * proportions[depth]
  #gamma: relevant if decrease == "exponential", a tuning parameter for exponential decay of SOC. Positive gammas mean SOC increases with depth 
  #proportions: relevenat if decrease == "manual", explicitly define proportions at lower depth. The first proportion (the surface) should always be 1
  stopifnot(decrease == "exponential" | decrease == "manual")
  
  surface_SOC <- surface %>% pull(z)
  if(decrease == "exponential"){
    index <- 0:(increments-1)
    proportions <- exp(gamma * index)
  }
  proportions_matrix <- matrix(proportions, ncol = increments, nrow = length(surface_SOC), byrow = TRUE)
  z_depth <- as.data.frame(surface_SOC * proportions_matrix)
  colnames(z_depth) <- paste("z", 1:increments, sep = "")
  surface_depth <- surface %>% 
    dplyr::select(x,y) %>%
    bind_cols(z_depth) %>%
    pivot_longer(cols = starts_with("z"), names_to = "depth", names_prefix = "z", values_to = "z")
}





########## function to collect samples from a simulated surface ###########
collect_sample <- function(surface, design = "transect", n_samp = 9){
  #inputs:
    #surface: a simulated surface, i.e. a matrix or array with SOC concentrations at each point
    #design: the sampling design, currently:
      #"transect" to sample along a transect starting randomly in either the lower left or right corner
      #"simple random sample" for a simple random sample from the entire surface
    #n_samp: the number of samples to draw
  #outputs: samples collected along random transects, these are the true values (i.e. with no measurement error)  
  if(design == "transect"){
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
  } else if(design == "simple random sample"){
    y_samples <- sample(1:max(surface$y), size = n_samp, replace = FALSE)
    x_samples <- sample(1:max(surface$x), size = n_samp, replace = FALSE)
    grid <- data.frame("x" = x_samples, "y" = y_samples)
  } else{
    stop("need to specify a valid sampling design")
  }
  
  true_samples <- surface %>% 
    inner_join(grid, by = c("x","y")) %>% 
    pull(z)
    
  true_samples
}





########### function to perturb measurements of samples ##########
perturbed_measurements <- function(true_samples, error_type, error_bounds){
  #adds independent, uniform measurement error to samples
  #inputs:
    #true_samples: true pct carbon (unknown) to be measured with error
    #error_type: how should the error perturb the true vale?
      #"additive": the errors are just added to the true values (unbiased implies centered at 0)
      #"multiplicative": the true values are dilated by the errors (unbiased implies centered at 1)
    #error_bounds: a length-2 vector specifying the lower and upper bounds of the (uniform error)
  
  error <- runif(n = true_samples, min = error_bounds[1], max = error_bounds[2])
  if(error_type == "multiplicative"){
    measured_samples <- true_samples * error
  } else if(error_type == "additive"){
    measured_samples <- true_samples + error
  } else{
    stop("input a valid error_type")
  }
  measured_samples
}




