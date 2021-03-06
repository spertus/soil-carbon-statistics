library(tidyverse)
library(data.table)
library(gstat)
library(BalancedSampling)

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
  colnames(z_depth) <- paste("z", 0:(increments-1), sep = "")
  surface_depth <- surface %>% 
    dplyr::select(x,y) %>%
    bind_cols(z_depth) %>%
    pivot_longer(cols = starts_with("z"), names_to = "depth", names_prefix = "z", values_to = "z")
}

########## function to draw a stratified sample from a vector ##########
get_stratified_sample <- function(population_vector, within_strata_sample_size, strata_endpoints){
  #inputs:
  #population_vector: a vector of numbers to be sampled
  #within_strata_sample_size: how many samples to draw from each strata? Either a single number (equal sample sizes) or a vector of length(strata_endpoints) + 1
  #strata_endpoints: the endpoints of the strata in terms of the entries (not indices) of the population_vector. 
  #The smallest entry is non-inclusive, e.g. if you want to include 1 the smallest endpoint should be 0.
  #outputs:
  #a vector of samples, the names of the entries indicate the strata
  strata <- cut(population_vector, strata_endpoints)
  if(anyNA(strata)){
    stop("Not all elements of population_vector are in a strata! strata_endpoints do not make a partition. Note that left strata endpoint is non-inclusive.")
  }
  if(length(unique(table(strata))) != 1){
    warning(paste("Not all strata are the same size. Sizes are:", paste(table(strata), collapse = " ")))
  }
  samples <- tapply(X = population_vector, INDEX = strata, FUN = sample, size = within_strata_sample_size, replace = FALSE)
  samples_vec <- unlist(samples)
  names(samples_vec) <- str_sub(names(samples_vec), end = -2)
  samples_vec
}



########## function to collect samples from a simulated surface ###########
collect_sample <- function(surface, design = "transect", n_samp = 9, n_strata = NULL){
  #inputs:
    #surface: a simulated surface, i.e. a matrix or array with SOC concentrations at each point
    #design: the sampling design, currently:
      #"transect" to sample along a transect starting randomly in either the lower left or right corner
      #"simple random sample" for a simple random sample from the entire surface
      #"well-spread" for a spatially balanced sample as per Grafström, A. Lisic, J (2018). BalancedSampling: Balanced and Spatially Balanced Sampling
      #"stratified random sample" for a stratified random sample. As of now, the strata are defined only in the y direction and are of equal size max(surface$y) / n_strata
    #n_samp: the number of samples to draw
  #outputs: samples collected along random transects, these are the true values (i.e. with no measurement error)  
  if(design == "transect"){
    #start randomly in lower left corner
    start_x <- sample(1:floor(max(surface$x)/n_samp), size = 1)
    start_y <- sample(1:floor(max(surface$y)/n_samp), size = 1)
    x_increment <- floor(max(surface$x) / n_samp)
    y_increment <- floor(max(surface$y) / n_samp)
    x_grid <- start_x + (0:(n_samp-1)) * x_increment
    y_grid <- start_y + (0:(n_samp-1)) * y_increment
    #also randomly start from lower left or right corner
    if(sample(c(TRUE, FALSE), size = 1)){
      x_grid <- rev(x_grid)
    }
    grid <- data.frame("x" = x_grid, "y" = y_grid, "strata" = paste("(", min(surface$y), ",", max(surface$y), "]", sep = ""))
  } else if(design == "simple random sample"){
    y_samples <- sample(1:max(surface$y), size = n_samp, replace = FALSE)
    x_samples <- sample(1:max(surface$x), size = n_samp, replace = FALSE)
    grid <- data.frame("x" = x_samples, "y" = y_samples, "strata" = paste("(", min(surface$y), ",", max(surface$y), "]", sep = ""))
  } else if(design == "stratified random sample"){
    x_samples <- sample(1:max(surface$x), size = n_samp, replace = FALSE)
    if((n_samp %% n_strata) != 0){stop("sample size is not a multiple of the number of strata -> unequal sampling across strata (not currently supported)")}
    strata_endpoints <- round(seq(0, max(surface$y), length.out = n_strata + 1))
    y_samples <- get_stratified_sample(1:max(surface$y), within_strata_sample_size = n_samp / n_strata, strata_endpoints = strata_endpoints)
    
    grid <- data.frame("x" = x_samples, "y" = y_samples, "strata" = names(y_samples))
  } else if(design == "well-spread"){
    inclusion_probs <- rep(n_samp / nrow(surface), nrow(surface))
    sample_rows <- SamplingBigData::lpm2_kdtree(prob = inclusion_probs, x = cbind(surface$x, surface$y))
    grid <- data.frame("x" = surface$x[sample_rows], "y" = surface$y[sample_rows], "strata" = paste("(", min(surface$y), ",", max(surface$y), "]", sep = ""))
  } else{
    stop("need to specify a valid sampling design")
  }
  
  true_samples <- surface %>% 
    inner_join(grid, by = c("x","y")) 
  true_samples
}


######### function to plot a surface along with the samples taken from it ######3
plot_surface_samples <- function(surface, samples){
  ggplot(data = surface, aes(x = x, y = y, fill = z)) + 
    geom_raster() + 
    geom_point(data = samples, aes(x = x, y = y), size = 1.3, colour = "red")
}


########### function to composite samples #############
composite_samples <- function(samples, k = 5){
  #this function reduces a set of samples to a smaller set wherein each element is the average over groups of the original sample
  #assumes perfect mixing and equal weights of composited material
  #groupings are determined by samples, with n/k samples in each group (except for possibly the last, see below)
  #if (length(samples) %% k) != 0, then the last group will accomodate all the additional samples
  #input:
    #samples: a dataframe of samples from a plot, as output by collect_sample()
    #k: the number of samples we want out ("composite from n down to k") per strata. For simple or transect sampling, this is the total number of samples we want out.
  #output: a dataframe of composite samples.
  strata_size <- table(samples$strata)[1]
  num_strata <- nrow(samples) / strata_size
  if((strata_size %% k) != 0){stop("n is not divisible by k!")}
  index <- rep(1:k, num_strata)
  
  composite_frame <- samples %>%
    mutate(composite_group = index) %>%
    group_by(strata, composite_group) %>%
    summarize(composite_size = n(), z = mean(z))
}


########### function to measureme samples with optional error ##########
measure_samples <- function(true_samples, error_type = "multiplicative", error_bounds, error_sd, replicates = 1){
  #corrupts samples with independent, symmetric, beta distributed measurement error 
  #to measure samples without error, set error_sd = 0
  #inputs:
    #true_samples: dataframe of samples as output by composite_samples
    #error_type: how should the error perturb the true vale?
      #"additive": the errors are just added to the true values (unbiased implies centered at 0)
      #"multiplicative": the true values are dilated by the errors (unbiased implies centered at 1)
    #error_bounds: a length-2 vector specifying the lower and upper bounds of the error, the mean of these (halfway between the left and right bound) is the expected value of the measurement error
    #replicates: the number of times to measure each sample (duplicate, triplicate, etc)
  #output: 
    #the measured samples, a dataframe with column for sample #, replicate, and measurement value
  if(error_sd^2 > (1/4)*(error_bounds[2] - error_bounds[1])^2){
    stop("error variance is too big given the error bounds!")
  }
  alpha <- (error_bounds[1] - error_bounds[2])^2 / (8 * error_sd^2) - 1/2
  delta_star <- rbeta(nrow(true_samples)*replicates, shape1 = alpha, shape2 = alpha)
  delta <- (delta_star - 1/2) * abs(error_bounds[1] - error_bounds[2]) + mean(error_bounds)
  samples_frame <- expand.grid("sample" = rep(1:nrow(true_samples)), "measurement_replicate" = rep(1:replicates))
  
  if(error_type == "multiplicative"){
    measured_samples <- rep(true_samples$z, replicates) * delta
  } else if(error_type == "additive"){
    measured_samples <- rep(true_samples$z, replicates) + delta
  } else{
    stop("input a valid error_type")
  }
  measured_samples_frame <- samples_frame %>% 
    mutate(measurement = measured_samples) %>%
    mutate(composite_size = rep(true_samples$composite_size, replicates)) %>%
    mutate(strata = rep(true_samples$strata, replicates))
  measured_samples_frame 
}


################ function to bundle measured samples taken from a number of treatment and control plots ################
bundle_samples <- function(samples, treatment_indicator, time = NULL){
  #inputs: 
    #measured samples: a list of data frames. Each entry is samples as output by perturb_measurements()
    #treatment_indicator: a binary vector of length length(measured_samples), each entry indicates whether the plot is in the treatment (1) or control (0) group
    #time: an optional argument specifying the time of measurement of each sample. If there is no time it's simply excluded. Comparing differences is likely to be considerably more powerful then a cross-sectional analysis, (e.g. not measuring at baseline)
  #outputs: 
    #a dataframe of samples with columns for composited sample number, number of samples in composite, plot number, treatment indicator, and time (if given)
  num_replicates <- length(samples[treatment_indicator == 1])
  replicates_frame <- data.frame(treatment = treatment_indicator, plot = 1:length(samples))
  
  sample_frame <- samples %>%
    reduce(bind_rows) %>% 
    mutate(plot = rep(replicates_frame$plot, each = nrow(samples[[1]]))) %>%
    mutate(treatment = rep(replicates_frame$treatment, each = nrow(samples[[1]]))) 
  sample_frame
}




############# estimation: mean, variance, confidence interval for individual plot total ###############
estimate_plot_carbon <- function(plot_samples, confidence_level = .95){
  #goal is to estimate the mean, variance, and a confidence interval for the plot (population) average given any kind of sampling scheme
  #inputs: 
    #plot_samples: a dataframe with samples as output by measure_samples()
    #confidence_level: a desired confidence level to construct an approximate confidence interval using CLT
  #outputs:
    #mean, variance, and confidence interval for plot population mean SOC concentration
  alpha <- 1-confidence_level
  plot_mean <- plot_samples %>%
    group_by(sample) %>% 
    summarize(strata = first(strata), measurement_mean = mean(measurement)) %>%
    group_by(strata) %>%
    summarize(strata_mean = mean(measurement_mean), strata_variance = var(measurement_mean)) %>%
    summarize(plot_mean = mean(strata_mean), plot_variance = 1/n()^2 * sum(strata_variance)) %>%
    mutate(lower_CI = plot_mean + qnorm(alpha/2)*sqrt(plot_variance), upper_CI = plot_mean + qnorm(1-alpha/2)*sqrt(plot_variance))
  plot_mean
}



########### run stratification ANOVA #############
#run an analysis of a stratification (how succesful was the stratification?)
#Reference: page 81 of Lark and Webster "Field sampling for environmental science..."
analyze_stratification <- function(plot_samples){
  #report the intraclass correlation and relative variance of a stratified samples
  #a succesful stratification has intraclass correlation near 1 and relative variance near 0
  #inputs: 
    #plot_samples: a dataframe of samples as output by measure_samples()
  #outputs:
    #a dataframe with one row, columns are intraclass correlation and relative variance
  
  
  #as a first step, take mean across measurement replicates
  plot_samples <- plot_samples %>%
    group_by(sample) %>%
    summarize(strata = first(strata), measurement_mean = mean(measurement))
  
  within_strata_variance <- plot_samples %>%
    group_by(strata) %>%
    summarize(strata_sum_squares = sum((measurement_mean - mean(measurement_mean))^2), strata_DOF = n() - 1) %>%
    summarize(within_strata_variance = sum(strata_sum_squares) / sum(strata_DOF)) %>%
    pull(within_strata_variance)
  
  total_variance <- plot_samples %>%
    summarize(total_variance = var(measurement_mean)) %>%
    pull(total_variance)
  
  B <- plot_samples %>%
    group_by(strata) %>%
    summarize(strata_mean = mean(measurement_mean), strata_size = n()) %>%
    summarize(B = (1 / (n() - 1)) * sum(strata_size * (strata_mean - mean(strata_mean))^2)) %>%
    pull(B)
  
  n_w_star <- plot_samples %>%
    group_by(strata) %>%
    summarize(strata_size = n()) %>%
    summarize(K = n(), N = sum(strata_size), sum_size_squared = sum(strata_size^2)) %>%
    mutate(n_w_star = (1/(K-1)) * (N - sum_size_squared / N)) %>%
    pull(n_w_star)
  
  between_strata_variance <- (B - within_strata_variance) / n_w_star
  
  intraclass_correlation <- between_strata_variance / (within_strata_variance + between_strata_variance)
  
  relative_variance <- within_strata_variance / total_variance
  
  data.frame(intraclass_correlation = intraclass_correlation, relative_variance = relative_variance)
}


########### estimate treatment effect and confidence interval ###########
run_ANOVA <- function(samples, include_plot_variance = TRUE, detail = FALSE){
  #analyze an experiment at a given time point using ANOVA (normal theory inference)
  #input:
    #samples: a dataframe of samples as output by bundle_samples()
    #include_plot_variance: a boolean. If TRUE, return full ANOVA accounting for multiple samples within plots. If FALSE, conduct inference on plot means
    #detail: if FALSE, return only the F-statistic (e.g. for simulations). If TRUE, return full linear model output (e.g. for data analysis)
  #output:
    #either a full linear model (lm class) or the F-statistic from a linear model
  measurement_avg_samples <- samples %>%
    group_by(sample, plot, treatment) %>% 
    summarize(measurement_mean = mean(measurement)) %>%
    ungroup() %>%
    mutate(plot = as_factor(plot), treatment = as_factor(treatment))
  
  if(include_plot_variance){
    linear_model <- aov(measurement_mean ~ treatment + Error(plot), data = measurement_avg_samples)
    p_value <- summary(linear_model)[["Error: plot"]][[1]][["Pr(>F)"]][1] 
  } else {
    measurement_plot_avg_samples <- measurement_avg_samples %>%
      group_by(plot, treatment) %>%
      summarize(measurement_plot_mean = mean(measurement_mean))
    linear_model <- lm(measurement_plot_mean ~ treatment, data = measurement_plot_avg_samples)
    p_value <- anova(linear_model)[["Pr(>F)"]][1]
  }
  if(detail){
    linear_model
  } else {
    p_value
  }    
}

############### randomization inference ###############

shuffle_treatment_assignment <- function(samples){
  #auxiliary function for permutation test: shuffle plots in a dataframe
  #inputs:
    #samples: a dataframe of samples as output by bundle_samples()
  #outputs:
    #the same dataframe with the treatment indicator shuffled
  shuffled_plots_treatment <- samples %>% 
    select(plot, treatment) %>%
    distinct() %>%
    mutate(treatment = sample(treatment, size = length(treatment), replace = FALSE))
  
  shuffled_samples <- samples %>%
    select(-treatment) %>%
    inner_join(shuffled_plots_treatment, by = "plot")
  
  shuffled_samples
}

get_difference_in_means <- function(samples){
  #function for permutation test. Given a dataframe compute difference in means between trt and control
  #inputs:
    #samples: a dataframe of samples as output by bundle_samples()
  #outputs:
    #the difference in means between treatment and control plots (averaged over plots and cores)
  difference_in_means <- samples %>%
    group_by(treatment) %>%
    summarize(trt_mean = mean(measurement_mean)) %>%
    pivot_wider(names_from = treatment, values_from = trt_mean) %>%
    mutate(difference_in_means = `1` - `0`) %>%
    pull(difference_in_means)
}


run_permutation_analysis <- function(samples, B = 1000, bootstrap_cores = FALSE, plot = FALSE){
  #take average over measurement for now
  measurement_avg_samples <- samples %>%
    group_by(sample, plot, treatment, strata) %>% 
    summarize(measurement_mean = mean(measurement)) %>%
    ungroup() %>%
    mutate(plot = as_factor(plot), treatment = as_factor(treatment))
  
  #observed value of test statistic
  difference_in_means <- get_difference_in_means(measurement_avg_samples)
  
  #permutation distribution
  if(!bootstrap_cores){
    permutation_distribution <- replicate(n = B, expr = get_difference_in_means(shuffle_treatment_assignment(measurement_avg_samples)))
  } else {
    permutation_distribution <- replicate(n = B, expr = get_difference_in_means(shuffle_treatment_assignment(bootstrap_cores(measurement_avg_samples))))
  }
  
  if(plot == TRUE){
    hist(permutation_distribution, breaks = 30, xlim = c(min(permutation_distribution, difference_in_means), max(permutation_distribution, difference_in_means)))
    lines(c(difference_in_means, difference_in_means), c(0,B), lwd = 3, col = "red")
  }
  
  #p-value for two sided test
  data.frame(difference_in_means = difference_in_means, p_value = 1 - mean(abs(permutation_distribution) > abs(difference_in_means)))
}


################ optimization of samples over budget ###############
get_variance <- function(n, k, sigma_p, mu, sigma_delta){
  #function to get the variance given parameters. sample size (n), number of measured samples (k), plot heterogeneity (sigma_p), measurement variance (sigma_delta), and average carbon concentration in a plot (mu)
  #input: 
    #n: sample size
    #k: number of measured samples (after compositing)
    #sigma_p: plot heterogeneity (standard deviation of carbon concentration)
    #mu: average concentration of carbon in the plot
    #sigma_delta: measurement variance
  #output:
    #the theoretical variance of the empirical mean of k equally and perfectly composited samples (composited from n total samples) measured with multiplicative error
    
  variance <- (sigma_p^2 * (1 + sigma_delta^2)) / n + (mu^2 * sigma_delta^2) / k
  variance
}



get_minimum_error <- function(sigma_p, sigma_delta, mu, C_0, cost_c, cost_P, cost_A, B, measurement_error = TRUE){
  #solve optimization problem (in closed form, by lagrange multiplier) for simple random sampling and a fixed measurement method that determines sigma_delta
  #input: 
    #sigma_p: the plot variance
    #mu: the average carbon concentration in the plot
    #sigma_delta: the variance of the (multiplicative) measurement error
    #C_0: the fixed cost(s) of the survey
    #cost_c: the cost of collecting a single core (sample) from the plot
    #cost_P: the cost of prepping a single (composited) sample
    #cost_A: the cost of measuring a single (composited) sample from the plot
    #B: the total budget for sampling and measurement
    #measurement_error: option to compute n without measurement error or cost (e.g. to lower bound error as sampling error)
  #output:
    #a dataframe with 4 elements: 
      #n_star: the optimum number of samples to take from the field
      #k_star: the optimum number of samples to measure after compositing
      #total_cost: the total cost of collecting n_star and measuring k_star samples (if not equal to B, we have a problem)
      #optimum_variance: the variance attained when n_star samples are collected and k_star are measured (given parameters sigma_p, sigma_delta, and mu)
  
  #if there is no measurement error and cost of measurement, the solution is very simple
  if(!measurement_error){
    cost_P <- 0
    cost_A <- 0
    sigma_delta <- 0
    n_star <- (B - C_0) / cost_c
    k_star <- n_star #this doesnt actually matter if there is no measurement error or cost of measurement
    
    
  } else {
    #proceed to compute optimization solution
    n_star <- (B - C_0) * sigma_p * sqrt(1 + sigma_delta^2) / (sigma_p * sqrt((1 + sigma_delta^2)) * cost_c + mu * sigma_delta * sqrt((cost_P + cost_A) * cost_c))
    k_star <- (B - C_0) * mu * sigma_delta / ((sigma_p * sqrt((1 + sigma_delta^2) * cost_c * (cost_P + cost_A))  + mu * sigma_delta * (cost_P + cost_A)))
    
    #check boundary conditions and correct if violated
    n_star[which(k_star < 1)] <- ((B - C_0 - cost_P - cost_A) / cost_c)[which(k_star < 1)]
    k_star[which(k_star < 1)] <- 1
    k_star_greater_nstar <- which(k_star > n_star)
    n_star[k_star_greater_nstar] <- ((B - C_0) / (cost_c + cost_P + cost_A))[k_star_greater_nstar]
    k_star[k_star_greater_nstar] <- ((B - C_0) / (cost_c + cost_P + cost_A))[k_star_greater_nstar]
    
    #stop if there's not enough money to do at least one sample and measurement
    max_cost_P <- max(cost_P)
    max_cost_A <- max(cost_A)
    max_cost_c <- max(cost_c)
    if(any(n_star < 1 | k_star < 1)){
      stop(paste("Not enough budget to do at least one sample and measurement! Increase minimum budget to at least", C_0 + max_cost_c + max_cost_A + max_cost_P))
    }
  }
  
  optimum_variance <- get_variance(n = n_star, k = k_star, sigma_p = sigma_p, sigma_delta = sigma_delta, mu = mu)
  
  total_cost <- C_0 + n_star * cost_c + k_star * (cost_P + cost_A)
  
  data.frame(n = n_star, k = k_star, total_cost = total_cost, optimum_variance = optimum_variance)
}

get_composite_error_grid <- function(n, sigma_p, sigma_delta, mu, C_0, cost_c, cost_P, cost_A){
  #for a range of composite sizes 1 to n, return the achieved standard error
  #input: 
    #n: the number of samples
    #sigma_p: the plot variance
    #mu: the average carbon concentration in the plot
    #sigma_delta: the variance of the (multiplicative) measurement error
    #cost_c: the cost of collecting a single core (sample) from the plot
    #cost_P: the cost of prepping a single (composited) sample
    #cost_A: the cost of measuring a single (composited) sample from the plot
  #output:
    #a 2 column dataframe with all tenable composite sizes and achieved standard errors
  possible_k <- 1:n
  k <- possible_k[(n %% possible_k) == 0]
  composite_sizes <- n / k
  variances <- get_variance(n = n, k = k, sigma_p = sigma_p, mu = mu, sigma_delta = sigma_delta)
  costs <- C_0 + n * cost_c + k * (cost_P + cost_A)
  data.frame(composite_size = composite_sizes, std_error = sqrt(variances), cost = costs)
}

get_optimal_composite_size <- function(sigma_p, sigma_delta, mu, cost_c, cost_P, cost_A){
  #get the optimal size of composites, which doesn't depend on the budget or fixed costs
  #input: 
    #sigma_p: the plot variance
    #mu: the average carbon concentration in the plot
    #sigma_delta: the variance of the (multiplicative) measurement error
    #cost_c: the cost of collecting a single core (sample) from the plot
    #cost_P: the cost of prepping a single (composited) sample
    #cost_A: the cost of measuring a single (composited) sample from the plot
  #output:
    #a scalar, the optimal composite size
  optimal_composite_size <- sigma_p * sqrt(1 + sigma_delta^2) * sqrt(cost_P + cost_A) / (mu * sigma_delta * sqrt(cost_c))
  optimal_composite_size
}

get_relative_efficiency <- function(sigma_p, mu, cost_c, sigma_delta_1, cost_P_1, cost_A_1, sigma_delta_2, cost_P_2, cost_A_2){
  #input: 
    #sigma_p: the plot variance
    #mu: the average carbon concentration in the plot
    #cost_c: the cost of collecting a single core (sample) from the plot
    #sigma_delta_1: the measurement error variance of method 1
    #cost_P_1: the sample prep cost of method 1
    #cost_A_1: the measurement cost of method 1
    #sigma_delta_2: the measurement error variance of method 2
    #cost_P_2: the sample prep cost of method 2
    #cost_A_2: the measurement cost of method 2
  #output:
    #a scalar, the relative efficiency of measurement 1 over measurement 2
  numerator <- sigma_p * sqrt((1 + sigma_delta_1) * cost_c) + mu * sigma_delta_1 * sqrt(cost_P_1 + cost_A_1)
  denominator <- sigma_p * sqrt((1 + sigma_delta_2) * cost_c) + mu * sigma_delta_2 * sqrt(cost_P_2 + cost_A_2)
  numerator / denominator
}



get_minimum_cost <- function(sigma_p, sigma_delta, mu, C_0, cost_c, cost_P, cost_A, V, integer_outputs = TRUE){
  #given a fixed precision (variance) that we would like to achieve, what is the minimum total cost of the design input. 
  #input: 
  #sigma_p: the plot variance
  #mu: the average carbon concentration in the plot
  #sigma_delta: the variance of the (multiplicative) measurement error
  #C_0: fixed cost of the survey 
  #cost_c: the cost of collecting a single core (sample) from the plot
  #cost_P: the cost of prepping a single (composited) sample
  #cost_A: the cost of measuring a single (composited) sample from the plot after prep
  #V: the maximum tolerable variance
  #integer_outputs: should outputs be appropriately bounded and rounded to produce integer solutions? If FALSE, outputs are not rounded and bounds (k < n, n,k > 1) are not respected.
  #output:
  #a dataframe with columns for the optimal n and k, the variance attained (which is currently subject to minor floating point errors and so not exactly equal to V), and the minimum cost to attain that variance
  n_star <- (sigma_p^2 * (1+sigma_delta^2)) / (V * (1 - 1/(sqrt(cost_c * sigma_p^2 * (1 + sigma_delta^2) / (cost_P + cost_A)) + 1)))
  k_star <- (sqrt(cost_c * sigma_p^2 * (1 + sigma_delta^2) / (cost_P + cost_A)) + 1) * (mu^2 * sigma_delta^2 / V)
  
  if(integer_outputs){
    #boundary conditions
    #if k is greater than n, then set k = n and solve for n
    k_star_greater_nstar <- which(k_star > n_star)
    n_star[k_star_greater_nstar] <- ((sigma_p^2 * (1 + sigma_delta^2) + mu^2 * sigma_delta^2) / V)[k_star_greater_nstar]
    k_star[k_star_greater_nstar] <- ((sigma_p^2 * (1 + sigma_delta^2) + mu^2 * sigma_delta^2) / V)[k_star_greater_nstar]
    #if k is less than 1 set k = 1 and then solve for n
    n_star[which(k_star < 1)] <- ((sigma_p^2 * (1 + sigma_delta^2)) / (V - mu^2 * sigma_delta^2))[which(k_star < 1)]
    k_star[which(k_star < 1)] <- 1
    #if n is less than 1, set n = 1
    n_star[which(n_star < 1)] <- 1
    
    #take to ceiling to ensure variance constraint is met
    n_star <- ceiling(n_star)
    k_star <- ceiling(k_star)
  }
  
  
  variance <- sigma_p^2 * (1+sigma_delta^2) / n_star + mu^2 * sigma_delta^2 / k_star
  minimum_cost <- C_0 + n_star * cost_c + k_star * (cost_P + cost_A)
  
  
  data.frame(n = n_star, k = k_star, variance = variance, minimum_cost = minimum_cost)
}

############### power analysis #######################
get_power_two_sample <- function(n_1 = NULL, k_1 = NULL, n_2 = NULL, k_2 = NULL, mu_1, mu_2, sigma_p_1, sigma_p_2, sigma_delta, alpha = .05, beta = NULL){
  #return the (asymptotic) power of a two sample t test given sample sizes and plot parameters (no cost model). The null hypothesis is no difference.
  #inputs:
    #n_1: sample size for plot 1
    #k_1: number of assays of samples from plot 1
    #n_2: sample size for plot 2
    #k_2: number of assays of samples from plot 2
    #mu_1: true population mean of plot 1
    #mu_2: true popoulation mean of plot 2 (the "effect size" is mu_1 - mu_2)
    #sigma_p_1: the population heterogeneity (sd) of plot 1
    #sigma_p_2: the population heterogeneity (sd) of plot 2
  #output:
    #the standard deviation of the difference-in-means estimator and the power 
  if(!is.null(n_1) & !is.null(n_2)){
    if(is.null(k_1)){k_1 <- n_1}
    if(is.null(k_2)){k_2 <- n_2}
    var_1 <- get_variance(n = n_1, k = k_1, sigma_p = sigma_p_1, mu = mu_1, sigma_delta = sigma_delta)
    var_2 <- get_variance(n = n_2, k = k_2, sigma_p = sigma_p_2, mu = mu_2, sigma_delta = sigma_delta)
    delta <- mu_1 - mu_2
    std_error_diff <- sqrt(var_1 + var_2)
    dof <- (var_1 + var_2)^2 / (var_1^2 / (n_1 - 1) + var_2^2 / (n_2 - 1))
    critical_value <- qt(p = alpha/2, df = dof, lower.tail = FALSE)
    power <- pt(critical_value, df = dof, ncp = delta/std_error_diff, lower.tail = FALSE) + pt(-critical_value, df = dof, ncp = delta/std_error_diff)
    power
  } else if(!is.null(beta)){
    var_1 <- sigma_p_1^2 * (1 + sigma_delta^2) + mu_1^2 * sigma_delta^2
    var_2 <- sigma_p_2^2 * (1 + sigma_delta^2) + mu_2^2 * sigma_delta^2
    pooled_sd <- sqrt(var_1 + var_2)
    delta <- abs(mu_1 - mu_2)
    z_alpha <- qnorm(1-alpha/2)
    z_beta <- qnorm(1-beta)
    sample_size <- 2 * ((z_alpha + z_beta) / (delta/pooled_sd))^2
    sample_size
  } else {
    stop("Supply either sample sizes (n_1 and n_2) or a targeted type 2 error rate (beta)!")
  }
}
