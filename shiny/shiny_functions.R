library(tidyverse)
library(data.table)
library(gstat)
#library(BalancedSampling)

############### power analysis #######################
get_power_two_sample <- function(n_1 = NULL, k_1 = NULL, n_2 = NULL, k_2 = NULL, mu_1, mu_2, sigma_p_1, sigma_p_2, sigma_delta, alpha = .05, beta = NULL){
  #return the (asymptotic) power of a two sample t test given sample sizes and plot parameters (no cost model). The null hypothesis is no difference.
  #inputs:
  #n_1: sample size for plot 1
  #k_1: number of assays of samples from plot 1
  #n_2: sample size for plot 2
  #k_2: number of assays of samples from plot 2
  #mu_1: true population mean of plot 1
  #mu_2: true population mean of plot 2 (the "effect size" is mu_1 - mu_2)
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


############## minimum cost ###############
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



################ optimization of samples over budget ###############

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




