library(tidyverse)
library(readxl)
library(sampling)

#necessary functions
#stratification helper function to make vector an integer vector with the sum preserved
#input: 
#n_strata: the number of samples to take from each stratum, which could be non-integer, but sums to overall sample size
#output:
#a rounded vector with the sum equal to the sum of the original vector
round_strata_sizes <- function(n_strata){
  rounded_n_strata <- floor(n_strata)
  indices <- tail(order(n_strata-rounded_n_strata), round(sum(n_strata)) - sum(rounded_n_strata))
  rounded_n_strata[indices] <- rounded_n_strata[indices] + 1
  rounded_n_strata
}

#function to get mean and standard error estimate from a stratified sample used alongside 'sampling' package
#input: 
#sample: a dataframe of stratified samples as output by sampling::strata()
#output:
#length 2 vector with the estimate of the population mean and the estimated standard error of that estimate
get_mean_se_stratified <- function(sample, N_strata){
  N <- sum(N_strata)
  strata_weights <- N_strata / N
  n_strata <- as.numeric(table(sample$strata))
  strata_means <- tapply(sample$TC, sample$strata, mean)
  strata_vars <- tapply(sample$TC, sample$strata, var)
  var_estimate <- N^(-2) * sum(N_strata^2 * strata_vars / n_strata)
  c(sum(strata_weights * strata_means), sqrt(var_estimate))
}



#Learned-Miller and Thomas bound https://arxiv.org/pdf/1905.06208.pdf
#Originally proposed by Gaffke 
#inputs:
#x: a vector of iid random samples
#alpha: a desired level for the confidence interval
#B: the number of Monte Carlo iterations to estimate the interval
#side: the side of the confidence interval, either "upper" or "lower"
#outputs:
#an upper or lower confidence bound
gaffke_CI <- function(x, alpha = .05, B = 10000, side = "upper"){
  n <- length(x)
  if(side == "lower"){
    x <- 1 - x
  }
  z <- sort(x, decreasing = FALSE)
  ms <- rep(NA, B)
  s <- c(diff(z), 1 - z[n])
  #u_matrix <- matrix(runif(n * B), nrow = B, ncol = n)
  #u_matrix <- t(apply(u_matrix, 1, sort, decreasing = FALSE))
  #ms <- 1 - u_matrix %*% s
  for(i in 1:B){
    u <- sort(runif(n), decreasing = FALSE)
    ms[i] <- 1 - sum(u * s)
  }
  ms_alpha <- quantile(ms, 1 - alpha)
  if(side == "lower"){
    1 - ms_alpha
  } else if(side == "upper"){
    ms_alpha
  }
}

run_two_sample_t_test <- function(sample_size, pop_1, pop_2){
  sample_1 <- sample(pop_1, size = sample_size, replace = TRUE)
  sample_2 <- sample(pop_2, size = sample_size, replace = TRUE)
  diff_mean <- mean(sample_2) - mean(sample_1)
  std_error <- sqrt(var(sample_1)/sample_size + var(sample_2)/sample_size)
  #this is Welch's t-test
  dof <- std_error^4 / ((var(sample_1)/sample_size)^2 / (sample_size - 1) + (var(sample_2)/sample_size)^2 / (sample_size - 1))
  pt(diff_mean/std_error, df = dof, lower.tail = FALSE)
}


two_sample_martingale <- function(sample_1, sample_2, bounds, alternative = NULL, d = NULL, sequential = FALSE){
  if(length(sample_1) != length(sample_2)){
    stop("Unbalanced sample sizes not yet supported.")
  }
  #assume bounds are equal for the two populations; this makes the difference an RV on [0,1]
  Z <- (sample_2 - sample_1 + bounds[2]) / (2 * diff(bounds)) 
  n <- length(Z)
  #the null is that the means of the lists are equal
  null_mean <- 1/2
  if(is.null(alternative) | is.null(d)){
    eta_j <- cummean(Z) 
  } else{
    scaled_alt_mean <- (alternative + diff(bounds)) / (2 * diff(bounds))
    eta_j <- pmin(pmax((d * scaled_alt_mean + cumsum(Z)) / (d + 1:n - 1), null_mean + ((scaled_alt_mean - null_mean) / 2) / sqrt(d + 1:n - 1)), 1)
  }
  terms <- c(1, (Z/null_mean) * (eta_j - null_mean)/(1 - null_mean) + (1 - eta_j)/(1 - null_mean))
  T_j <- cumprod(terms)
  if(sequential){
    pval <- 1 / T_j
  } else{
    pval <- 1 / max(T_j)
  }
  pval
}

#a function to run a one-sided two-sample nonparametric test based on the Gaffke (LMT) test
#the hypothesis H_0: mu_2 <= mu_1 is tested
#populations are assumed to be bounded between [0,1] but can be pre-processed to accomodate any bounded distribution 
#subtract the lower bound, divide by the upper bound to get to [0,1]
#inputs:
#sample_1: a vector, a random sample from a larger population
#sample_2: a vector, a random sample from a second population. 
#alpha: a double in (0,1), the desired level of the test
#B: a positive integer, the number of Monte Carlo iterations to be used in running the test
#method: a string, the method of combining the two one-sample tests into a two-sample test of equality
#Sidak: compute a 1-sqrt(1-alpha) upper and lower bound (respectively) and see if they overlap
#Fisher: compute a combined p-value using Fisher's combining function of whether both means are equal to a particular mu_0, then maximize this over possible values of mu_0
#Liptak: compute a combined p-value using Liptak's combining function
#pval: a boolean, if TRUE return a pvalue (only works for Fisher or Liptak), if FALSE return a boolean indicating rejection at level alpha 
two_sample_gaffke_test <- function(sample_1, sample_2, alpha, B = 1000, method = "sidak", bounds = c(0,1), pval = FALSE){
  sample_1 <- (sample_1 - bounds[1]) / diff(bounds)
  sample_2 <- (sample_2 - bounds[1]) / diff(bounds)
  if(method == "sidak"){
    upper_1 <- gaffke_CI(x = sample_1, alpha = 1-sqrt(1-alpha), B = B, side = "upper")
    lower_2 <- gaffke_CI(x = sample_2, alpha = 1-sqrt(1-alpha), B = B, side = "lower")
    reject <- ifelse(upper_1 < lower_2, TRUE, FALSE)
    reject
  } else if(method %in% c("fisher","liptak","tippett")){
    n_1 <- length(sample_1)
    n_2 <- length(sample_2)
    x_1 <- sample_1
    x_2 <- sample_2
    ms_1 <- rep(NA, B)
    ms_2 <- rep(NA, B)
    for(b in 1:B){
      Z_1 <- rexp(n_1 + 1)
      D_1 <- Z_1 / sum(Z_1)
      Z_2 <- rexp(n_2 + 1)
      D_2 <- Z_2 / sum(Z_2)
      ms_1[b] <- sum(D_1 * c(x_1, 1))
      ms_2[b] <- sum(D_2 * c(x_2, 0))
    }
    if(method == "fisher"){
      combined_p <- function(mu_0){
        p_1 <- mean(c(ms_1 >= mu_0, TRUE))
        p_2 <- mean(c(ms_2 <= mu_0, TRUE))
        combined_test_stat <- -2 * (log(p_1) + log(p_2))
        p_val <- pchisq(q = combined_test_stat, df = 4, lower.tail = FALSE)
        p_val
      }
    } else if(method == "liptak"){
      combined_p <- function(mu_0){
        p_1 <- mean(c(ms_1 >= mu_0, TRUE))
        p_2 <- mean(c(ms_2 <= mu_0, TRUE))
        combined_test_stat <- qnorm(1 - p_1) + qnorm(1 - p_2)
        p_val <- pnorm(q = combined_test_stat, sd = 2, lower.tail = FALSE)
        p_val
      }
    } else{
      combined_p <- function(mu_0){
        p_1 <- mean(c(ms_1 >= mu_0, TRUE))
        p_2 <- mean(c(ms_2 <= mu_0, TRUE))
        combined_test_stat <- min(p_1, p_2)
        p_val <- pbeta(q = combined_test_stat, shape1 = 1, shape2 = 2, lower.tail = TRUE)
        p_val
      }
    }
    
    max_p_val <- optimize(combined_p, lower = min(c(ms_1, ms_2)), upper = max(c(ms_1, ms_2)), maximum = TRUE)$objective
    if(pval){
      max_p_val
    } else{
      ifelse(max_p_val < alpha, TRUE, FALSE)
    }
  } else{
    stop("Input a valid method: sidak, fisher, liptak, or tippett")
  }
}



#read in data
rangeland_master <- read_excel("R_Heterogeneity_Master_PS_04282021.xlsx", sheet = "Rangeland_All_soliTOC") %>%
  mutate(TOC = ifelse(TOC == "NA", NA, TOC)) %>%
  mutate(TIC = ifelse(TIC == "NA", NA, TIC)) %>%
  mutate(TC = ifelse(TC == "NA", NA, TC)) %>%
  mutate(TOC = as.numeric(TOC), TIC = as.numeric(TIC), TC = as.numeric(TC), sample_number = as.numeric(sample_number)) %>%
  mutate(depth_long = dplyr::recode(depth, a = "0-10 cm", b = "10-30 cm", c = "30-50 cm", d = "50-75 cm", e = "75-100 cm"))
cropland_master <- read_excel("R_Heterogeneity_Master_PS_04282021.xlsx", sheet = "Cropland_All_Costech") %>%
  mutate(depth_long = dplyr::recode(depth, a = "0-15 cm", b = "15-30 cm", c = "30-60 cm", d = "60-100 cm"))



# power of t-test (unstratified and stratified) or nonparametric (unstratified) to detect topsoil change 
# run tests on topsoil from rangeland and cropland
topsoil_rangeland <- rangeland_master %>% 
  filter(depth == "a") %>%
  filter(!is.na(TC)) %>%
  arrange(transect) %>%
  select(transect, TC)
topsoil_cropland <- cropland_master %>% 
  filter(depth == "a") %>%
  filter(site == "CROP5") %>%
  filter(!is.na(TC)) 

topsoil_TC_rangeland <- topsoil_rangeland %>% pull(TC)
topsoil_TC_cropland <- topsoil_cropland %>% pull(TC)


# par(mfrow = c(2,1))
# hist(topsoil_TC_rangeland - mean(topsoil_TC_rangeland) + 3, breaks = 30, xlab = "Percent Carbon", main = "", cex.lab = 1.5, cex.axis = 2, xlim = c(0,10))
# hist(-(topsoil_TC_cropland - mean(topsoil_TC_cropland)) + 3, breaks = 30, xlab = "Percent Carbon", main = "", cex.lab = 1.5, cex.axis = 2, xlim = c(0,10))
# par(mfrow = c(1,1))


###### VALIDITY SIMULATIONS #####

run_validity_simulations <- function(type, n_sims = 500){
  N <- 1000
  
  if(type == "skewed"){
    #means are exactly 3 in both populations. Population 1 is highly skewed so that high values (which make means equal) are rarely sampled
    pop_1 <- c(rnorm(N-10, mean = 3, sd = .05), rep(20, 10))
    pop_1 <- pop_1 - mean(pop_1) + 3
    pop_2 <- rep(3, N) + rnorm(N, sd = .05)
    pop_2 <- pop_2 - mean(pop_2) + 3
  } else if(type == "rangeland_to_cropland"){
    pop_1 <- topsoil_TC_rangeland - mean(topsoil_TC_rangeland) + 3
    pop_2 <- topsoil_TC_cropland - mean(topsoil_TC_cropland) + 3
  } else if(type == "rangeland_to_negcropland"){
    pop_1 <-  topsoil_TC_rangeland - mean(topsoil_TC_rangeland) + 3
    pop_2 <- mean(topsoil_TC_cropland) - topsoil_TC_cropland + 3
  } else if(type == "symmetric_reduced_spread"){
    pop_1 <- rnorm(N, mean = 3, sd = .5)
    pop_2 <- rnorm(N, mean = 3, sd = .1)
    pop_1 <- pop_1 - mean(pop_1) + 3
    pop_2 <- pop_2 - mean(pop_2) + 3
  } else if(type == "gaussian_gaussian"){
    pop_1 <- rnorm(N, mean = 3, sd = .5)
    pop_2 <- rnorm(N, mean = 3, sd = .5)
    pop_1 <- pop_1 - mean(pop_1) + 3
    pop_2 <- pop_2 - mean(pop_2) + 3
  } else if(type == "rangeland_to_gaussian"){
    pop_1 <-  topsoil_TC_rangeland - mean(topsoil_TC_rangeland) + 3
    pop_2 <- rnorm(N, mean = 3, sd = .5)
    pop_2 <- pop_2 - mean(pop_2) + 3
  } else if(type == "cropland_to_gaussian"){
    pop_1 <-  topsoil_TC_cropland - mean(topsoil_TC_cropland) + 3
    pop_2 <- rnorm(N, mean = 3, sd = .5)
    pop_2 <- pop_2 - mean(pop_2) + 3
  } 
  else if(type == "rangeland_to_negrangeland"){
    pop_1 <-  topsoil_TC_rangeland - mean(topsoil_TC_rangeland) + 6
    pop_2 <- mean(topsoil_TC_rangeland) - topsoil_TC_rangeland + 6
  } else if(type == "cropland_to_negcropland"){
    pop_1 <-  topsoil_TC_rangeland - mean(topsoil_TC_rangeland) + 6
    pop_2 <- mean(topsoil_TC_rangeland) - topsoil_TC_rangeland + 6
  } else{
    stop("Supply valid argument to type")
  }
  
  n_grid <- seq(6,200, by = 2)
  
  t_test_rejection_rate <- rep(0, length(n_grid))
  gaffke_rejection_rate <- rep(0, length(n_grid))
  
  for(i in 1:length(n_grid)){
    t_test_p_values <- replicate(n = n_sims, run_two_sample_t_test(sample_size = n_grid[i], pop_1, pop_2))
    gaffke_reject <- replicate(n = n_sims, two_sample_gaffke_test(sample_1 = sample(pop_1/20, size = n_grid[i], replace = T), sample_2 = sample(pop_2/20, size = n_grid[i], replace = T), alpha = .05, B = 100, method = "fisher"))
    t_test_rejection_rate[i] <- mean(t_test_p_values < .05)
    gaffke_rejection_rate[i] <- mean(gaffke_reject)
  }
  data.frame("sample_size" = n_grid, "t_test_rejections" = t_test_rejection_rate, "gaffke_rejections" = gaffke_rejection_rate, "population" = type)
}


# pop_list <- c("skewed", "rangeland_to_cropland", "rangeland_to_negcropland", "rangeland_to_negrangeland", "rangeland_to_gaussian", "cropland_to_negcropland", "cropland_to_gaussian", "gaussian_gaussian")
# validity_simulations <- lapply(pop_list, run_validity_simulations, n_sims = 5000)
# save(validity_simulations, file = "validity_simulations")

# par(mar = c(5.1,4.3, 4.3, 2.1))
# plot(y = t_test_rejection_rate, x = n_grid, type ='l', ylim = c(0,1), xlab = "Sample size", ylab = "Simulated significance level", col = 'darkorange3', lwd = 4, cex.axis = 1.8, cex.lab = 1.8)
# points(y = gaffke_rejection_rate, x = n_grid, type = 'l', col = 'steelblue', lwd = 4)
# #points(y = hedged_rejection_rate, x = n_grid, type = 'l', col = 'darkorange3', lwd = 2, lty = "dashed")
# legend(x = 100, y = 0.8, legend = c("Nonparametric test","t-test"), lty = c( "solid","solid"), col = c("steelblue","darkorange3"), lwd = 4, bty = "n", cex = 1.5)
# abline(a = 0.05, b = 0, lty = 'dashed', col = 'black', lwd = 2)
# par(mar = c(5.1,4.1, 4.1, 2.1))






###### POWER SIMULATIONS ######
run_twosample_sims <- function(land_use, sample_size, n_sims = 300, effect = "shift"){
  if(land_use == "cropland"){
    x <- topsoil_TC_cropland
  } else if(land_use == "rangeland"){
    x <- topsoil_TC_rangeland
  } else{
    stop("input valid land_use")
  }
  effect_grid <- seq(0,.6,by=.025)
  mu <- mean(x)
  n_x <- length(x)
  shift <- effect_grid * mu
  scaling <- 1+effect_grid
  spike <- cbind(matrix(0, ncol = floor(n_x * .8), nrow = length(effect_grid)), 5*matrix(rep(shift, each = ceiling(n_x * .2)), ncol = ceiling(n_x * .2), byrow = TRUE))
  
  t_test_p_values <- matrix(NA, nrow = n_sims, ncol = length(effect_grid))
  stratified_t_test_p_values <- matrix(NA, nrow = n_sims, ncol = length(effect_grid))
  gaffke_rejections_10 <- matrix(NA, nrow = n_sims, ncol = length(effect_grid))
  gaffke_rejections_20 <- matrix(NA, nrow = n_sims, ncol = length(effect_grid))
  
  
  for(i in 1:n_sims){
    for(j in 1:length(shift)){
      sample_1 <- sample(x, size = sample_size, replace = TRUE)
      if(effect == "shift"){
        pop_2 <- x + shift[j]
      }
      if(effect == "scaling"){
        pop_2 <- x * scaling[j]
      }
      if(effect == "spike"){
        pop_2 <- x + spike[j,]
      }
      if(effect == "to_gaussian"){
        pop_2 <- rnorm(n = length(x), mean = 0, sd = .25/diff(bounds)) 
        pop_2 <- pop_2 - mean(pop_2) + mu + shift[j]
      }
      sample_2 <- sample(pop_2, size = sample_size, replace = TRUE)
      diff_mean <- mean(sample_2) - mean(sample_1)
      std_error <- sqrt(var(sample_1)/sample_size + var(sample_2)/sample_size)
      #this is for Welch's t-test
      dof <- std_error^4 / ((var(sample_1)/sample_size)^2 / (sample_size - 1) + (var(sample_2)/sample_size)^2 / (sample_size - 1))
      t_test_p_values[i,j] <- pt(diff_mean/std_error, df = dof, lower.tail = FALSE)
      gaffke_rejections_10[i,j] <- two_sample_gaffke_test(sample_1 = sample_1, sample_2 = sample_2, B = 200, bounds = c(0,10), pval = TRUE, method = "fisher")
      gaffke_rejections_20[i,j] <- two_sample_gaffke_test(sample_1 = sample_1, sample_2 = sample_2, B = 200, bounds = c(0,20), pval = TRUE, method = "fisher")
      #mart_p_values[i,j] <- two_sample_martingale(sample_1 = sample_1, sample_2 = sample_2, bounds = bounds, d = 5, alternative = shift[j])
      if(land_use == "rangeland"){
        strata <- topsoil_rangeland$transect
        N_strata <- as.numeric(table(topsoil_rangeland$transect))
        strata_weights_prop <- N_strata / length(strata)
        sigma_strata <- tapply(topsoil_TC_rangeland, strata, sd)
        strata_weights_opt <- N_strata * sigma_strata / sum(N_strata * sigma_strata) 
        
        #dataframe format needed to work with sampling package
        K <- length(unique(strata))
        pop_1_frame <- data.frame(TC = x, strata = strata)
        #n_strata_opt <- round_strata_sizes(sample_size * strata_weights_opt)
        n_strata_prop <- round_strata_sizes(sample_size * strata_weights_prop)
    
        pop_2_frame <- data.frame(TC = pop_2, strata = strata)
        strat_sample_1 <- strata(data = pop_1_frame, stratanames = "strata", size = n_strata_prop, method = "srswr")
        strat_sample_2 <- strata(data = pop_2_frame, stratanames = "strata", size = n_strata_prop, method = "srswr")
        stratified_estimate_1 <- get_mean_se_stratified(sample = getdata(pop_1_frame, strat_sample_1), N_strata = table(pop_1_frame$strata))
        stratified_estimate_2 <- get_mean_se_stratified(sample = getdata(pop_2_frame, strat_sample_2), N_strata = table(pop_2_frame$strata))
        difference_estimate <- stratified_estimate_2[1] - stratified_estimate_1[1]
        combined_se <- sqrt(stratified_estimate_1[2]^2 + stratified_estimate_2[2]^2)
        stratified_t_test_p_values[i,j] <- 1 - pt(q = difference_estimate / combined_se, df = 2 * (sample_size - K))
      }
    }
  }
  t_test_power <- colMeans(t_test_p_values < .05)
  stratified_t_test_power <- colMeans(stratified_t_test_p_values < .05)
  #martingale_power <- colMeans(mart_p_values < .1)
  gaffke_power_bound10_alpha05 <- colMeans(gaffke_rejections_10 < .05)
  gaffke_power_bound10_alpha10 <- colMeans(gaffke_rejections_10 < .10)
  gaffke_power_bound20_alpha05 <- colMeans(gaffke_rejections_20 < .05)
  gaffke_power_bound20_alpha10 <- colMeans(gaffke_rejections_20 < .10)
  data.frame("t_test" = t_test_power, "gaffke_bound10_alpha05"= gaffke_power_bound10_alpha05, "gaffke_bound10_alpha10" = gaffke_power_bound10_alpha10, "gaffke_bound20_alpha10" = gaffke_power_bound20_alpha10, "gaffke_bound20_alpha05" = gaffke_power_bound20_alpha05, "stratified_t_test" = stratified_t_test_power, "effect_size" = effect_grid, "land_use" = land_use, "sample_size" = sample_size, effect = effect)
}

alternative_list <- c("shift", "scaling")

power_simulations_cropland_10 <- lapply(alternative_list, run_twosample_sims, land_use = "cropland", sample_size = 10, n_sims = 500)
power_simulations_cropland_30 <- lapply(alternative_list, run_twosample_sims, land_use = "cropland", sample_size = 30, n_sims = 500)
power_simulations_cropland_90 <- lapply(alternative_list, run_twosample_sims, land_use = "cropland", sample_size = 90, n_sims = 500)
power_simulations_cropland_200 <- lapply(alternative_list, run_twosample_sims, land_use = "cropland", sample_size = 200, n_sims = 500)
power_simulations_rangeland_10 <- lapply(alternative_list, run_twosample_sims, land_use = "rangeland", sample_size = 10, n_sims = 500)
power_simulations_rangeland_30 <- lapply(alternative_list, run_twosample_sims, land_use = "rangeland", sample_size = 30, n_sims = 500)
power_simulations_rangeland_90 <- lapply(alternative_list, run_twosample_sims, land_use = "rangeland", sample_size = 90, n_sims = 500)
power_simulations_rangeland_200 <- lapply(alternative_list, run_twosample_sims, land_use = "rangeland", sample_size = 200, n_sims = 500)

power_simulations <- list(
  power_simulations_cropland_10,
  power_simulations_cropland_30,
  power_simulations_cropland_90,
  power_simulations_cropland_200,
  power_simulations_rangeland_10,
  power_simulations_rangeland_30,
  power_simulations_rangeland_90,
  power_simulations_rangeland_200
)
save(power_simulations, file = "power_simulations")
