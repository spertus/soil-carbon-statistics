library(tidyverse)
library(readxl)
library(permuter)
library(sampling)
library(geoR)
#permuter is used for permutation tests. It can be installed from Github: https://github.com/statlab/permuter
source("functions.R")

########### read in data from excel spreadsheet(s) ##############
#data from Paicines Ranch
#from the original data sent to me by Paige Stanley, I deleted a few free-floating cells on the Rangeland_BD sheet. I also changed 'dulk_density_g/cm3' to 'bulk_density_g/cm3' on Rangeland_BD and Cropland_BD.

paicines_master <- read_excel("../Data/R_Heterogeneity_Master_PS_03112021.xlsx", sheet = "Rangeland_All_soliTOC") %>%
  rename(TOC = TOC) %>%
  mutate(TOC = ifelse(TOC == "NA", NA, TOC)) %>%
  mutate(TIC = ifelse(TIC == "NA", NA, TIC)) %>%
  mutate(TC = ifelse(TC == "NA", NA, TC)) %>%
  mutate(TOC = as.numeric(TOC), TIC = as.numeric(TIC), TC = as.numeric(TC), sample_number = as.numeric(sample_number)) %>%
  mutate(depth_long = dplyr::recode(depth, a = "0-10 cm", b = "10-30 cm", c = "30-50 cm", d = "50-75 cm", e = "75-100 cm"))

paicines_solitoc_reps <- read_excel("../Data/R_Heterogeneity_Master_PS_03112021.xlsx", sheet = "Rangeland_Reps_soliTOC") %>%
  mutate(TOC = ifelse(TOC == "NA", NA, TOC)) %>%
  mutate(TOC = as.numeric(TOC), sample_number = as.numeric(sample_number)) %>%
  mutate(machine = "solitoc")

paicines_costech_reps <- read_excel("../Data/R_Heterogeneity_Master_PS_03112021.xlsx", sheet = "Rangeland_Reps_Costech") %>%
  mutate(TC = ifelse(TC == "NA", NA, TC), TOC = ifelse(TOC == "NA", NA, TOC)) %>%
  mutate(TC = as.numeric(TC), TOC = as.numeric(TOC), sample_number = as.numeric(sample_number)) %>%
  mutate(machine = "costech")

paicines_BD <- read_excel("../Data/R_Heterogeneity_Master_PS_03112021.xlsx", sheet = "Rangeland_BD") %>%
  rename('bd' = 'bulk_density_g/cm3') %>%
  mutate(bd = as.numeric(bd)) %>%
  mutate(depth_long = dplyr::recode(depth, a = "0-10 cm", b = "10-30 cm", c = "30-50 cm", d = "50-75 cm", e = "75-100 cm"))

#data from various croplands around California
cropland_master <- read_excel("../Data/R_Heterogeneity_Master_PS_03112021.xlsx", sheet = "Cropland_All_Costech") 

cropland_solitoc_reps <- read_excel("../Data/R_Heterogeneity_Master_PS_03112021.xlsx", sheet = "Cropland_Reps_soliTOC") %>%
  mutate(sample_number = as.numeric(sample_number)) %>%
  mutate(machine = "solitoc")

cropland_costech_reps <- read_excel("../Data/R_Heterogeneity_Master_PS_03112021.xlsx", sheet = "Cropland_Reps_Costech") %>%
  mutate(TC = ifelse(TC == "NA", NA, TC), TOC = ifelse(TOC == "NA", NA, TOC)) %>%
  mutate(TC = as.numeric(TC), TOC = as.numeric(TOC), sample_number = as.numeric(sample_number)) %>%
  mutate(machine = "costech")


cropland_BD <- read_excel("../Data/R_Heterogeneity_Master_PS_03112021.xlsx", sheet = "Cropland_BD")  %>%
  rename('bd' = 'bulk_density_g/cm3') %>%
  mutate(bd = as.numeric(bd))

replicates_comparison <- read_excel("../Data/R_Heterogeneity_Master_PS_03112021.xlsx", sheet = "Replicate_Comparison") %>% 
  rename(TC_solitoc = TC_soliTOC, TOC_solitoc = TOC_soliTOC) %>%
  mutate_at(vars(starts_with(c("TC_", "TOC_"))), as.numeric)


###################### TOC concentration and BD in space ######################

#stacked histograms for Paicines data
ggplot(paicines_master, aes(TC)) +
  geom_histogram() +
  facet_grid(depth_long ~ transect) +
  xlab("% Carbon") +
  ylab("Number of Samples") +
  theme_bw() +
  theme(text = element_text(size = 16))


#stacked histograms for cropland data
ggplot(cropland_master, aes(TC)) +
  geom_histogram() +
  facet_grid(depth ~ site) +
  xlab("% Carbon") +
  ylab("Number of Samples") +
  theme_bw() +
  theme(text = element_text(size = 16))


#summary tables
paicines_summary <- paicines_master %>% 
  group_by(depth, transect) %>%
  summarize(mean_TC = mean(TC, na.rm = T), sd_TC = sd(TC, na.rm = T))

cropland_summary <- cropland_master %>% 
  group_by(depth, site) %>%
  summarize(mean_TC = mean(TC, na.rm = T), sd_TC = sd(TC, na.rm = T))


#linear models
paicines_model <- lm(TC ~ transect*depth, data = paicines_master)
cropland_model <- lm(TC ~ site*depth, data = cropland_master)

anova(paicines_model)
anova(cropland_model)


#bulk densities
#NOTE: histograms are not disaggregated by location, which drives variation
ggplot(paicines_BD, aes(bd)) +
  geom_histogram() +
  facet_grid(depth_long ~ .) +
  xlab("Bulk Density") +
  ylab("Number of samples") +
  xlim(.8, 2.5) +
  theme_bw() +
  theme(text = element_text(size = 16))


ggplot(cropland_BD, aes(bd)) +
  geom_histogram() +
  facet_grid(depth ~ .) +
  xlab("Bulk Density") +
  ylab("Number of Samples") +
  xlim(.8, 2.5) +
  theme_bw() +
  theme(text = element_text(size = 16))

paicines_summary_bd <- paicines_BD %>%
  group_by(depth, transect) %>%
  summarize(mean_bd = mean(bd, na.rm = TRUE), sd_bd = sd(bd, na.rm = TRUE))

cropland_summary_bd <- cropland_BD %>%
  group_by(depth, site) %>%
  summarize(mean_bd = mean(bd, na.rm = TRUE), sd_bd = sd(bd, na.rm = TRUE))



############## replicates and assay error ##############
assay_error <- replicates_comparison %>%
  group_by(site, sample_number) %>%
  summarize(
    sigma_delta_TC_solitoc = sqrt(var(TC_solitoc, na.rm = T) / (mean(TC_solitoc,  na.rm = T)^2 - var(TC_solitoc, na.rm = T)/n())),
    TC_solitoc = mean(TC_solitoc, na.rm = T),
    sigma_delta_TOC_solitoc = sqrt(var(TOC_solitoc, na.rm = T) / (mean(TOC_solitoc,  na.rm = T)^2 - var(TOC_solitoc, na.rm = T)/n())),
    TOC_solitoc = mean(TOC_solitoc, na.rm = T),
    sigma_delta_TC_costech = sqrt(var(TC_costech, na.rm = T) / (mean(TC_costech,  na.rm = T)^2 - var(TC_costech, na.rm = T)/n())),
    TC_costech = mean(TC_costech, na.rm = T),
    sigma_delta_TOC_costech = sqrt(var(TOC_costech, na.rm = T) / (mean(TOC_costech,  na.rm = T)^2 - var(TOC_costech, na.rm = T)/n())),
    TOC_costech = mean(TOC_costech, na.rm = T),
  ) %>%
  ungroup() 


assay_error_long <- assay_error %>%
  select(site, sample_number, sigma_delta_TC_solitoc, sigma_delta_TC_costech) %>%
  pivot_longer(cols = c("sigma_delta_TC_solitoc", "sigma_delta_TC_costech"), names_to = "machine", names_prefix = "sigma_delta_TC_", values_to = "sigma_delta_TC")
  

#compare average errors on paicines samples
ggplot(assay_error_long, aes(sigma_delta_TC*100, fill = machine)) +
  geom_density(alpha = .5) +
  xlim(0,10) +
  xlab("Percent Error")

median_sigma_delta <- assay_error_long %>%
  group_by(machine) %>%
  summarize(sigma_delta = median(sigma_delta_TC))



#non-parametric permutation test for differences in measurement

#nonparametric analysis of no difference in labs/machines:
reps_long_TC <- replicates_comparison %>% 
  select(site, sample_number, TC_solitoc, TC_costech) %>%
  pivot_longer(cols = c("TC_solitoc", "TC_costech"), names_prefix = "TC_", names_to = "machine", values_to = "TC") %>%
  na.omit()

#there are only a few measurements for TOC on the costech
reps_long_TOC <- replicates_comparison %>% 
  select(site, sample_number, TOC_solitoc, TOC_costech) %>%
  group_by(sample_number) %>%
  filter(sum(!is.na(TOC_costech)) >= 2) %>%
  ungroup() %>%
  pivot_longer(cols = c("TOC_solitoc", "TOC_costech"), names_prefix = "TOC_", names_to = "machine", values_to = "TOC") %>%
  na.omit()


B <- 10000
strata <- as_factor(reps_long_TOC$sample_number)
num_strata <- length(unique(strata))
#test statistic is difference in means between labs
diff_means <- rep(0, num_strata)
pooled_ses <- rep(0, num_strata)
null_distributions <- matrix(0, ncol = num_strata, nrow = B)
p_values <- rep(0, num_strata)
#iterate across strata (samples)
for(k in 1:num_strata){
  solitoc_TOC <- reps_long_TOC$TOC[reps_long_TOC$machine == "solitoc" & strata == levels(strata)[k]]
  costech_TOC <- reps_long_TOC$TOC[reps_long_TOC$machine == "costech" & strata == levels(strata)[k]]
  diff_means[k] <- mean(solitoc_TOC) - mean(costech_TOC)
  pooled_ses[k] <- sqrt(var(solitoc_TOC) / length(solitoc_TOC) + var(costech_TOC) / length(solitoc_TOC))
  null_distributions[,k] <- two_sample(x = solitoc_TOC, y = costech_TOC, reps = B)
  p_values[k] <- mean(abs(diff_means[k]) <= abs(null_distributions[,k])) 
}

npc_p_value <- npc(diff_means, distr = null_distributions, combine=  "fisher", alternatives = "two-sided")



############# spatial correlation of TOC concentrations at Paicines ##########
#first just do for a single transect
T_topsoil_TOC <- paicines_master %>% 
  filter(transect == "T", depth == "a") %>%
  arrange(sample_number) %>%
  pull(TOC)

#dataframe for use with geodata
T_topsoil_geodata <- as.geodata(data.frame(x = 1:length(T_topsoil_TOC), y = 1, TOC = T_topsoil_TOC), coords.col = 1:2, data.col = 3)
T_topsoil_variogram <- variog(T_topsoil_geodata, option = "bin")
plot(T_topsoil_variogram)

#function to plot a variogram for any depth and transect
#inputs:
  #paicines_depth: a string denoting the depth, either "a", "b", "c", "d", or "e"
  #paicines_transect: a string denoting the transect, eitehr "T", "Mx", "My", "Bx", or "By"
#output:
  #a plot of the empirical variogram generated by binning variances between points at a range of distances
plot_variogram <- function(paicines_depth = "a", paicines_transect = "T"){
  if(!(paicines_depth %in% c("a","b","c","d","e")) | !(paicines_transect %in% c("T","Mx","My","Bx","By"))){
    stop("Depth or transect is invalid! See comments.")
  }
  TOC <- paicines_master %>% 
    filter(transect == paicines_transect, depth == paicines_depth) %>%
    arrange(sample_number) %>%
    pull(TOC)
  
  #dataframe for use with geodata
  geodataframe <- as.geodata(data.frame(x = 3 * 1:length(TOC), y = 1, TOC = TOC), coords.col = 1:2, data.col = 3)
  variogram <- variog(geodataframe, option = "bin")
  plot(variogram, xlab = "Distance (m)", ylab = "Semivariance", main = paste("Variogram for transect", paicines_transect, "and depth", paicines_depth))
}

par(mfrow = c(5,1))
plot_variogram(paicines_depth = "a", paicines_transect = "T")
plot_variogram(paicines_depth = "a", paicines_transect = "Mx")
plot_variogram(paicines_depth = "a", paicines_transect = "My")
plot_variogram(paicines_depth = "a", paicines_transect = "Bx")
plot_variogram(paicines_depth = "a", paicines_transect = "By")
par(mfrow = c(1,1))





############### costs and optimal compositing #############
#Costs based on 'Cost of Carbon Analysis.xlsx' assembled by Jessica Chiartas
#in sampling, about 20-24 samples (4 transects * 5-6 composites) costs 400$ at 25$/hr and 16hrs, so we'll take 20 USD per sample.
#the commercial labs had costs of 19.5, 7.75, 22, 10, and 8.75 USD per sample to assess TOC. So we'll take the cost of measurement with dry combustion as 13.60, including sample prep.
#the in-house cost of measurement was 2.78


#sample size
sample_size <- 60
#max budget, i.e. the budget to sample and measure sample_size samples, with no compositing
max_budget <- 20 * sample_size + 13.6 * sample_size
#sample size if we measure once and devote the rest of the budget to sampling
max_sample_size <- floor((max_budget - 13.6) / 20)
#measurement error of a costech
sigma_delta_costech <- median_sigma_delta %>% 
  filter(machine == "costech") %>%
  pull(sigma_delta)

#precision of the sample mean 
precision_paicines <- paicines_master %>%
  group_by(depth) %>%
  summarize(mu = mean(TC, na.rm = T), sigma_p = sd(TC, na.rm = T)) %>%
  mutate(se_no_compositing = sqrt(get_variance(n = sample_size, k = sample_size, sigma_p = sigma_p, mu = mu, sigma_delta = sigma_delta_costech))) %>%
  mutate(se_full_compositing = sqrt(get_variance(n = max_sample_size, k = 1, sigma_p = sigma_p, mu = mu, sigma_delta = sigma_delta_costech))) %>%
  mutate(se_optimal_compositing = sqrt(get_minimum_error(sigma_p = sigma_p, mu = mu, sigma_delta = sigma_delta_costech, C_0 = 0, cost_c = 20, cost_P = 0, cost_A = 13.60, B = max_budget)$optimum_variance)) %>%
  mutate(composite_efficiency_ratio = se_optimal_compositing / se_no_compositing) %>%
  mutate(optimal_composite_size_commercial = get_optimal_composite_size(sigma_p = sigma_p, mu = mu, sigma_delta = sigma_delta_costech, cost_c = 20, cost_P = 0, cost_A = 13.60)) %>%
  mutate(optimal_composite_size_inhouse = get_optimal_composite_size(sigma_p = sigma_p, mu = mu, sigma_delta = sigma_delta_costech, cost_c = 20, cost_P = 0, cost_A = 2.78))

#power of two-sample t-test to detect a half-percentage point change in TOC 
power_change_paicines <- paicines_master %>%
  group_by(depth) %>%
  summarize(mu = mean(TC, na.rm = T), sigma_p = sd(TC, na.rm = T)) %>%
  mutate(power_no_compositing = get_power_two_sample(n_1 = sample_size, k_1 = sample_size, n_2 = sample_size, k_2 = sample_size, mu_1 = mu, sigma_p_1 = sigma_p, mu_2 = mu + 0.5, sigma_p_2 = sigma_p, sigma_delta = sigma_delta_costech)) %>%
  mutate(power_full_compositing = get_power_two_sample(n_1 = max_sample_size, k_1 = 1, n_2 = max_sample_size, k_2 = 1, mu_1 = mu, sigma_p_1 = sigma_p, mu_2 = mu + 0.5, sigma_p_2 = sigma_p, sigma_delta = sigma_delta_costech))
  

#power of a permutation test to detect topsoil change
topsoil_TOC_paicines <- paicines_master %>% 
  filter(depth == "a" & !is.na(TOC)) %>%
  pull(TOC) 
  
#number of simulations to run
n_sims <- 200
#number of samples to take from each plot
sample_size <- 60
#fixed effect
shift <- seq(0, 2, by=.1)
#container for p values
perm_p_values <- matrix(NA, nrow = n_sims, ncol = length(shift))
normal_p_values <- matrix(NA, nrow = n_sims, ncol = length(shift))
for(i in 1:n_sims){
  for(j in 1:length(shift)){
    sample_1 <- sample(topsoil_TOC_paicines, size = sample_size, replace = TRUE)
    sample_2 <- sample(topsoil_TOC_paicines + shift[j], size = sample_size, replace = TRUE)
    diff_mean <- mean(sample_1) - mean(sample_2)
    normal_p_values[i,j] <- t.test(x = sample_1, y = sample_2, alternative = "two.sided")$p.value
    perm_p_values[i,j] <- t2p(diff_mean, two_sample(x = sample_1, y = sample_2, reps = 200), alternative = "two-sided")
  }
}

normal_power_shift <- colMeans(normal_p_values < .05)
perm_power_shift <- colMeans(perm_p_values < .05)

#the normal theory and permutation power functions are almost identical
plot(x = shift, y = perm_power_shift, type = 'l', col = 'blue', lwd = 3, ylim = c(0,1))
points(x = shift, y = normal_power_shift, type = 'l', col = 'red', lwd = 3)


############## analyze advantages of stratified versus uniform random sampling ###########
#we will empirically investigate the advantages of stratified sampling in these settings
#the idea is similar to a bootstrap. 
#take the samples as a population and then simulate taking subsamples

#first we analyze stratified sampling at a small scale by using the Paicines samples
#the strata are defined by transects, which roughly correspond to different locations on the ranch
paicines_data_topsoil <- paicines_master %>% 
  filter(depth == "a") %>%
  filter(!is.na(TOC))

#the sample size for simulations
n <- 100
#function to return the estimated mean and standard error given a population and sample index
get_mean_se <- function(population, sample_index){
  c(mean(population[sample_index]), sd(population[sample_index])/sqrt(length(sample_index)))
}
#function to return estimated mean and standard error given a population, sampling index, and stratification information
get_mean_se_stratified <- function(population, strata_sizes, sample_index, strata){
  strata_sample_sizes <- as.numeric(table(strata))
  strata_means <- tapply(population[sample_index], strata, mean)
  strata_vars <- tapply(population[sample_index], strata, var)
  var_estimate <- sum((strata_sizes / length(population))^2 * strata_vars / strata_sample_sizes)
  c(mean(strata_means), sqrt(var_estimate))
}

#function to run a single simulation on the Paicines data
run_paicines_sim <- function(data_frame){
  N_strata <- as.numeric(table(data_frame$transect))
  #proportional allocation to strata
  n_strata <- round(n * N_strata / sum(N_strata))
  if(sum(n_strata) < n){
    n_strata[length(n_strata)] <- n_strata[length(n_strata)] + 1
  } 
  if(sum(n_strata) > n){
    n_strata[length(n_strata)] <- n_strata[length(n_strata)] - 1
  }
  proportional_stratified_sample <- strata(data = data_frame, stratanames = "transect", size = n_strata, method = "srswr")
  
  random_sample <- sample(1:nrow(data_frame), size = n, replace = TRUE)
  
  
  stratified_estimates <- get_mean_se_stratified(data_frame$TOC, N_strata, proportional_stratified_sample$ID_unit, proportional_stratified_sample$Stratum)
  random_estimates <- get_mean_se(data_frame$TOC, random_sample)
  #local_pivotal_estimates <- get_mean_se(population, local_pivotal_sample)
  
  cbind(stratified_estimates, random_estimates)
}
#compute the estimand, i.e. the true mean in the Paicines topsoil data
true_paicines_mean <- mean(paicines_data_topsoil$TOC)

#run simulations, replicated 1000 times
paicines_sims <- replicate(n = 1000, run_paicines_sim(data_frame = paicines_data_topsoil))

#compute properties of the samples
#the empirical coverage of 95% normal theory confidence intervals (should be 95%)
coverage_paicines <- apply(paicines_sims[1,,] - qnorm(p = .975) * paicines_sims[2,,] <= true_paicines_mean & true_paicines_mean <= paicines_sims[1,,] + qnorm(p = .975) * paicines_sims[2,,], 1, mean)
#the width of a confidence interval, basically 4 times the SE
ci_width_paicines <- apply(2 * qnorm(p = .975) * paicines_sims[2,,], 1, mean)
#the actual standard error over simulations
se_paicines <- apply(paicines_sims, c(1,2), sd)[1,]
#the average estimated standard error
se_hat_paicines <- apply(paicines_sims, c(1,2), mean)[2,]
#the mean squared error
rmse_paicines <- sqrt(apply((paicines_sims - true_paicines_mean)^2, c(1,2), mean)[1,])
mad_paicines <- apply(abs(paicines_sims - true_paicines_mean), c(1,2), median)[1,]
#ratio of mse compared to uniform independent random sampling
#rmse_paicines <- rmse_paicines / rmse_paicines[2]
#compile results
paicines_results_frame <- rbind(coverage_paicines, ci_width_paicines, rmse_paicines, mad_paicines) %>%
  as_tibble() %>%
  mutate(property = c("coverage", "ci width", "RMSE", "MAD")) %>%
  select(property, random_estimates, stratified_estimates)


#run simulations on the cropland data, considering stratification by site compared to randomly resampling from all points.
crop_data_topsoil <- cropland_master %>%
  filter(depth == "a") %>%
  filter(!is.na(TC))

run_cropland_sim <- function(data_frame){
  N_strata <- as.numeric(table(data_frame$site))
  #sample size in each strata should be about proportional to the size of each strata
  n_strata <- c(5, 5, 9, 4, 4, 3)
  proportional_stratified_sample <- strata(data = data_frame, stratanames = "site", size = n_strata, method = "srswr")
  
  random_sample <- sample(1:nrow(data_frame), size = n, replace = FALSE)
  
  stratified_estimates <- get_mean_se_stratified(data_frame$TC, N_strata, proportional_stratified_sample$ID_unit, proportional_stratified_sample$Stratum)
  random_estimates <- get_mean_se(data_frame$TC, random_sample)
  #local_pivotal_estimates <- get_mean_se(population, local_pivotal_sample)
  
  cbind(stratified_estimates, random_estimates)
}


true_cropland_mean <- mean(crop_data_topsoil$TC)

crop_sims <- replicate(n = 1000, run_cropland_sim(data_frame = crop_data_topsoil))
coverage_crop <- apply(crop_sims[1,,] - qnorm(p = .975) * crop_sims[2,,] <= true_cropland_mean & true_cropland_mean <= crop_sims[1,,] + qnorm(p = .975) * crop_sims[2,,], 1, mean)
ci_width_crop <- apply(2 * qnorm(p = .975) * crop_sims[2,,], 1, mean)
se_crop <- apply(crop_sims, c(1,2), sd)[1,]
se_hat_crop <- apply(crop_sims, c(1,2), mean)[2,]
mse_crop <- apply((crop_sims - true_cropland_mean)^2, c(1,2), mean)[1,]
#ratio of mse compared to simple random sampling
mse_crop <- mse_crop / mse_crop[2]

crop_results_frame <- rbind(coverage_crop, ci_width_crop, mse_crop) %>%
  as_tibble() %>%
  mutate(property = c("coverage", "ci width", "relative MSE")) %>%
  select(property, random_estimates, stratified_estimates)




