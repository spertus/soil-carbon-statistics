library(tidyverse)
library(readxl)
library(permuter)
library(sampling)
#permuter is used for permutation tests. It can be installed from Github: https://github.com/statlab/permuter
source("functions.R")

########### read in data from excel spreadsheet(s) ##############
#data from Paicines Ranch
#from the original data sent to me by Paige Stanley, I deleted a few free-floating cells on the Rangeland_BD sheet. I also changed 'dulk_density_g/cm3' to 'bulk_density_g/cm3' on Rangeland_BD and Cropland_BD.
paicines_master <- read_excel("../Data/Heterogeneity_Master_PS_02032021.xlsx", sheet = "Rangeland_All_soliTOC") %>%
  rename(SOC = TOC) %>%
  mutate(SOC = ifelse(SOC == "NA", NA, SOC)) %>%
  mutate(TIC = ifelse(TIC == "NA", NA, TIC)) %>%
  mutate(TC = ifelse(TC == "NA", NA, TC)) %>%
  mutate(SOC = as.numeric(SOC), TIC = as.numeric(TIC), TC = as.numeric(TC), sample_number = as.numeric(sample_number))

paicines_solitoc_reps <- read_excel("../Data/Heterogeneity_Master_PS_02032021.xlsx", sheet = "Rangeland_Reps_soliTOC") %>%
  rename(SOC = TOC) %>%
  mutate(SOC = ifelse(SOC == "NA", NA, SOC)) %>%
  mutate(SOC = as.numeric(SOC), sample_number = as.numeric(sample_number)) %>%
  select(-TOC_replicate_avg, -TOC_replicate_stdev) %>%
  mutate(machine = "solitoc")

paicines_costech_reps <- read_excel("../Data/Heterogeneity_Master_PS_02032021.xlsx", sheet = "Rangeland_Reps_Costech") %>%
  rename(TC = TOC) %>%
  mutate(TC = ifelse(TC == "NA", NA, TC)) %>%
  mutate(TC = as.numeric(TC), sample_number = as.numeric(sample_number)) %>%
  select(-replicate_avg, -replicate_stdev) %>%
  mutate(machine = "costech")

paicines_BD <- read_excel("../Data/Heterogeneity_Master_PS_02032021.xlsx", sheet = "Rangeland_BD") %>%
  rename('bd' = 'bulk_density_g/cm3') %>%
  mutate(bd = as.numeric(bd))

#data from various croplands around California
cropland_master <- read_excel("../Data/Heterogeneity_Master_PS_02032021.xlsx", sheet = "Cropland_All_Costech") %>%
  rename(TC = TOC) %>%
  mutate(TC = ifelse(TC == "NA", NA, TC)) %>%
  mutate(TC = as.numeric(TC), sample_number = as.numeric(sample_number))

cropland_solitoc_reps <- read_excel("../Data/Heterogeneity_Master_PS_02032021.xlsx", sheet = "Cropland_Reps_soliTOC") %>%
  rename(SOC = TOC) %>%
  mutate(SOC = ifelse(SOC == "NA", NA, SOC)) %>%
  mutate(SOC = as.numeric(SOC), sample_number = as.numeric(sample_number)) %>%
  select(-TOC_replicate_avg, -TOC_replicate_stdev) %>% 
  mutate(machine = "solitoc")

cropland_costech_reps <- read_excel("../Data/Heterogeneity_Master_PS_02032021.xlsx", sheet = "Cropland_Reps_Costech") %>%
  rename(TC = TOC) %>%
  mutate(TC = ifelse(TC == "NA", NA, TC)) %>%
  mutate(TC = as.numeric(TC), sample_number = as.numeric(sample_number)) %>%
  select(-replicate_avg, -replicate_stdev) %>% 
  mutate(machine = "costech")

cropland_BD <- read_excel("../Data/Heterogeneity_Master_PS_02032021.xlsx", sheet = "Cropland_BD")  %>%
  rename('bd' = 'bulk_density_g/cm3') %>%
  mutate(bd = as.numeric(bd))



###################### SOC concentration and BD in space ######################

#stacked histograms for Paicines data
ggplot(paicines_master, aes(TC)) +
  geom_histogram() +
  facet_grid(depth ~ transect)


#stacked histograms for cropland data
ggplot(cropland_master, aes(TC)) +
  geom_histogram() +
  facet_grid(depth ~ site)


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
  facet_grid(depth ~ .)

ggplot(cropland_BD, aes(bd)) +
  geom_histogram() +
  facet_grid(depth ~ .)

paicines_summary_bd <- paicines_BD %>%
  group_by(depth, transect) %>%
  summarize(mean_bd = mean(bd, na.rm = TRUE), sd_bd = sd(bd, na.rm = TRUE))

cropland_summary_bd <- cropland_BD %>%
  group_by(depth, site) %>%
  summarize(mean_bd = mean(bd, na.rm = TRUE), sd_bd = sd(bd, na.rm = TRUE))

############## replicates and assay error ##############
assay_error_solitoc <- paicines_solitoc_reps %>% 
  bind_rows(cropland_solitoc_reps) %>%
  group_by(site, sample_number) %>%
  summarize(
    machine = first(machine),
    depth = first(depth), 
    sigma_delta = sqrt(var(TC) / (mean(TC)^2 - var(TC)/n())),
    TC = mean(TC)
  ) %>%
  ungroup() %>%
  na.omit()

assay_error_costech <- paicines_costech_reps %>% 
  bind_rows(cropland_costech_reps) %>%
  group_by(site, sample_number) %>%
  summarize(
    machine = first(machine),
    depth = first(depth), 
    sigma_delta = sqrt(var(TC) / (mean(TC)^2 - var(TC)/n())),
    TC = mean(TC)
  ) %>%
  na.omit()

assay_error_long <- bind_rows(assay_error_solitoc, assay_error_costech)

#compare average errors on paicines samples
ggplot(assay_error_long, aes(sigma_delta*100, fill = machine)) +
  geom_density(alpha = .5) +
  xlim(0,10) +
  xlab("Percent Error")

sigma_delta_solitoc <- median(assay_error_solitoc$sigma_delta)
sigma_delta_costech <- median(assay_error_costech$sigma_delta)

#compare measurements
assay_error_wide <- assay_error_long %>%
  pivot_wider(names_from = "machine", values_from = c("sigma_delta", "TC"))

ggplot(assay_error_wide, aes(x = TC_solitoc, y = TC_costech, color = site)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  xlab("%TC SoliTOC") +
  ylab("%TC Costech")
cor(assay_error_wide$TC_solitoc, assay_error_wide$TC_costech, method = "spearman")


ggplot(assay_error_wide, aes(x = sigma_delta_solitoc, y = sigma_delta_costech, color = site)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  xlab("% Error SoliTOC") +
  ylab("% Error Costech") +
  ylim(0,.25) +
  xlim(0,.25)

median(assay_error_solitoc$TC)
median(assay_error_costech$TC)


#non-parametric permutation test for differences in measurement
B <- 10000

#nonparametric analysis of no difference in labs/machines:
reps_long <- bind_rows(paicines_solitoc_reps, cropland_solitoc_reps) %>%
  select(-SOC, -TIC) %>%
  bind_rows(paicines_costech_reps %>% rename(TC = TC), cropland_costech_reps %>% rename(TC = TC)) %>%
  na.omit() %>%
  mutate(sample_ID = paste(site, sample_number, sep = "_"))


strata <- as_factor(reps_long$sample_ID)
num_strata <- length(unique(strata))
#test statistic is difference in means between labs
diff_means <- rep(0, num_strata)
pooled_ses <- rep(0, num_strata)
null_distributions <- matrix(0, ncol = num_strata, nrow = B)
p_values <- rep(0, num_strata)
#iterate across strata (samples)
for(k in 1:num_strata){
  solitoc_TC <- reps_long$TC[reps_long$machine == "solitoc" & strata == levels(strata)[k]]
  costech_TC <- reps_long$TC[reps_long$machine == "costech" & strata == levels(strata)[k]]
  diff_means[k] <- mean(solitoc_TC) - mean(costech_TC)
  pooled_ses[k] <- sqrt(var(solitoc_TC) / length(solitoc_TC) + var(costech_TC) / length(solitoc_TC))
  null_distributions[,k] <- two_sample(x = solitoc_TC, y = costech_TC, reps = B)
  p_values[k] <- mean(abs(diff_means[k]) <= abs(null_distributions[,k])) 
}
fisher_p_value <- pchisq(-2 * sum(log(p_values)), df = 2 * num_strata, lower.tail = F)
npc_p_value <- npc(diff_means, distr = null_distributions, combine=  "fisher", alternatives = "two-sided")


#average inorganic carbon is interesting as a predictor of difference between costech and solitoc
avg_TIC <- bind_rows(paicines_solitoc_reps, cropland_solitoc_reps) %>%
  mutate(sample_ID = paste(site, sample_number, sep = "_")) %>%
  group_by(sample_ID) %>%
  summarize(avg_TIC = mean(TIC), avg_TC = mean(TC), site = first(site))

p_value_frame <- data.frame(sample_ID = levels(strata), diff_mean = diff_means, t_stat = diff_means/pooled_ses, permutation_p_value = p_values) %>%
  left_join(avg_TIC, by = "sample_ID") %>%
  as_tibble()

ggplot(p_value_frame, aes(x = log10(avg_TC), y = t_stat, shape = permutation_p_value < .05, color = site)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = qt(p = .975, df = 9), linetype = 'dashed') + 
  geom_hline(yintercept = qt(p = .025, df = 9), linetype = 'dashed') +
  geom_point()

############### costs and optimal compositing #############
#Costs based on 'Cost of Carbon Analysis.xlsx' assembled by Jessica Chiartas
#in sampling, about 20-24 samples (4 transects * 5-6 composites) costs 400$ at 25$/hr and 16hrs, so we'll take 20 USD per sample.
#the commercial labs had costs of 19.5, 7.75, 22, 10, and 8.75 USD per sample to assess SOC. So we'll take the cost of measurement with dry combustion as 13.60, including sample prep.
#the in-house cost of measurement was 2.78


#sample size
sample_size <- 60
#max budget, i.e. the budget to sample and measure sample_size samples, with no compositing
max_budget <- 20 * sample_size + 13.6 * sample_size 
#sample size if we measure once and devote the rest of the budget to sampling
max_sample_size <- floor((max_budget - 13.6) / 20)

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

#power to detect a half-percentage point change in SOC
power_change_paicines <- paicines_master %>%
  group_by(depth) %>%
  summarize(mu = mean(TC, na.rm = T), sigma_p = sd(TC, na.rm = T)) %>%
  mutate(power_no_compositing = get_power_two_sample(n_1 = sample_size, k_1 = sample_size, n_2 = sample_size, k_2 = sample_size, mu_1 = mu, sigma_p_1 = sigma_p, mu_2 = mu + 0.5, sigma_p_2 = sigma_p, sigma_delta = sigma_delta_costech)) %>%
  mutate(power_full_compositing = get_power_two_sample(n_1 = max_sample_size, k_1 = 1, n_2 = max_sample_size, k_2 = 1, mu_1 = mu, sigma_p_1 = sigma_p, mu_2 = mu + 0.5, sigma_p_2 = sigma_p, sigma_delta = sigma_delta_costech))
  


############## analyze advantages of stratified versus uniform random sampling ###########
#we will empirically investigate the advantages of stratified sampling in these settings
#the idea is similar to a bootstrap. 
#take the samples as a population and then simulate taking subsamples

#first we analyze stratified sampling at a small scale by using the Paicines samples
#the strata are defined by transects, which roughly correspond to different locations on the ranch
paicines_data_topsoil <- paicines_master %>% 
  filter(depth == "a") %>%
  filter(!is.na(SOC))

#the sample size for simulations
n <- 30
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
  #technically there are slightly fewer samples in transect Bx
  n_strata <- c(6, 6, 6, 6, 6)
  proportional_stratified_sample <- strata(data = data_frame, stratanames = "transect", size = n_strata, method = "srswr")
  
  random_sample <- sample(1:nrow(data_frame), size = n, replace = TRUE)
  
  
  stratified_estimates <- get_mean_se_stratified(data_frame$SOC, N_strata, proportional_stratified_sample$ID_unit, proportional_stratified_sample$Stratum)
  random_estimates <- get_mean_se(data_frame$SOC, random_sample)
  #local_pivotal_estimates <- get_mean_se(population, local_pivotal_sample)
  
  cbind(stratified_estimates, random_estimates)
}
#compute the estimand, i.e. the true mean in the Paicines topsoil data
true_paicines_mean <- mean(paicines_data_topsoil$SOC)

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
mse_paicines <- apply((paicines_sims - true_paicines_mean)^2, c(1,2), mean)[1,]
#ratio of mse compared to uniform independent random sampling
mse_paicines <- mse_paicines / mse_paicines[2]
paicines_results_frame <- rbind(coverage_paicines, ci_width_paicines, mse_paicines) %>%
  as_tibble() %>%
  mutate(property = c("coverage", "ci width", "relative MSE")) %>%
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




