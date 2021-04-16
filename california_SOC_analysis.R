library(tidyverse)
library(readxl)
library(permuter)
library(sampling)
library(geoR)
#permuter is used for permutation tests. It can be installed from Github: https://github.com/statlab/permuter
source("functions.R")

########### read in data from excel spreadsheet(s) ##############
#note that some of the conversions to numeric cause NA warnings, I have spot checked these against the excel files, there does not appear to be a real issue.
rangeland_master <- read_excel("R_Heterogeneity_Master_PS_04132021.xlsx", sheet = "Rangeland_All_soliTOC") %>%
  mutate(TOC = ifelse(TOC == "NA", NA, TOC)) %>%
  mutate(TIC = ifelse(TIC == "NA", NA, TIC)) %>%
  mutate(TC = ifelse(TC == "NA", NA, TC)) %>%
  mutate(TOC = as.numeric(TOC), TIC = as.numeric(TIC), TC = as.numeric(TC), sample_number = as.numeric(sample_number)) %>%
  mutate(depth_long = dplyr::recode(depth, a = "0-10 cm", b = "10-30 cm", c = "30-50 cm", d = "50-75 cm", e = "75-100 cm"))

rangeland_solitoc_reps <- read_excel("R_Heterogeneity_Master_PS_04132021.xlsx", sheet = "Rangeland_Reps_soliTOC") %>%
  mutate(TOC = ifelse(TOC == "NA", NA, TOC)) %>%
  mutate(TOC = as.numeric(TOC), sample_number = as.numeric(sample_number)) %>%
  mutate(machine = "solitoc")

rangeland_costech_reps <- read_excel("R_Heterogeneity_Master_PS_04132021.xlsx", sheet = "Rangeland_Reps_Costech") %>%
  mutate(TC = ifelse(TC == "NA", NA, TC), TOC = ifelse(TOC == "NA", NA, TOC)) %>%
  mutate(TC = as.numeric(TC), TOC = as.numeric(TOC), sample_number = as.numeric(sample_number)) %>%
  mutate(machine = "costech")

rangeland_BD <- read_excel("R_Heterogeneity_Master_PS_04132021.xlsx", sheet = "Rangeland_BD") %>%
  rename('bd' = 'bulk_density_g/cm3') %>%
  mutate(bd = as.numeric(bd)) %>%
  mutate(depth_long = dplyr::recode(depth, a = "0-10 cm", b = "10-30 cm", c = "30-50 cm", d = "50-75 cm", e = "75-100 cm"))

#data from various croplands around California
cropland_master <- read_excel("R_Heterogeneity_Master_PS_04132021.xlsx", sheet = "Cropland_All_Costech") 

cropland_solitoc_reps <- read_excel("R_Heterogeneity_Master_PS_04132021.xlsx", sheet = "Cropland_Reps_soliTOC") %>%
  mutate(sample_number = as.numeric(sample_number)) %>%
  mutate(machine = "solitoc")

cropland_costech_reps <- read_excel("R_Heterogeneity_Master_PS_04132021.xlsx", sheet = "Cropland_Reps_Costech") %>%
  mutate(TC = ifelse(TC == "NA", NA, TC), TOC = ifelse(TOC == "NA", NA, TOC)) %>%
  mutate(TC = as.numeric(TC), TOC = as.numeric(TOC), sample_number = as.numeric(sample_number)) %>%
  mutate(machine = "costech")

cropland_BD <- read_excel("R_Heterogeneity_Master_PS_04132021.xlsx", sheet = "Cropland_BD")  %>%
  rename('bd' = 'bulk_density_g/cm3') %>%
  mutate(bd = as.numeric(bd))

replicates_comparison <- read_excel("R_Heterogeneity_Master_PS_04132021.xlsx", sheet = "Replicate_Comparison") %>% 
  rename(TC_solitoc = TC_soliTOC, TOC_solitoc = TOC_soliTOC) %>%
  mutate_at(vars(starts_with(c("TC_", "TOC_"))), as.numeric)

standards_comparison <- read_excel("R_Heterogeneity_Master_PS_04132021.xlsx", sheet = "Standards_Comparison") %>%
  rename(TC = `TC%`, N = `N%`)


combined_master <- rangeland_master %>% 
  select(site, soil_type, depth, TC) %>%
  bind_rows(cropland_master %>% select(site, soil_type, depth, TC)) %>%
  mutate(land_use = ifelse(site == "PAIC", "Rangeland", "Cropland"))

###################### TOC concentration and BD in space ######################

#stacked histograms for rangeland data
ggplot(rangeland_master, aes(TC)) +
  geom_histogram(binwidth = 0.25) +
  facet_grid(depth ~ transect) +
  xlab("% Total Carbon (TC)") +
  ylab("Number of Samples") +
  ylim(0,23) +
  theme_bw() +
  theme(text = element_text(size = 16))


#stacked histograms for cropland data
ggplot(cropland_master, aes(TC)) +
  geom_histogram(binwidth = 0.25) +
  facet_grid(depth ~ site) +
  xlab("% Total Carbon (TC)") +
  ylab("Number of Samples") +
  theme_bw() +
  xlim(0,8) +
  ylim(0,23) +
  theme(text = element_text(size = 16))


#summary tables
rangeland_summary <- rangeland_master %>% 
  group_by(depth, transect) %>%
  summarize(mean_TOC = mean(TOC, na.rm = T), sd_TOC = sd(TOC, na.rm = T)) %>%
  mutate(sd_TOC = paste("(", round(sd_TOC,2), ")", sep = "")) %>%
  mutate(mean_TOC = round(mean_TOC, 2)) %>%
  unite(col = "mean (SD)", mean_TOC, sd_TOC, sep = " ") %>%
  pivot_wider(names_from = transect, values_from = "mean (SD)")

#rangeland data combined across transects
rangeland_summary_overall <- rangeland_master %>% 
  group_by(depth) %>%
  summarize(mean_TOC = mean(TOC, na.rm = T), sd_TOC = sd(TOC, na.rm = T)) %>%
  mutate(sd_TOC = paste("(", round(sd_TOC,2), ")", sep = "")) %>%
  mutate(mean_TOC = round(mean_TOC, 2)) %>%
  unite(col = "mean (SD)", mean_TOC, sd_TOC, sep = " ")

#note cropland is total carbon not total organic carbon
cropland_summary <- cropland_master %>% 
  group_by(depth, site) %>%
  summarize(mean_TC = mean(TC, na.rm = T), sd_TC = sd(TC, na.rm = T))


#linear models
rangeland_model <- lm(TC ~ transect*depth, data = rangeland_master)
cropland_model <- lm(TC ~ site*depth, data = cropland_master)

anova(rangeland_model)
anova(cropland_model)


#bulk densities
#NOTE: histograms are not disaggregated by location, which drives variation
ggplot(rangeland_BD, aes(bd)) +
  geom_histogram(binwidth = 0.05) +
  facet_grid(depth ~ .) +
  xlab("Bulk Density") +
  ylab("Number of samples") +
  scale_x_continuous(breaks = c(0.8, 1.2, 1.6, 2.0, 2.4), limits = c(0.8,2.4)) +
  scale_y_continuous(breaks = c(0,2,4,6,8,10), limits = c(0,10)) +
  theme_bw() +
  theme(text = element_text(size = 16), panel.grid.minor = element_blank())


ggplot(cropland_BD, aes(bd)) +
  geom_histogram(binwidth = 0.05) +
  facet_grid(depth ~ .) +
  xlab("Bulk Density") +
  ylab("Number of Samples") +
  scale_x_continuous(breaks = c(0.8, 1.2, 1.6, 2.0, 2.4), limits = c(0.8,2.4)) +
  scale_y_continuous(breaks = c(0,2,4,6,8,10), limits = c(0,10)) +
  theme_bw() +
  theme(text = element_text(size = 16))

rangeland_summary_bd <- rangeland_BD %>%
  group_by(depth) %>%
  summarize(mean_bd = mean(bd, na.rm = TRUE), sd_bd = sd(bd, na.rm = TRUE)) %>%
  mutate(cv_bd = sd_bd / mean_bd)

#averages within sites
cropland_summary_bd <- cropland_BD %>%
  group_by(depth, site) %>%
  summarize(mean_bd = mean(bd, na.rm = TRUE), sd_bd = sd(bd, na.rm = TRUE)) %>%
  mutate(cv_bd = sd_bd / mean_bd) %>% 
  group_by(depth) %>% 
  summarize(mean_bd = mean(mean_bd, na.rm = TRUE), sd_bd = sd(sd_bd, na.rm = TRUE), cv_bd = mean(cv_bd, na.rm = TRUE))



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
  

#compare average errors on rangeland samples
ggplot(assay_error_long, aes(sigma_delta_TC*100, fill = machine)) +
  geom_density(alpha = .5) +
  xlim(0,10) +
  xlab("Percent Error") +
  ylab("Density") +
  theme_bw() +
  scale_fill_discrete(name = "Machine") +
  theme(text = element_text(size = 16))

median_sigma_delta <- assay_error_long %>%
  group_by(machine) %>%
  summarize(sigma_delta = median(sigma_delta_TC))



#non-parametric permutation test for differences in measurement
#nonparametric analysis of no difference in labs/machines:
reps_long_TC <- replicates_comparison %>% 
  select(site, sample_number, TC_solitoc, TC_costech) %>%
  pivot_longer(cols = c("TC_solitoc", "TC_costech"), names_prefix = "TC_", names_to = "machine", values_to = "carbon") %>%
  na.omit() %>% 
  mutate(carbon_type = "TC")

#there are only a few measurements for TOC on the costech
reps_long_TOC <- replicates_comparison %>% 
  select(site, sample_number, TOC_solitoc, TOC_costech) %>%
  group_by(sample_number) %>%
  filter(sum(!is.na(TOC_costech)) >= 2) %>%
  ungroup() %>%
  pivot_longer(cols = c("TOC_solitoc", "TOC_costech"), names_prefix = "TOC_", names_to = "machine", values_to = "carbon") %>%
  na.omit() %>% 
  mutate(carbon_type = "TOC")

#combined TC and TOC reps
reps_long_both <- reps_long_TOC %>%
  bind_rows(reps_long_TC) %>% 
  mutate(identifier = paste(site, sample_number, carbon_type,  sep = "_"))



B <- 10000
strata <- as_factor(reps_long_both$identifier)
num_strata <- length(unique(strata))
#test statistic is difference in means between labs
diff_means <- rep(0, num_strata)
pooled_ses <- rep(0, num_strata)
null_distributions <- matrix(0, ncol = num_strata, nrow = B)
p_values <- rep(0, num_strata)
#iterate across strata (samples)
for(k in 1:num_strata){
  solitoc_carbon <- reps_long_both$carbon[reps_long_both$machine == "solitoc" & strata == levels(strata)[k]]
  costech_carbon <- reps_long_both$carbon[reps_long_both$machine == "costech" & strata == levels(strata)[k]]
  diff_means[k] <- mean(solitoc_carbon) - mean(costech_carbon)
  pooled_ses[k] <- sqrt(var(solitoc_carbon) / length(solitoc_carbon) + var(costech_carbon) / length(solitoc_carbon))
  null_distributions[,k] <- two_sample(x = solitoc_carbon, y = costech_carbon, reps = B)
  p_values[k] <- mean(abs(diff_means[k]) <= abs(null_distributions[,k])) 
}

#this is a non-parametric test of whether *any* of the measurements are significantly different, using nonparametric combination of tests
npc_p_value <- npc(diff_means, distr = null_distributions, combine = "fisher", alternatives = "two-sided")

reps_long_both <- reps_long_both %>%
  mutate(rejected = ifelse(identifier %in% levels(strata)[p_values < .05], 1, 0))

#plot assay densities stratified by sample number
assay_density_plot_accepted <- ggplot(reps_long_both %>% filter(rejected == 0), aes(carbon, fill = machine)) +
  geom_density(alpha = .5) +
  facet_wrap(~ identifier, scales = "free_y") +
  xlim(0,6) +
  xlab("Percent Organic Carbon") +
  ggtitle("Similar Assays") +
  theme_bw() +
  scale_fill_discrete(name = "Machine") +
  theme(text = element_text(size = 16), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())
assay_density_plot_rejected <- ggplot(reps_long_both %>% filter(rejected == 1, carbon_type == "TOC"), aes(carbon, fill = machine)) +
  geom_density(alpha = .5) +
  facet_wrap(~ identifier, scales = "free_y") +
  xlab("Percent Organic Carbon") +
  ggtitle("Significantly Different Assays") +
  theme_bw() +
  scale_fill_discrete(name = "Machine") +
  theme(text = element_text(size = 16), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())


#standards comparison
known_standards <- standards_comparison %>%
  group_by(sample_ID) %>%
  summarize(known_TC = first(known_TC))
standards_comparison <- standards_comparison %>%
  mutate(known_TC = ifelse(sample_ID == "EML", 1.86, 0.926))

standards_density_plot <- ggplot(standards_comparison, aes(TC, fill = machine)) +
  geom_density(alpha = .5) +
  geom_vline(data = standards_comparison %>% filter(sample_ID == "EML"), aes(xintercept = 1.86)) +
  geom_vline(data = standards_comparison %>% filter(sample_ID == "LECO"), aes(xintercept = 0.926)) +
  facet_grid(~ sample_ID, scales = "free") +
  xlim(0.8, 2.2) +
  xlab("% Total Carbon") +
  ylab("") +
  theme_bw()


############# spatial correlation of TOC concentrations at rangeland ##########
#first just do for a single transect
T_topsoil_TOC <- rangeland_master %>% 
  filter(transect == "T", depth == "a") %>%
  arrange(sample_number) %>%
  pull(TOC)

#dataframe for use with geodata
T_topsoil_geodata <- as.geodata(data.frame(x = 1:length(T_topsoil_TOC), y = 1, TOC = T_topsoil_TOC), coords.col = 1:2, data.col = 3)
T_topsoil_variogram <- variog(T_topsoil_geodata, option = "bin")
plot(T_topsoil_variogram)

#function to plot a variogram for any depth and transect
#inputs:
  #rangeland_depth: a string denoting the depth, either "a", "b", "c", "d", or "e"
  #rangeland_transect: a string denoting the transect, eitehr "T", "Mx", "My", "Bx", or "By"
#output:
  #a plot of the empirical variogram generated by binning variances between points at a range of distances
plot_variogram <- function(rangeland_depth = "a", rangeland_transect = "T", plot = TRUE){
  if(!(rangeland_depth %in% c("a","b","c","d","e")) | !all(rangeland_transect %in% c("T","Mx","My","Bx","By"))){
    stop("Depth or transect is invalid! See comments.")
  }
  TOC <- rangeland_master %>% 
    filter(transect == rangeland_transect, depth == rangeland_depth) %>%
    arrange(sample_number) %>%
    pull(TOC)
  
  #dataframe for use with geodata
  geodataframe <- as.geodata(data.frame(x = 3 * 1:length(TOC), y = 1, TOC = TOC), coords.col = 1:2, data.col = 3)
  variogram <- variog(geodataframe, option = "bin")
  if(plot){
    plot(variogram, xlab = "Distance (meters)", ylab = "Semivariance", main = paste("Variogram for transect", rangeland_transect, "and depth", rangeland_depth))
  } else{
    variogram
  }
}


variog_T <- plot_variogram(rangeland_depth = "a", rangeland_transect = "T", plot = FALSE)
variog_Mx <- plot_variogram(rangeland_depth = "a", rangeland_transect = "Mx", plot = FALSE)
variog_My <- plot_variogram(rangeland_depth = "a", rangeland_transect = "My", plot = FALSE)
variog_Bx <- plot_variogram(rangeland_depth = "a", rangeland_transect = "Bx", plot = FALSE)
variog_By <- plot_variogram(rangeland_depth = "a", rangeland_transect = "By", plot = FALSE)

distance <- variog_T$u
avg_variogram <- colMeans(rbind(variog_T$v, variog_Mx$v, variog_My$v, variog_Bx$v, variog_By$v))


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
#measurement errors of costech and solitoc
sigma_delta_costech <- median_sigma_delta %>% 
  filter(machine == "costech") %>%
  pull(sigma_delta)
sigma_delta_solitoc <- median_sigma_delta %>% 
  filter(machine == "solitoc") %>%
  pull(sigma_delta)

#precision of the sample mean 
precision_combined <- combined_master %>%
  group_by(depth, site, land_use) %>%
  summarize(mu_site = mean(TC, na.rm = T), sigma_p_site = sd(TC, na.rm = T)) %>%
  group_by(depth, land_use) %>%
  summarize(mu = mean(mu_site, na.rm = T), sigma_p = mean(sigma_p_site, na.rm = T)) %>%
  mutate(se_no_compositing = sqrt(get_variance(n = sample_size, k = sample_size, sigma_p = sigma_p, mu = mu, sigma_delta = sigma_delta_costech))) %>%
  mutate(se_full_compositing = sqrt(get_variance(n = max_sample_size, k = 1, sigma_p = sigma_p, mu = mu, sigma_delta = sigma_delta_costech))) %>%
  mutate(se_optimal_compositing = sqrt(get_minimum_error(sigma_p = sigma_p, mu = mu, sigma_delta = sigma_delta_costech, C_0 = 0, cost_c = 20, cost_P = 0, cost_A = 13.60, B = max_budget)$optimum_variance)) %>%
  mutate(composite_efficiency_ratio = se_optimal_compositing / se_no_compositing) %>%
  mutate(optimal_composite_size_commercial = get_optimal_composite_size(sigma_p = sigma_p, mu = mu, sigma_delta = sigma_delta_costech, cost_c = 20, cost_P = 0, cost_A = 13.60)) %>%
  mutate(optimal_composite_size_inhouse = get_optimal_composite_size(sigma_p = sigma_p, mu = mu, sigma_delta = sigma_delta_costech, cost_c = 20, cost_P = 0, cost_A = 2.78)) %>%
  arrange(land_use, depth)


#power of two-sample t-test to detect a range of changes in rangeland topsoil

mu_0_rangeland <- precision_combined %>% filter(depth == "a", land_use == "Rangeland") %>% pull(mu)
sigma_p_rangeland <- precision_combined %>% filter(depth == "a", land_use == "Rangeland") %>% pull(sigma_p)
mu_0_cropland <- precision_combined %>% filter(depth == "a", land_use == "Cropland") %>% pull(mu)
sigma_p_cropland <- precision_combined %>% filter(depth == "a", land_use == "Cropland") %>% pull(sigma_p)

rangeland_grid <- expand.grid(
  land_use = "Rangeland",
  sample_size = c(5,10,30,90),
  mu_0 = mu_0_rangeland, 
  sigma_p = sigma_p_rangeland, 
  sigma_delta = c(sigma_delta_costech, sigma_delta_solitoc),
  delta = seq(0,1.5,by=.01)
)
cropland_grid <- expand.grid(
  land_use = "Cropland",
  sample_size = c(5,10,30,90),
  mu_0 = mu_0_cropland, 
  sigma_p = sigma_p_cropland, 
  sigma_delta = c(sigma_delta_costech, sigma_delta_solitoc),
  delta = seq(0,1.5,by=.01)
)


power_change_topsoil <- rangeland_grid %>%
  bind_rows(cropland_grid) %>%
  mutate(Machine = ifelse(sigma_delta == sigma_delta_costech, "Costech", "SoliTOC")) %>%
  mutate(budget = 20 * sample_size + 13.6 * sample_size) %>%
  mutate(max_sample_size = floor((budget - 13.6) / 20)) %>%
  mutate(opt_n = get_minimum_error(sigma_p = sigma_p, mu = mu_0, sigma_delta = sigma_delta, C_0 = 0, cost_c = 20, cost_P = 0, cost_A = 13.60, B = budget)$n) %>%
  mutate(opt_k = get_minimum_error(sigma_p = sigma_p, mu = mu_0, sigma_delta = sigma_delta, C_0 = 0, cost_c = 20, cost_P = 0, cost_A = 13.60, B = budget)$k) %>%
  mutate(power_no_compositing = get_power_two_sample(n_1 = sample_size, k_1 = sample_size, n_2 = sample_size, k_2 = sample_size, mu_1 = mu_0, sigma_p_1 = sigma_p, mu_2 = mu_0 + delta, sigma_p_2 = sigma_p, sigma_delta = sigma_delta)) %>%
  mutate(power_full_compositing = get_power_two_sample(n_1 = max_sample_size, k_1 = 1, n_2 = max_sample_size, k_2 = 1, mu_1 = mu_0, sigma_p_1 = sigma_p, mu_2 = mu_0 + delta, sigma_p_2 = sigma_p, sigma_delta = sigma_delta)) %>%
  mutate(power_optimal_compositing = get_power_two_sample(n_1 = opt_n, k_1 = opt_k, n_2 = opt_n, k_2 = opt_k, mu_1 = mu_0, sigma_p_1 = sigma_p, mu_2 = mu_0 + delta, sigma_p_2 = sigma_p, sigma_delta = sigma_delta)) %>%
  pivot_longer(cols = starts_with("power_"), names_to = "Compositing", values_to = "power", names_prefix = "power_") %>%
  mutate(Compositing = ifelse(Compositing == "full_compositing", "Full", ifelse(Compositing == "no_compositing", "None", "Optimal"))) 

ggplot(power_change_topsoil, aes(x = delta, y = power, color = Compositing)) +
  geom_line(size = 1.5) +
  facet_grid(sample_size ~ land_use + Machine) +
  xlab("%TC Change") +
  ylab("Power of two-sample t-test") +
  scale_y_continuous(labels = scales::percent, limits = c(0,1), breaks = c(0,.25,.5,.75,1)) +
  xlim(0,1.5) +
  theme_bw() +
  theme(text = element_text(size = 16))

########### power of a permutation test to detect topsoil change #########
  
effect_grid = seq(0,2,by=.1)
run_twosample_sims <- function(sample_size, n_sims = 400){
 shift <- effect_grid
 perm_p_values <- matrix(NA, nrow = n_sims, ncol = length(shift))
 normal_p_values <- matrix(NA, nrow = n_sims, ncol = length(shift))
 for(i in 1:n_sims){
   for(j in 1:length(shift)){
     sample_1 <- sample(topsoil_TOC_rangeland, size = sample_size, replace = TRUE)
     sample_2 <- sample(topsoil_TOC_rangeland + shift[j], size = sample_size, replace = TRUE)
     diff_mean <- mean(sample_1) - mean(sample_2)
     normal_p_values[i,j] <- t.test(x = sample_1, y = sample_2, alternative = "two.sided")$p.value
     perm_p_values[i,j] <- t2p(diff_mean, two_sample(x = sample_1, y = sample_2, reps = 500), alternative = "two-sided")
   }
 }
 normal_power_shift <- colMeans(normal_p_values < .05)
 perm_power_shift <- colMeans(perm_p_values < .05)
 cbind("normal" = normal_power_shift, "permutation" = perm_power_shift)
}

#these take a while to run, we can save them as an object
#power_5 <- run_twosample_sims(sample_size = 5)
#power_10 <- run_twosample_sims(sample_size = 10)
#power_30 <- run_twosample_sims(sample_size = 30)
#power_90 <- run_twosample_sims(sample_size = 90)
# power_frame <- bind_rows(
#   as.data.frame(power_5) %>% mutate(sample_size = 5, effect_size = effect_grid),
#   as.data.frame(power_10) %>% mutate(sample_size = 10, effect_size = effect_grid),
#   as.data.frame(power_30) %>% mutate(sample_size = 30, effect_size = effect_grid),
#   as.data.frame(power_90) %>% mutate(sample_size = 90, effect_size = effect_grid)) %>%
#   rename("t test" = normal, "Permutation test" = permutation) %>%
#   pivot_longer(cols = c("t test", "Permutation test"), names_to = "Test", values_to = "Power")
#save(power_frame, file = "power_frame")
load("power_frame")

ggplot(power_frame, aes(x = effect_size, y = Power, color = Test, linetype = Test)) +
  geom_line(size = 1.5) +
  facet_wrap(~ sample_size) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent) +
  xlab("Effect size (additional % TOC)") +
  theme(text = element_text(size = 16))

############## analyze advantages of stratified versus uniform random sampling ###########
#we will empirically investigate the advantages of stratified sampling in these settings
#the idea is similar to a bootstrap. 
#take the samples as a population and then simulate taking subsamples

#first we analyze stratified sampling at a small scale by using the rangeland samples
#the strata are defined by transects, which roughly correspond to different locations on the ranch
rangeland_data_topsoil <- rangeland_master %>% 
  filter(depth == "a") %>%
  select(transect, sample_number, TOC) %>%
  filter(!is.na(TOC)) %>%
  arrange(transect)



#the sample size for simulations
n <- 90
#function to return the estimated mean and standard error given a population and sample index
get_mean_se <- function(population, sample_index){
  c(mean(population[sample_index]), sd(population[sample_index])/sqrt(length(sample_index)))
}
#function to return estimated mean and standard error given a population, sampling index, and stratification information
get_mean_se_stratified <- function(sample, N_strata){
  N <- sum(N_strata)
  strata_weights <- N_strata / N
  n_strata <- as.numeric(table(sample$transect))
  strata_means <- tapply(sample$TOC, sample$transect, mean)
  strata_vars <- tapply(sample$TOC, sample$transect, var)
  var_estimate <- N^(-2) * sum(N_strata^2 * strata_vars / n_strata)
  c(sum(strata_weights * strata_means), sqrt(var_estimate))
}

N_strata <- as.numeric(table(rangeland_data_topsoil$transect))
#helper function to make integer sample sizes with the sum preserved
round_strata_sizes <- function(n_strata){
  rounded_n_strata <- floor(n_strata)
  indices <- tail(order(n_strata-rounded_n_strata), round(sum(n_strata)) - sum(rounded_n_strata))
  rounded_n_strata[indices] <- rounded_n_strata[indices] + 1
  rounded_n_strata
}

#proportional allocation to strata
strata_weights_prop <- N_strata / length(rangeland_data_topsoil$transect)
n_strata_prop <- round_strata_sizes(n * strata_weights_prop)
sigma_strata <- tapply(rangeland_data_topsoil$TOC, rangeland_data_topsoil$transect, sd)
strata_weights_opt <- N_strata * sigma_strata / sum(N_strata * sigma_strata) 
n_strata_opt <- round_strata_sizes(n * strata_weights_opt)




  
#function to run a single simulation on the rangeland data
run_rangeland_sim <- function(data_frame){
  proportional_stratified_sample <- strata(data = data_frame, stratanames = "transect", size = n_strata_prop, method = "srswr")
  optimal_stratified_sample <- strata(data = data_frame, stratanames = "transect", size = n_strata_opt, method = "srswr")
  
  random_sample <- sample(1:nrow(data_frame), size = n, replace = TRUE)
  
  prop_stratified_estimates <- get_mean_se_stratified(sample = getdata(data_frame, proportional_stratified_sample), N_strata = table(data_frame$transect))
  opt_stratified_estimates <- get_mean_se_stratified(sample = getdata(data_frame, optimal_stratified_sample), N_strata = table(data_frame$transect))
  
  random_estimates <- get_mean_se(data_frame$TOC, random_sample)
  #local_pivotal_estimates <- get_mean_se(population, local_pivotal_sample)
  
  cbind(random_estimates, prop_stratified_estimates, opt_stratified_estimates)
}
#compute the estimand, i.e. the true mean in the rangeland topsoil data
true_rangeland_mean <- mean(rangeland_data_topsoil$TOC)

#run simulations, replicated 1000 times
rangeland_sims <- replicate(n = 2000, run_rangeland_sim(data_frame = rangeland_data_topsoil))

#compute properties of the samples
#the empirical coverage of 95% normal theory confidence intervals (should be 95%)
coverage_rangeland <- apply(rangeland_sims[1,,] - qnorm(p = .975) * rangeland_sims[2,,] <= true_rangeland_mean & true_rangeland_mean <= rangeland_sims[1,,] + qnorm(p = .975) * rangeland_sims[2,,], 1, mean)
#the width of a confidence interval, basically 4 times the SE
ci_width_rangeland <- apply(2 * qnorm(p = .975) * rangeland_sims[2,,], 1, mean)
#the actual standard error over simulations
se_rangeland <- apply(rangeland_sims, c(1,2), sd)[1,]
#the average estimated standard error
se_hat_rangeland <- apply(rangeland_sims, c(1,2), mean)[2,]
#the mean squared error
rmse_rangeland <- sqrt(apply((rangeland_sims - true_rangeland_mean)^2, c(1,2), mean)[1,])
mad_rangeland <- apply(abs(rangeland_sims - true_rangeland_mean), c(1,2), median)[1,]
#ratio of mse compared to uniform independent random sampling
#rmse_rangeland <- rmse_rangeland / rmse_rangeland[2]
#compile results
rangeland_results_frame <- rbind(coverage_rangeland, ci_width_rangeland, rmse_rangeland, mad_rangeland) %>%
  as_tibble() %>%
  mutate(property = c("coverage", "ci width", "RMSE", "MAD")) %>%
  select(property, random_estimates, prop_stratified_estimates, opt_stratified_estimates)



