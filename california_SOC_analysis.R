library(tidyverse)
library(readxl)
library(permuter)
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
null_distributions <- matrix(0, ncol = num_strata, nrow = B)
p_values <- rep(0, num_strata)
#iterate across strata (samples)
for(k in 1:num_strata){
  solitoc_TC <- reps_long$TC[reps_long$machine == "solitoc" & strata == levels(strata)[k]]
  costech_TC <- reps_long$TC[reps_long$machine == "costech" & strata == levels(strata)[k]]
  diff_means[k] <- mean(solitoc_TC) - mean(costech_TC)
  null_distributions[,k] <- two_sample(x = solitoc_TC, y = costech_TC, reps = B)
  p_values[k] <- mean(abs(diff_means[k]) <= abs(null_distributions[,k])) 
}
combined_p_value <- npc(diff_means, distr = null_distributions, combine=  "fisher", alternatives = "two-sided")


#average inorganic carbon is interesting as a predictor of difference between costech and solitoc
avg_TIC <- bind_rows(paicines_solitoc_reps, cropland_solitoc_reps) %>%
  mutate(sample_ID = paste(site, sample_number, sep = "_")) %>%
  group_by(sample_ID) %>%
  summarize(avg_TIC = mean(TIC), site = first(site))

p_value_frame <- data.frame(sample_ID = levels(strata), diff_mean = diff_means, p_value = p_values) %>%
  left_join(avg_TIC, by = "sample_ID") %>%
  as_tibble()

ggplot(p_value_frame, aes(x = log10(avg_TIC), y = diff_mean, shape = p_value < .05, color = site)) +
  geom_hline(yintercept = 0) +
  geom_point()

############### costs and optimal compositing #############
#precision of the sample mean 
precision_paicines <- paicines_master %>%
  group_by(depth) %>%
  summarize(mu = mean(TC, na.rm = T), sigma_p = sd(TC, na.rm = T)) %>%
  mutate(se_no_compositing = sqrt(get_variance(n = 60, k = 60, sigma_p = sigma_p, mu = mu, sigma_delta = sigma_delta_costech))) %>%
  mutate(se_full_compositing = sqrt(get_variance(n = 60, k = 1, sigma_p = sigma_p, mu = mu, sigma_delta = sigma_delta_costech)))

#power to detect a half-percentage point change in SOC
power_change_paicines <- paicines_master %>%
  group_by(depth) %>%
  summarize(mu = mean(TC, na.rm = T), sigma_p = sd(TC, na.rm = T)) %>%
  mutate(power_no_compositing = get_power_two_sample(n_1 = 60, k_1 = 60, n_2 = 60, k_2 = 60, mu_1 = mu, sigma_p_1 = sigma_p, mu_2 = mu + .5, sigma_p_2 = sigma_p, sigma_delta = sigma_delta_costech)) %>%
  mutate(power_full_compositing = get_power_two_sample(n_1 = 30, k_1 = 1, n_2 = 30, k_2 = 1, mu_1 = mu, sigma_p_1 = sigma_p, mu_2 = mu + .5, sigma_p_2 = sigma_p, sigma_delta = sigma_delta_costech))
  

