library(tidyverse)
library(readxl)
#permuter is used for permutation tests. It can be installed from Github: https://github.com/statlab/permuter
library(permuter)
library(sampling)
library(geoR)

source("functions.R")

########### read in data from excel spreadsheet(s) ##############
#note that some of the conversions to numeric cause NA warnings, I have spot checked these against the excel files, there does not appear to be a real issue.
#note, what we refer to as ECS 4010 in the paper is referred to as Costech in most of this code.
rangeland_master <- read_excel("R_Heterogeneity_Master_PS_04282021.xlsx", sheet = "Rangeland_All_soliTOC") %>%
  mutate(TOC = ifelse(TOC == "NA", NA, TOC)) %>%
  mutate(TIC = ifelse(TIC == "NA", NA, TIC)) %>%
  mutate(TC = ifelse(TC == "NA", NA, TC)) %>%
  mutate(TOC = as.numeric(TOC), TIC = as.numeric(TIC), TC = as.numeric(TC), sample_number = as.numeric(sample_number)) %>%
  mutate(depth_long = dplyr::recode(depth, a = "0-10 cm", b = "10-30 cm", c = "30-50 cm", d = "50-75 cm", e = "75-100 cm"))

rangeland_solitoc_reps <- read_excel("R_Heterogeneity_Master_PS_04282021.xlsx", sheet = "Rangeland_Reps_soliTOC") %>%
  mutate(TOC = ifelse(TOC == "NA", NA, TOC)) %>%
  mutate(TOC = as.numeric(TOC), sample_number = as.numeric(sample_number)) %>%
  mutate(machine = "solitoc")

rangeland_costech_reps <- read_excel("R_Heterogeneity_Master_PS_04282021.xlsx", sheet = "Rangeland_Reps_Costech") %>%
  mutate(TC = ifelse(TC == "NA", NA, TC), TOC = ifelse(TOC == "NA", NA, TOC)) %>%
  mutate(TC = as.numeric(TC), TOC = as.numeric(TOC), sample_number = as.numeric(sample_number)) %>%
  mutate(machine = "costech")

rangeland_BD <- read_excel("R_Heterogeneity_Master_PS_04282021.xlsx", sheet = "Rangeland_BD") %>%
  rename('bd' = 'bulk_density_g/cm3') %>%
  mutate(bd = as.numeric(bd)) %>%
  mutate(depth_long = dplyr::recode(depth, a = "0-10 cm", b = "10-30 cm", c = "30-50 cm", d = "50-75 cm", e = "75-100 cm"))

#data from various croplands around California
cropland_master <- read_excel("R_Heterogeneity_Master_PS_04282021.xlsx", sheet = "Cropland_All_Costech") %>%
  mutate(depth_long = dplyr::recode(depth, a = "0-15 cm", b = "15-30 cm", c = "30-60 cm", d = "60-100 cm"))

cropland_solitoc_reps <- read_excel("R_Heterogeneity_Master_PS_04282021.xlsx", sheet = "Cropland_Reps_soliTOC") %>%
  mutate(sample_number = as.numeric(sample_number)) %>%
  mutate(machine = "solitoc")

cropland_costech_reps <- read_excel("R_Heterogeneity_Master_PS_04282021.xlsx", sheet = "Cropland_Reps_Costech") %>%
  mutate(TC = ifelse(TC == "NA", NA, TC), TOC = ifelse(TOC == "NA", NA, TOC)) %>%
  mutate(TC = as.numeric(TC), TOC = as.numeric(TOC), sample_number = as.numeric(sample_number)) %>%
  mutate(machine = "costech")

cropland_BD <- read_excel("R_Heterogeneity_Master_PS_04282021.xlsx", sheet = "Cropland_BD")  %>%
  rename('bd' = 'bulk_density_g/cm3') %>%
  mutate(bd = as.numeric(bd))  %>%
  mutate(depth_long = dplyr::recode(depth, a = "0-15 cm", b = "15-30 cm", c = "30-60 cm", d = "60-100 cm"))

replicates_comparison <- read_excel("R_Heterogeneity_Master_PS_04282021.xlsx", sheet = "Replicate_Comparison") %>% 
  rename(TC_solitoc = TC_soliTOC, TOC_solitoc = TOC_soliTOC) %>%
  mutate_at(vars(starts_with(c("TC_", "TOC_"))), as.numeric)

standards_comparison <- read_excel("R_Heterogeneity_Master_PS_04282021.xlsx", sheet = "Standards_Comparison") %>%
  rename(TC = `TC%`, N = `N%`)


combined_master <- rangeland_master %>% 
  select(site, soil_type, depth, depth_long, TC) %>%
  bind_rows(cropland_master %>% select(site, soil_type, depth, depth_long, TC)) %>%
  mutate(land_use = ifelse(site == "RANGE", "Rangeland", "Cropland"))

###################### TOC concentration and BD in space ######################

#stacked TC histograms for rangeland data
ggplot(rangeland_master, aes(TC)) +
  geom_histogram(binwidth = 0.25, alpha = 1) +
  facet_grid(depth_long ~ transect) +
  xlab("% Total Carbon (TC)") +
  ylab("Number of Samples") +
  #ylim(0,23) +
  theme_bw() +
  theme(text = element_text(size = 20))


#violin plot for rangeland data
ggplot(rangeland_master, aes(y = TC, x = transect, fill = transect, color = transect)) +
  geom_violin() +
  coord_flip() +
  facet_grid(depth_long ~ .) +
  scale_fill_manual(values = c("firebrick","forestgreen","darkorchid4","darkorange3","steelblue")) +
  scale_color_manual(values = c("firebrick","forestgreen","darkorchid4","darkorange3", "steelblue")) +
  xlab("Transect") +
  ylab("% Total Carbon (TC)") +
  theme_bw() +
  theme(text = element_text(size = 20), legend.position = "none")


#TOC histogram for rangeland
ggplot(rangeland_master, aes(TOC)) +
  geom_histogram(binwidth = 0.25) +
  facet_grid(depth_long ~ transect) +
  xlab("% Total Organic Carbon (TOC)") +
  ylab("Number of Samples") +
  theme_bw() +
  theme(text = element_text(size = 20))


#histograms for cropland data
ggplot(cropland_master, aes(TC)) +
  geom_histogram(binwidth = 0.25) +
  facet_grid(depth_long ~ site) +
  xlab("% Total Carbon (TC)") +
  ylab("Number of Samples") +
  theme_bw() +
  #xlim(0,8) +
  #ylim(0,23) +
  theme(text = element_text(size = 20))

#violin plot for cropland data
ggplot(cropland_master, aes(y = TC, x = site, fill = site, color = site)) +
  geom_violin() +
  coord_flip() +
  facet_grid(depth ~ .) +
  scale_fill_manual(values = c("hotpink","gold3","firebrick","forestgreen","darkorchid4","darkorange3","steelblue")) +
  scale_color_manual(values = c("hotpink","gold3","firebrick","forestgreen","darkorchid4","darkorange3", "steelblue")) +
  xlab("Transect") +
  ylab("% Total Carbon (TC)") +
  theme_bw() +
  theme(text = element_text(size = 20), legend.position = "none")




#summary tables
rangeland_summary <- rangeland_master %>% 
  group_by(depth, transect) %>%
  summarize(mean_TC = mean(TC, na.rm = T), cv_TC = sd(TC, na.rm = T) / mean(TC, na.rm = T)) %>%
  mutate(cv_TC = paste("(", round(cv_TC,2), ")", sep = "")) %>%
  mutate(mean_TC = round(mean_TC, 2)) %>%
  unite(col = "mean (SD)", mean_TC, cv_TC, sep = " ") %>%
  pivot_wider(names_from = transect, values_from = "mean (SD)")

#rangeland data combined across transects
rangeland_summary_overall <- rangeland_master %>% 
  group_by(depth) %>%
  summarize(mean_TC = mean(TC, na.rm = T), sd_TC = sd(TC, na.rm = T))

#note cropland is total carbon not total organic carbon
cropland_summary <- cropland_master %>% 
  group_by(depth, site) %>%
  summarize(mean_TC = mean(TC, na.rm = T), sd_TC = sd(TC, na.rm = T))
cropland_summary_overall <- cropland_summary %>%
  group_by(depth) %>%
  summarize(mean_TC = mean(mean_TC, na.rm = T), sd_TC = mean(sd_TC, na.rm = T))

#average and CV total carbon by depth, site for table 1
table_1 <- combined_master %>%
  group_by(site, depth, land_use) %>%
  summarize(mean_TC = mean(TC, na.rm = T), cv_TC = sd(TC, na.rm = T) / mean(TC, na.rm = T)) %>% 
  group_by(depth, land_use) %>%
  summarize(mean_TC = mean(mean_TC, na.rm = T), cv_TC = mean(cv_TC, na.rm = T)) %>%
  mutate(cv_TC = paste("(", round(cv_TC,2), ")", sep = "")) %>%
  mutate(mean_TC = round(mean_TC, 2)) %>%
  unite(col = "mean (CV)", mean_TC, cv_TC, sep = " ") %>%
  pivot_wider(names_from = land_use, values_from = "mean (CV)")



#permutation ANOVAs 
#need to drop missing TC observations and only use complete data
rangeland_master_complete <- rangeland_master %>% 
  filter(!is.na(TC)) %>%
  mutate(transect_id = as.numeric(as.factor(transect))) %>%
  mutate(depth_id = as.numeric(as.factor(depth)))
cropland_master_complete <- cropland_master %>% 
  filter(!is.na(TC)) %>%
  mutate(site_id = as.numeric(as.factor(site))) %>%
  mutate(depth_id = as.numeric(as.factor(depth)))

#Rangeland permutation ANOVA results
rangeland_transect_ANOVA <- t2p(
  tst = sum(table(rangeland_master_complete$transect_id) * tapply(rangeland_master_complete$TC, rangeland_master_complete$transect_id, mean)^2),
  distr = k_sample(x = rangeland_master_complete$TC, group = rangeland_master_complete$transect_id, stat = "oneway_anova", reps = 10000),
  alternative = "two-sided")
rangeland_depth_ANOVA <- t2p(
  tst = sum(table(rangeland_master_complete$depth_id) * tapply(rangeland_master_complete$TC, rangeland_master_complete$depth_id, mean)^2),
  distr = k_sample(x = rangeland_master_complete$TC, group = rangeland_master_complete$depth_id, stat = "oneway_anova", reps = 10000),
  alternative = "two-sided")

#Cropland permutation ANOVA results
cropland_site_ANOVA <- t2p(
  tst = sum(table(cropland_master_complete$site_id) * tapply(cropland_master_complete$TC, cropland_master_complete$site_id, mean)^2),
  distr = k_sample(x = cropland_master_complete$TC, group = cropland_master_complete$site_id, stat = "oneway_anova", reps = 10000),
  alternative = "two-sided")
cropland_depth_ANOVA <- t2p(
  tst = sum(table(cropland_master_complete$depth_id) * tapply(cropland_master_complete$TC, cropland_master_complete$depth_id, mean)^2),
  distr = k_sample(x = cropland_master_complete$TC, group = cropland_master_complete$depth_id, stat = "oneway_anova", reps = 10000),
  alternative = "two-sided")



#bulk densities
#NOTE: histograms are not disaggregated by location, which drives variation
ggplot(rangeland_BD, aes(bd)) +
  geom_histogram(binwidth = 0.05) +
  facet_grid(depth_long ~ .) +
  xlab("Bulk Density") +
  ylab("Number of samples") +
  coord_cartesian(xlim=c(0.8,1.8)) +
  scale_x_continuous(breaks = c(0.8, 1.0, 1.2, 1.4, 1.6, 1.8)) +
  #scale_y_continuous(breaks = c(0,2,4,6), limits = c(0,6)) +
  theme_bw() +
  theme(text = element_text(size = 20), panel.grid.minor = element_blank())

#cropland bulk density for site 7
ggplot(cropland_BD %>% filter(site == "CROP7"), aes(bd)) +
  geom_histogram(binwidth = 0.05) +
  facet_grid(depth_long ~ .) +
  xlab("Bulk Density") +
  ylab("Number of samples") +
  coord_cartesian(xlim=c(0.8,1.8)) +
  scale_x_continuous(breaks = c(0.8, 1.0, 1.2, 1.4, 1.6, 1.8)) +
  scale_y_continuous(breaks = c(0,2,4,6), limits = c(0,6)) +
  theme_bw() +
  theme(text = element_text(size = 20), panel.grid.minor = element_blank())


#all cropland bulk density
ggplot(cropland_BD, aes(bd)) +
  geom_histogram(binwidth = 0.05) +
  facet_grid(depth ~ .) +
  xlab("Bulk Density") +
  ylab("Number of Samples") +
  scale_x_continuous(breaks = c(0.8, 1.2, 1.6, 2.0, 2.4), limits = c(0.8,2.4)) +
  scale_y_continuous(breaks = c(0,2,4,6,8,10), limits = c(0,10)) +
  theme_bw() +
  theme(text = element_text(size = 20))

rangeland_summary_bd <- rangeland_BD %>%
  group_by(depth) %>%
  summarize(mean_bd = mean(bd, na.rm = TRUE), sd_bd = sd(bd, na.rm = TRUE)) %>%
  mutate(cv_bd = sd_bd / mean_bd)

#averages within sites
cropland_summary_bd <- cropland_BD %>%
  group_by(depth, site) %>%
  summarize(mean_bd = mean(bd, na.rm = TRUE), sd_bd = sd(bd, na.rm = TRUE)) %>%
  mutate(cv_bd = sd_bd / mean_bd)

cropland_site7_bd <- cropland_summary_bd %>% 
  filter(site == "CROP7")

average_cropland_bd <- cropland_summary_bd %>% 
  group_by(depth) %>% 
  summarize(mean_bd = mean(mean_bd, na.rm = TRUE), sd_bd = sd(sd_bd, na.rm = TRUE), cv_bd = mean(cv_bd, na.rm = TRUE))


#stock totals on rangeland 
#first compute volumetric concentration (vc) at each sampled point
#then compute rangeland stocks within depth
#assume that the entire volume is soil 
rangeland_stocks <- rangeland_master %>% 
  select(transect, depth, TC) %>%
  mutate(TC = TC / 100) %>% #since TC is currently in percent
  left_join(
    rangeland_BD %>% 
      group_by(transect, depth) %>% 
      summarize(mean_BD = mean(bd, na.rm = TRUE), max_BD = max(bd, na.rm = TRUE), min_BD = min(bd, na.rm = TRUE)) %>% 
      mutate(max_BD = ifelse(max_BD == -Inf, NA, max_BD), min_BD = ifelse(min_BD == Inf, NA, min_BD))
  ) %>%
  mutate(mean_stock = TC * mean_BD, max_stock = TC * max_BD, min_stock = TC * min_BD) %>%
  mutate(length = as.numeric(recode(depth, a = "10", b = "20", c = "20", d = "25", e = "25"))) %>%
  group_by(depth) %>%
  summarize(mean_TC = mean(TC, na.rm = T), mean_stock_g_cm3 = mean(mean_stock, na.rm = T), se_stock_g_cm3 = sd(mean_stock, na.rm = T)/sqrt(n()), max_stock_g_cm3 = mean(max_stock, na.rm = T), min_stock_g_cm3 = mean(min_stock, na.rm = T), length = first(length)) %>%
  mutate(mean_stock_Mg_ha = mean_stock_g_cm3 * length * 1e8 * 1e-6) %>% #1e8 cm2 in a hectare, 1e-6 metric tons in a gram
  mutate(se_stock_Mg_ha = se_stock_g_cm3 * length * 1e8 * 1e-6) %>%
  mutate(max_stock_Mg_ha = max_stock_g_cm3 * length * 1e8 * 1e-6) %>%
  mutate(min_stock_Mg_ha = min_stock_g_cm3 * length * 1e8 * 1e-6) 

#whole profile stock on rangeland
rangeland_wp_stock <- rangeland_stocks %>%
  summarize(wp_stock_Mg_ha = sum(mean_stock_Mg_ha), wp_stock_se_Mg_ha = sum(se_stock_Mg_ha), wp_stock_max_BD = sum(max_stock_Mg_ha), wp_stock_min_BD = sum(min_stock_Mg_ha))

#cropland stocks
cropland_stocks <- cropland_master %>% 
  select(site, depth, TC) %>%
  mutate(TC = TC / 100) %>% #since TC is currently in percent
  left_join(
    cropland_BD %>% 
      group_by(site, depth) %>% 
      summarize(mean_BD = mean(bd, na.rm = TRUE), max_BD = max(bd, na.rm = TRUE), min_BD = min(bd, na.rm = TRUE)) %>% 
      mutate(max_BD = ifelse(max_BD == -Inf, NA, max_BD), min_BD = ifelse(min_BD == Inf, NA, min_BD))
  ) %>%
  mutate(mean_stock = TC * mean_BD, max_stock = TC * max_BD, min_stock = TC * min_BD) %>%
  mutate(length = as.numeric(recode(depth, a = "15", b = "15", c = "30", d = "40"))) %>%
  group_by(site, depth) %>%
  summarize(mean_TC = mean(TC, na.rm = T), mean_stock_g_cm3 = mean(mean_stock, na.rm = T), se_stock_g_cm3 = sd(mean_stock, na.rm = T)/sqrt(n()), max_stock_g_cm3 = mean(max_stock, na.rm = T), min_stock_g_cm3 = mean(min_stock, na.rm = T), length = first(length)) %>%
  mutate(mean_stock_Mg_ha = mean_stock_g_cm3 * length * 1e8 * 1e-6) %>% #1e8 cm2 in a hectare, 1e-6 metric tons in a gram
  mutate(se_stock_Mg_ha = se_stock_g_cm3 * length * 1e8 * 1e-6) %>%
  mutate(max_stock_Mg_ha = max_stock_g_cm3 * length * 1e8 * 1e-6) %>%
  mutate(min_stock_Mg_ha = min_stock_g_cm3 * length * 1e8 * 1e-6) 

#whole profile cropland stock
cropland_wp_stock <- cropland_stocks %>%
  group_by(site) %>%
  summarize(wp_stock_Mg_ha = sum(mean_stock_Mg_ha), wp_stock_se_Mg_ha = sum(se_stock_Mg_ha), wp_stock_max_BD = sum(max_stock_Mg_ha), wp_stock_min_BD = sum(min_stock_Mg_ha))



############## replicates and assay error ##############
#compare variation due to sampling and variation due to assay
#note that samples were specifically selected along a grid of carbon concentrations
#this means that sampling heterogeneity is likely to be artificially inflated in the replicate data.
# replicates_long <- rangeland_solitoc_reps %>% 
#   bind_rows(cropland_solitoc_reps) %>%
#   bind_rows(rangeland_costech_reps) %>%
#   bind_rows(cropland_costech_reps)
# average_error <- replicates_long %>% 
#   group_by(machine, sample_number, site) %>%
#   summarize(error_estimate = var(TC, na.rm = T) / (mean(TC,  na.rm = T)^2 - var(TC, na.rm = T)/n())) %>%
#   group_by(machine) %>%
#   summarize(within_error = mean(error_estimate, na.rm = T), samples = n()) %>% 
#   pivot_wider(names_from = machine, values_from = within_error)



#compute percent error on each replicated sample
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
  #geom_density(alpha = .75) +
  geom_histogram(binwidth = 1, alpha = .75, position = "stack") +
  xlab("Relative Error") +
  ylab("Number of Samples") +
  theme_bw() +
  scale_x_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_manual(name = "Instrument", values = c("steelblue","darkorange3"), labels = c("ECS 4010", "SoliTOC")) +
  theme(text = element_text(size = 20))

median_sigma_delta <- assay_error_long %>%
  group_by(machine) %>%
  summarize(sigma_delta = median(sigma_delta_TC))

mean_sigma_delta <- assay_error_long %>%
  group_by(machine) %>%
  summarize(sigma_delta = mean(sigma_delta_TC))

upper_quartile_sigma_delta <- assay_error_long %>%
  group_by(machine) %>%
  summarize(sigma_delta = quantile(sigma_delta_TC, .75))

#proportions of assay heterogeneity as proportions of field heterogeneity
sample_size <- 90
assay_field_proportions <- combined_master %>%
  group_by(site, depth_long, land_use) %>%
  summarize(mean_TC = mean(TC, na.rm = T), sd_TC = sd(TC, na.rm = T)) %>% 
  group_by(depth_long, land_use) %>%
  summarize(sd_TC = mean(sd_TC, na.rm = T), mean_TC = mean(mean_TC, na.rm = T)) %>%
  mutate(SoliTOC = median_sigma_delta$sigma_delta[median_sigma_delta$machine == "solitoc"]) %>%
  mutate(Costech = median_sigma_delta$sigma_delta[median_sigma_delta$machine == "costech"]) %>%
  pivot_longer(cols = c(SoliTOC, Costech), names_to = "machine", values_to = "assay_error") %>%
  mutate(assay_contribution = (assay_error * mean_TC)^2) %>%
  mutate(variance_prop_nocompositing = assay_contribution / (sd_TC^2 + assay_contribution)) %>%
  mutate(variance_prop_fullcompositing =  sample_size * assay_contribution / (sd_TC^2 + sample_size * assay_contribution)) %>%
  pivot_longer(cols = c(variance_prop_nocompositing, variance_prop_fullcompositing), names_to = "compositing", names_prefix = "variance_prop_", values_to = "assay_proportion") %>%
  mutate(compositing = ifelse(compositing == "nocompositing", "No Compositing", "Full Compositing")) %>%
  mutate(compositing = factor(compositing, levels = c("No Compositing", "Full Compositing"))) %>%
  mutate(field_proportion = 1 - assay_proportion) %>%
  pivot_longer(cols = c(assay_proportion,field_proportion), names_to = "source", values_to = "proportion") %>%
  mutate(source = ifelse(source == "assay_proportion", "Analytical Variability", "Spatial Heterogeneity")) %>%
  mutate(land_use = factor(land_use, levels = c("Rangeland","Cropland"))) %>%
  mutate(machine = ifelse(machine == "Costech", "ECS 4010", machine))

#variance proportions plot
ggplot(assay_field_proportions, aes(x = depth_long, fill = source, y = proportion)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_grid(machine ~ compositing + land_use, scales = "free_x") +
  ylim(0,1) +
  guides(fill = guide_legend(title = "Source of Variance")) +
  scale_fill_manual(values = c("darkorange3","steelblue")) +
  ylab("Approximate Proportion of Estimation Variance") +
  xlab("Depth") +
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0))



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
#all samples that were not significantly different (TC and TOC)
assay_histogram_plot_accepted <- ggplot(reps_long_both %>% filter(rejected == 0), aes(carbon, fill = machine)) +
  geom_histogram( alpha = .75, position = "stack", binwidth = .2) +
  facet_wrap(~ identifier) +
  #xlim(0,6) +
  xlab("TOC% or TC%") +
  ylab("Number of Replicates") +
  ggtitle("Similar Assays") +
  scale_fill_manual(name = "Instrument", values = c("steelblue","darkorange3"), labels = c("ECS 4010", "SoliTOC")) +
  scale_y_continuous(breaks=seq(0,9)) +
  scale_x_continuous(breaks=seq(0,6), limits = c(0,6)) +
  theme_bw() +
  theme(text = element_text(size = 14), panel.grid.minor.y = element_blank())
#only TOC that were rejected (since rejected TC were re-run for TOC)
assay_histogram_plot_rejected <- ggplot(reps_long_both %>% filter(rejected == 1, carbon_type == "TOC"), aes(carbon, fill = machine)) +
  geom_histogram( alpha = .75, position = "stack", binwidth = .025) +
  facet_wrap(~ identifier) +
  #xlim(0,6) +
  xlab("TOC%") +
  ylab("Number of Replicates") +
  ggtitle("Significantly Different Analyses") +
  scale_fill_manual(name = "Instrument", values = c("steelblue","darkorange3"), labels = c("ECS 4010", "SoliTOC")) +
  scale_y_continuous(breaks=seq(0,9)) +
  scale_x_continuous(breaks=seq(0,3), limits = c(0,3.2)) +
  theme_bw() +
  theme(text = element_text(size = 20), panel.grid.minor.y = element_blank())


#standards comparison
known_standards <- standards_comparison %>%
  group_by(sample_ID) %>%
  summarize(known_TC = first(known_TC))
standards_comparison <- standards_comparison %>%
  mutate(known_TC = ifelse(sample_ID == "EML", 1.86, 0.926))

standards_histogram <- ggplot(standards_comparison, aes(TC, fill = machine)) +
  geom_histogram( alpha = .75, position = "stack", binwidth = .01) +
  geom_vline(data = standards_comparison %>% filter(sample_ID == "EML"), aes(xintercept = 1.86)) +
  geom_vline(data = standards_comparison %>% filter(sample_ID == "LECO"), aes(xintercept = 0.926)) +
  scale_fill_manual(name = "Instrument", values = c("steelblue","darkorange3"), labels = c("ECS 4010", "SoliTOC")) +
  facet_grid(~ sample_ID, scales = "free") +
  xlim(0.8, 2.2) +
  scale_y_continuous(breaks=seq(0,12)) +
  xlab("TC%") +
  ylab("Number of Replicates") +
  theme_bw() +
  theme(text = element_text(size = 20), panel.grid.minor.y = element_blank())


############# spatial correlation of TC concentrations on rangeland ##########
#first just do for a single transect
T_topsoil_TC <- rangeland_master %>% 
  filter(transect == "T", depth == "a") %>%
  arrange(sample_number) %>%
  pull(TC)

#dataframe for use with geodata
T_topsoil_geodata <- as.geodata(data.frame(x = 1:length(T_topsoil_TC), y = 1, TC = T_topsoil_TC), coords.col = 1:2, data.col = 3)
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
  TC <- rangeland_master %>% 
    filter(transect == rangeland_transect, depth == rangeland_depth) %>%
    arrange(sample_number) %>%
    pull(TC)
  
  #dataframe for use with geodata
  geodataframe <- as.geodata(data.frame(x = 3 * 1:length(TC), y = 1, TC = TC), coords.col = 1:2, data.col = 3)
  variogram <- variog(geodataframe, option = "bin")
  if(plot){
    plot(variogram, xlab = "Distance (meters)", ylab = "Semivariance", main = paste("Variogram for transect", rangeland_transect, "and depth", rangeland_depth))
  } else{
    variogram
  }
}



variog_T <- plot_variogram(rangeland_depth = "a", rangeland_transect = "T", plot = TRUE)
variog_Mx <- plot_variogram(rangeland_depth = "a", rangeland_transect = "Mx", plot = TRUE)
variog_My <- plot_variogram(rangeland_depth = "a", rangeland_transect = "My", plot = TRUE)
variog_Bx <- plot_variogram(rangeland_depth = "a", rangeland_transect = "Bx", plot = TRUE)
variog_By <- plot_variogram(rangeland_depth = "a", rangeland_transect = "By", plot = TRUE)




############### costs and optimal compositing #############
#Costs based on 'Cost of Carbon Analysis.xlsx' assembled by Jessica Chiartas
#in sampling, about 20-24 samples (4 transects * 5-6 composites) costs 400$ at 25$/hr and 16hrs, so we'll take 20 USD per sample.
#the commercial labs had costs of 19.5, 7.75, 22, 10, and 8.75 USD per sample to assess TOC. So we'll take the cost of measurement with dry combustion as 13.60, including sample prep.
#the in-house cost of measurement was 2.78


#sample size
sample_size <- 100
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
  mutate(se_no_compositing_costech = sqrt(get_variance(n = sample_size, k = sample_size, sigma_p = sigma_p, mu = mu, sigma_delta = sigma_delta_costech))) %>%
  mutate(se_full_compositing_costech = sqrt(get_variance(n = max_sample_size, k = 1, sigma_p = sigma_p, mu = mu, sigma_delta = sigma_delta_costech))) %>%
  mutate(se_optimal_compositing_costech = sqrt(get_minimum_error(sigma_p = sigma_p, mu = mu, sigma_delta = sigma_delta_costech, C_0 = 0, cost_c = 20, cost_P = 0, cost_A = 13.60, B = max_budget)$optimum_variance)) %>%
  mutate(se_no_compositing_solitoc = sqrt(get_variance(n = sample_size, k = sample_size, sigma_p = sigma_p, mu = mu, sigma_delta = sigma_delta_solitoc))) %>%
  mutate(se_full_compositing_solitoc = sqrt(get_variance(n = max_sample_size, k = 1, sigma_p = sigma_p, mu = mu, sigma_delta = sigma_delta_solitoc))) %>%
  mutate(se_optimal_compositing_solitoc = sqrt(get_minimum_error(sigma_p = sigma_p, mu = mu, sigma_delta = sigma_delta_solitoc, C_0 = 0, cost_c = 20, cost_P = 0, cost_A = 13.60, B = max_budget)$optimum_variance)) %>%
  arrange(land_use, depth)

precision_combined_long <- precision_combined %>%
  pivot_longer(cols = starts_with("se_"), names_to = "compositing", values_to = "std_error", names_prefix = "se_") %>%
  separate(compositing, into = c("Compositing", "Instrument"), sep = "_compositing_") %>%
  mutate(Instrument = ifelse(Instrument == "costech", "ECS 4010", "SoliTOC")) %>%
  mutate(Compositing = recode(Compositing, no = "None", full = "Full", optimal = "Optimal"))




#precision of sample mean in topsoil
mu_rangeland <- precision_combined %>% filter(depth == "a", land_use == "Rangeland") %>% pull(mu)
sigma_p_rangeland <- precision_combined %>% filter(depth == "a", land_use == "Rangeland") %>% pull(sigma_p)
mu_cropland <- precision_combined %>% filter(depth == "a", land_use == "Cropland") %>% pull(mu)
sigma_p_cropland <- precision_combined %>% filter(depth == "a", land_use == "Cropland") %>% pull(sigma_p)

rangeland_grid_precision <- expand.grid(
  land_use = "Rangeland",
  sample_size = seq(1,110,by=1),
  mu = mu_rangeland, 
  sigma_p = sigma_p_rangeland, 
  sigma_delta = c(sigma_delta_costech, sigma_delta_solitoc)
)
cropland_grid_precision <- expand.grid(
  land_use = "Cropland",
  sample_size = seq(1,110,by=1),
  mu = mu_cropland, 
  sigma_p = sigma_p_cropland, 
  sigma_delta = c(sigma_delta_costech, sigma_delta_solitoc)
)

precision_topsoil <- rangeland_grid_precision %>%
  bind_rows(cropland_grid_precision) %>%
  mutate(Machine = ifelse(sigma_delta == sigma_delta_costech, "ECS 4010", "SoliTOC")) %>%
  mutate(budget = 20 * sample_size + 13.6 * sample_size) %>%
  mutate(max_sample_size = floor((budget - 13.6) / 20)) %>%
  mutate(se_no_compositing = sqrt(get_variance(n = sample_size, k = sample_size, sigma_p = sigma_p, mu = mu, sigma_delta = sigma_delta))) %>%
  mutate(se_full_compositing = sqrt(get_variance(n = max_sample_size, k = 1, sigma_p = sigma_p, mu = mu, sigma_delta = sigma_delta))) %>%
  mutate(se_optimal_compositing = sqrt(get_minimum_error(sigma_p = sigma_p, mu = mu, sigma_delta = sigma_delta, C_0 = 0, cost_c = 20, cost_P = 0, cost_A = 13.60, B = budget)$optimum_variance)) %>%
  pivot_longer(cols = starts_with("se_"), names_to = "Compositing", values_to = "std_error", names_prefix = "se_") %>%
  mutate(Compositing = recode(Compositing, no_compositing = "None", full_compositing = "Full", optimal_compositing = "Optimal"))

precision_plot <- ggplot(precision_topsoil, aes(x = sample_size, y = std_error, color = Compositing)) +
  geom_line(size = 1.3) +
  facet_grid(land_use ~ Machine) + 
  theme_bw() +
  scale_color_manual(values = c("darkorange3", "forestgreen","steelblue")) +
  xlab("Sample Size") +
  ylab("Standard Error of Sample Mean") +
  coord_cartesian(xlim = c(5,100), ylim = c(0,.4)) +
  theme(text = element_text(size = 20))

#Power to detect change using t-test
rangeland_grid_power <- expand.grid(
  land_use = "Rangeland",
  sample_size = c(90),
  mu_0 = mu_rangeland,
  delta = seq(0, .5 * mu_rangeland, length.out = 30),
  sigma_p = sigma_p_rangeland, 
  sigma_delta = c(sigma_delta_costech, sigma_delta_solitoc)
)
cropland_grid_power <- expand.grid(
  land_use = "Cropland",
  sample_size = c(90),
  mu_0 = mu_cropland,
  delta = seq(0, .5 * mu_cropland, length.out = 100),
  sigma_p = sigma_p_cropland, 
  sigma_delta = c(sigma_delta_costech, sigma_delta_solitoc)
)

power_change_topsoil <- rangeland_grid_power %>%
  bind_rows(cropland_grid_power) %>%
  mutate(Machine = ifelse(sigma_delta == sigma_delta_costech, "ECS 4010", "SoliTOC")) %>%
  mutate(relative_delta = ifelse(land_use == "Cropland", delta / mu_cropland, delta / mu_rangeland)) %>%
  mutate(budget = 20 * sample_size + 13.6 * sample_size) %>%
  mutate(max_sample_size = floor((budget - 13.6) / 20)) %>%
  mutate(opt_n = get_minimum_error(sigma_p = sigma_p, mu = mu_0, sigma_delta = sigma_delta, C_0 = 0, cost_c = 20, cost_P = 0, cost_A = 13.60, B = budget)$n) %>%
  mutate(opt_k = get_minimum_error(sigma_p = sigma_p, mu = mu_0, sigma_delta = sigma_delta, C_0 = 0, cost_c = 20, cost_P = 0, cost_A = 13.60, B = budget)$k) %>%
  mutate(power_no_compositing = get_power_two_sample(n_1 = sample_size, k_1 = sample_size, n_2 = sample_size, k_2 = sample_size, mu_1 = mu_0, sigma_p_1 = sigma_p, mu_2 = mu_0 + delta, sigma_p_2 = sigma_p, sigma_delta = sigma_delta)) %>%
  mutate(power_full_compositing = get_power_two_sample(n_1 = max_sample_size, k_1 = 1, n_2 = max_sample_size, k_2 = 1, mu_1 = mu_0, sigma_p_1 = sigma_p, mu_2 = mu_0 + delta, sigma_p_2 = sigma_p, sigma_delta = sigma_delta)) %>%
  mutate(power_optimal_compositing = get_power_two_sample(n_1 = opt_n, k_1 = opt_k, n_2 = opt_n, k_2 = opt_k, mu_1 = mu_0, sigma_p_1 = sigma_p, mu_2 = mu_0 + delta, sigma_p_2 = sigma_p, sigma_delta = sigma_delta)) %>%
  pivot_longer(cols = starts_with("power_"), names_to = "Compositing", values_to = "power", names_prefix = "power_") %>%
  mutate(Compositing = ifelse(Compositing == "full_compositing", "Full", ifelse(Compositing == "no_compositing", "None", "Optimal"))) %>% 
  mutate(Budget = paste("$", formatC(budget, format = "d", big.mark = ","), sep = ""))

ggplot(power_change_topsoil, aes(x = 100*relative_delta, y = power, color = Compositing)) +
  geom_line(size = 1.2) +
  facet_grid(Machine ~ land_use) +
  xlab("Relative TC% Change") +
  ylab("Power of two-sample t-test") +
  scale_y_continuous(labels = scales::percent, limits = c(0,1), breaks = c(0,.25,.5,.75,1)) +
  scale_x_continuous(labels = function(x) paste0(x,"%")) +
  scale_color_manual(values = c("firebrick","forestgreen","steelblue")) +
  coord_cartesian(xlim = c(0,40)) +
  theme_bw() +
  theme(text = element_text(size = 20), axis.text.x = element_text(size = 16))




########### two-sample inference  #########
#an example where the t-test fails to give valid inference under the null (the mean does not change)
N <- 1000
#means are exactly 3 in both populations. Population 1 is highly skewed so that high values (which make means equal) are rarely sampled
pop_1 <- c(rnorm(N-10, mean = 3, sd = .05), rep(20, 10))
pop_1 <- pop_1 - mean(pop_1) + 3
pop_2 <- rep(3, N) + rnorm(N, sd = .05)
pop_2 <- pop_2 - mean(pop_2) + 3

par(mfrow = c(1,2))
hist(pop_1, breaks = 50, xlim = c(1,20), xlab = "%TC at time 1", freq = FALSE, main = "", cex.axis = 1.5, cex.lab = 1.5)
hist(pop_2, breaks = 50, xlim = c(1,20), xlab = "%TC at time 2", freq = FALSE, main = "", cex.axis = 1.5, cex.lab = 1.5)
par(mfrow = c(1,1))

#simulations are actually run in the script california_SOC_batch_simulations.R
par(mar = c(5.1,4.3, 4.3, 2.1))
plot(y = t_test_rejection_rate, x = n_grid, type ='l', ylim = c(0,1), xlab = "Sample size", ylab = "Simulated significance level", col = 'darkorange3', lwd = 4, cex.axis = 1.8, cex.lab = 1.8)
points(y = LMT_rejection_rate, x = n_grid, type = 'l', col = 'steelblue', lwd = 4)
#points(y = hedged_rejection_rate, x = n_grid, type = 'l', col = 'darkorange3', lwd = 2, lty = "dashed")
legend(x = 100, y = 0.8, legend = c("Nonparametric test","t-test"), lty = c( "solid","solid"), col = c("steelblue","darkorange3"), lwd = 4, bty = "n", cex = 1.5)
abline(a = 0.05, b = 0, lty = 'dashed', col = 'black', lwd = 2)
par(mar = c(5.1,4.1, 4.1, 2.1))


#plot simulations
load("validity_simulations")
validity_frame <- validity_simulations %>%
  reduce(bind_rows) %>%
  pivot_longer(cols = c("t_test_rejections", "LMT_rejections"), names_to = "test", values_to = "rejection_rate") %>%
  mutate(test = recode(test, "t_test_rejections" ="t-test", "LMT_rejections" = "Nonparametric test")) %>%
  mutate(population = recode(population,  "rangeland_to_cropland" = "R to C", "hotspots_to_deadspots" = "R to -C", "skewed" = "Major hotspot", "symmetric_reduced_spread" = "Gaussian unequal variance")) %>%
  group_by(population, test) %>%
  mutate(running_rejection_rate = cummin(rejection_rate)) %>%
  ungroup()

ggplot(validity_frame, aes(x = sample_size, y = running_rejection_rate, color = test, linetype = test)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = .05, linetype = "solid") +
  facet_wrap(~population) +
  coord_cartesian(xlim = c(5,100), ylim = c(0,.5)) +
  theme_bw() +
  scale_color_manual(values = c("darkorange3", "forestgreen")) +
  ylab("Simulated significance level") +
  xlab("Sample size") +
  theme(text = element_text(size = 20), axis.text.x = element_text(size = 16))



# power of t-test (unstratified and stratified) or nonparametric (unstratified) to detect topsoil change 
# run tests on topsoil from rangeland and cropland
topsoil_rangeland <- rangeland_master %>% 
  filter(depth == "a") %>%
  filter(!is.na(TC)) %>%
  arrange(transect) %>%
  select(transect, TC)
#cropland site 5 has the second most samples (20) and the lowest heterogeneity
topsoil_cropland <- cropland_master %>% 
  filter(depth == "a") %>%
  filter(site == "CROP5") %>%
  filter(!is.na(TC)) 

topsoil_TC_rangeland <- topsoil_rangeland %>% pull(TC)
topsoil_TC_cropland <- topsoil_cropland %>% pull(TC)

# set up stratification parameters
strata <- topsoil_rangeland$transect
N_strata <- as.numeric(table(topsoil_rangeland$transect))
strata_weights_prop <- N_strata / length(strata)
sigma_strata <- tapply(topsoil_TC_rangeland, strata, sd)
strata_weights_opt <- N_strata * sigma_strata / sum(N_strata * sigma_strata) 

effect_grid = seq(0,0.8,by=.025)
run_twosample_sims <- function(x, sample_size, n_sims = 300, stratified = FALSE){
  mu <- mean(x)
  n_x <- length(x)
  shift <- effect_grid * mu
  #scaling <- 1+effect_grid
  #spike <- cbind(matrix(0, ncol = floor(n_x * .8), nrow = length(shift)), 5*matrix(rep(shift, each = ceiling(n_x * .2)), ncol = ceiling(n_x * .2), byrow = TRUE))
  if(stratified){
    #dataframe format needed to work with sampling package
    K <- length(unique(strata))
    pop_1 <- data.frame(TC = x, strata = strata)
    #n_strata_opt <- round_strata_sizes(sample_size * strata_weights_opt)
    n_strata_prop <- round_strata_sizes(sample_size * strata_weights_prop)
    stratified_p_values <- matrix(NA, nrow = n_sims, ncol = length(shift))
  }
  normal_p_values <- matrix(NA, nrow = n_sims, ncol = length(shift))
  LMT_rejections <- matrix(NA, nrow = n_sims, ncol = length(shift))
  
  for(i in 1:n_sims){
    for(j in 1:length(shift)){
      sample_1 <- sample(x, size = sample_size, replace = TRUE)
      sample_2 <- sample(x + shift[j], size = sample_size, replace = TRUE)
      diff_mean <- mean(sample_2) - mean(sample_1)
      std_error <- sqrt(var(sample_1)/sample_size + var(sample_2)/sample_size)
      #this is for Welch's t-test
      dof <- std_error^4 / ((var(sample_1)/sample_size)^2 / (sample_size - 1) + (var(sample_2)/sample_size)^2 / (sample_size - 1))
      normal_p_values[i,j] <- pt(diff_mean/std_error, df = dof, lower.tail = FALSE)
      #assuming bounds are [0%,20%], as in mineral soils
      LMT_rejections[i,j] <- two_sample_LMT_test(sample_1 = sample_1/20, sample_2 = sample_2/20, B = 300, alpha = .05, method = "fisher")
      if(stratified){
        pop_2 <- data.frame(TC = x + shift[j], strata = strata)
        strat_sample_1 <- strata(data = pop_1, stratanames = "strata", size = n_strata_prop, method = "srswr")
        strat_sample_2 <- strata(data = pop_2, stratanames = "strata", size = n_strata_prop, method = "srswr")
        stratified_estimate_1 <- get_mean_se_stratified(sample = getdata(pop_1, strat_sample_1), N_strata = table(pop_1$strata))
        stratified_estimate_2 <- get_mean_se_stratified(sample = getdata(pop_2, strat_sample_2), N_strata = table(pop_2$strata))
        difference_estimate <- stratified_estimate_2[1] - stratified_estimate_1[1]
        combined_se <- sqrt(stratified_estimate_1[2]^2 + stratified_estimate_2[2]^2)
        stratified_p_values[i,j] <- 1 - pt(q = difference_estimate / combined_se, df = 2 * (sample_size - K))
      }
    }
  }
  normal_power_shift <- colMeans(normal_p_values < .05)
  LMT_power_shift <- colMeans(LMT_rejections)
  if(stratified){
    stratified_power_shift <- colMeans(stratified_p_values < 0.05)
    cbind("normal" = normal_power_shift, "stratified" = stratified_power_shift, "LMT"= LMT_power_shift)
  } else {
    cbind("normal" = normal_power_shift, "LMT"= LMT_power_shift)
  }
}


#since these take a while to run, I usually save them
# power_10_rangeland <- run_twosample_sims(topsoil_TC_rangeland, sample_size = 10, stratified = TRUE)
# power_30_rangeland <- run_twosample_sims(topsoil_TC_rangeland, sample_size = 30, stratified = TRUE)
# power_90_rangeland <- run_twosample_sims(topsoil_TC_rangeland, sample_size = 90, stratified = TRUE)
# power_200_rangeland <- run_twosample_sims(topsoil_TC_rangeland, sample_size = 200, stratified = TRUE)
# power_500_rangeland <- run_twosample_sims(topsoil_TC_rangeland, sample_size = 500, stratified = TRUE)
# power_10_cropland <- run_twosample_sims(topsoil_TC_cropland, sample_size = 10)
# power_30_cropland <- run_twosample_sims(topsoil_TC_cropland, sample_size = 30)
# power_90_cropland <- run_twosample_sims(topsoil_TC_cropland, sample_size = 90)
# power_200_cropland <- run_twosample_sims(topsoil_TC_cropland, sample_size = 200)
# power_500_cropland <- run_twosample_sims(topsoil_TC_cropland, sample_size = 500)
# power_frame <- bind_rows(
#   as.data.frame(power_10_rangeland) %>% mutate(sample_size = 10, relative_effect_size = effect_grid, land_use = "Rangeland"),
#   as.data.frame(power_30_rangeland) %>% mutate(sample_size = 30, relative_effect_size = effect_grid, land_use = "Rangeland"),
#   as.data.frame(power_90_rangeland) %>% mutate(sample_size = 90, relative_effect_size = effect_grid, land_use = "Rangeland"),
#   as.data.frame(power_200_rangeland) %>% mutate(sample_size = 200, relative_effect_size = effect_grid, land_use = "Rangeland"),
#   as.data.frame(power_500_rangeland) %>% mutate(sample_size = 500, relative_effect_size = effect_grid, land_use = "Rangeland"),
#   as.data.frame(power_10_cropland) %>% mutate(sample_size = 10, relative_effect_size = effect_grid, land_use = "Cropland"),
#   as.data.frame(power_30_cropland) %>% mutate(sample_size = 30, relative_effect_size = effect_grid, land_use = "Cropland"),
#   as.data.frame(power_90_cropland) %>% mutate(sample_size = 90, relative_effect_size = effect_grid, land_use = "Cropland"),
#   as.data.frame(power_200_cropland) %>% mutate(sample_size = 200, relative_effect_size = effect_grid, land_use = "Cropland"),
#   as.data.frame(power_500_cropland) %>% mutate(sample_size = 500, relative_effect_size = effect_grid, land_use = "Cropland")
#   ) %>%
#   rename("t test" = normal, "Nonparametric test" = LMT, "Stratified t test" = stratified) %>%
#   pivot_longer(cols = c("t test", "Nonparametric test", "Stratified t test"), names_to = "Test", values_to = "Power")
#save(power_frame, file = "power_frame_nonparametric_fisher_spike")
load("power_frame_nonparametric_fisher")

power_frame <- power_frame %>%
  mutate(sample_size_long = paste(sample_size, "Samples")) %>%
  mutate(sample_size_long = factor(sample_size_long, levels = paste(c(10, 30, 90, 200, 500), "Samples")))


ggplot(power_frame %>% filter(sample_size != 500), aes(x = relative_effect_size, y = Power, color = Test, linetype = Test)) +
  geom_line(size = 1.5) +
  facet_grid(sample_size_long ~ land_use) +
  theme_bw() +
  scale_color_manual(values = c("darkorange3","steelblue","forestgreen")) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  xlab("Relative effect size (additional percent TC%)") +
  coord_cartesian(xlim = c(0,.6)) +
  theme(text = element_text(size = 20))






