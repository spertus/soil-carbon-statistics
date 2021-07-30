#Hedgerow analysis
#Data from Jessica Chiartas
#Analysis by Jacob Spertus
library(tidyverse)
library(corrplot)
library(permuter)
library(car)
library(M3C)
library(RColorBrewer)
set.seed(1337)
source("functions.R")


############# Load data ##############
#this has only topsoil data, but includes physical indicators of soil health
#Site FBF is missing all biological data, we drop it.
hr_data <- read_csv("../Data/HR_Soil_Health_Updated.csv") %>%
  rename(site_name = `Site Name`, site = Site, hr_age = HR_age, compost = Compost, crops = Crops, cover_crop = Cover_Crop, fallow = Fallow, soil_type = Soil_Type, treatment = Treatment,  location = Location, upper_depth = Upper_Depth, profile_carbon = profileC, macroagg = Macroagg, microagg = Microagg, total_agg = Total_Agg, infiltration_dry = Infiltration_Dry, infiltration_wet = Infiltration_Wet, surface_hardness = SH_15, subsurface_hardness = SSH_30) %>% 
  filter(site_name != "FBF") %>%
  filter(!is.na(site_name)) #there's a tail of all NA values, which I'm deleting

carbon_data <- read_csv("../Data/HR_soil_carbon.csv") %>%
  rename(site_name = `Site Name`, site = Site, hr_age = HR_age, compost = Compost, crops = Crops, cover_crop = Cover_Crop, fallow = Fallow, soil_type = Soil_Type, treatment = Treatment,  location = Location, upper_depth = Upper_Depth, lower_depth = Lower_Depth, length = Length) %>% 
  filter(site_name != "FBF") 


topsoil_data <- hr_data %>% filter(upper_depth == 0)
subsoil_data <- hr_data %>% filter(upper_depth == 10)


#function to fill in the missing values in a vector with the mean of observed values
mean_impute <- function(column){
  column_mean <- mean(column, na.rm = TRUE)
  column[is.na(column)] <- column_mean
  column
}

################# Data handling ##############
#texture and pH should not differ between treatments but should differ between soil types
#GWC should not predict changes in soil health, microbial biomass should be controlled for GWC. Note that we do not see a significant difference in GWC between H and R
#these are names of the variables that are soil health indicators:
soil_health_vars <- c("per_C", "per_N", "profile_carbon", "POXc", "EOC", "EON", "MBC", "MBN", "glucam", "glucos", "cell", "BD")


#the site FBF is missing all biological information, so delete it for now.
#CLBL and Muller Co do not have infiltration data or hardness
soil_means <- hr_data %>%
  group_by(site_name, treatment, upper_depth) %>%
  summarize(
    HR_age = first(hr_age), 
    compost = first(compost), 
    crops = first(crops), 
    cover_crop = first(cover_crop),
    fallow = first(fallow),
    soil_type = first(soil_type),
    sand = mean(sand),
    clay = mean(clay),
    silt = mean(silt),
    silt_clay = mean(silt_clay),
    BD = mean(BD),
    pH = -log10(mean(10^(-pH))),
    per_N = mean(per_N),
    per_C = mean(per_C),
    POXc = mean(POXc),
    GWC = mean(GWC),
    EOC = mean(EOC),
    EON = mean(EON),
    MBC = mean(MBC), 
    MBN = mean(MBN),
    profile_carbon = mean(profile_carbon),
    glucam = mean(glucam),
    glucos = mean(glucos),
    cell = mean(cell),
    macroagg = mean(macroagg),
    microagg = mean(microagg),
    totalagg = mean(total_agg),
    infiltration_dry = mean(infiltration_dry, na.rm = T), #infiltration and surface hardness are missing for a few sampling locations
    infiltration_wet = mean(infiltration_wet, na.rm = T),
    surface_hardness = mean(surface_hardness, na.rm = T),
    subsurface_hardness = mean(subsurface_hardness, na.rm = T)
  ) %>%
  ungroup() %>%
  group_by(upper_depth) %>%
  mutate(infiltration_dry = mean_impute(infiltration_dry)) %>%
  mutate(infiltration_wet = mean_impute(infiltration_wet)) %>%
  mutate(surface_hardness = mean_impute(surface_hardness)) %>%
  mutate(subsurface_hardness = mean_impute(subsurface_hardness)) %>%
  ungroup() %>%
  mutate(soil_type = as_factor(soil_type)) 

carbon_means <- carbon_data %>% 
  group_by(site_name, soil_type, treatment, upper_depth) %>%
  summarize(C_stock = mean(C_Mg_ha, na.rm = TRUE), per_C = mean(per_C, na.rm = TRUE)) %>% 
  pivot_wider(names_from = upper_depth, values_from = c(C_stock, per_C)) %>%
  ungroup() %>%
  mutate(C_stock_75 = mean_impute(C_stock_75), per_C_75 = mean_impute(per_C_75))


topsoil_means <- soil_means %>% filter(upper_depth == 0)
subsoil_means <- soil_means %>% filter(upper_depth == 10)


#summarize profile carbon
carbon_means_long <- carbon_data %>%
  group_by(site_name, soil, treatment, upper_depth) %>%
  summarize(C_stock = mean(C_Mg_ha, na.rm = TRUE), per_C = mean(per_C, na.rm = TRUE)) %>%
  ungroup()

plot_profile_carbon_change <- ggplot(carbon_means_long, aes(x = C_stock, y = site_name, color = treatment)) +
  geom_line(aes(group = site_name), color = "black") +
  geom_point(size = 3) +
  facet_grid(upper_depth ~ . , scales = "free")

#estimate soil carbon change
carbon_table <- carbon_means_long %>%
  pivot_wider(names_from = treatment, values_from = c("C_stock", "per_C")) %>%
  group_by(upper_depth) %>%
  summarize(mean_diff_stock = mean(C_stock_H - C_stock_R, na.rm = T), mean_diff_per_C = mean(per_C_H - per_C_R, na.rm = T))

#plot carbon change
carbon_pct_change_boxplot <- ggplot(data = carbon_means_long %>% mutate(Treatment = ifelse(treatment == "H", "Hedgerow", "Row")), aes(x = as_factor(upper_depth), y = per_C, fill = Treatment)) +
  geom_boxplot() +
  facet_grid(~ soil) +
  xlab("Depth (cm)") +
  ylab("Percent SOC") +
  theme_bw() +
  theme(text = element_text(size = 16))


########### analysis of variance ###########
#total carbon ANOVAs for table
total_var <- carbon_data %>%
  group_by(upper_depth) %>%
  summarize(total_var = var(C_Mg_ha, na.rm = TRUE)) 
treatment_var <- carbon_data %>% 
  group_by(upper_depth, treatment) %>%
  summarize(within_treatment_variance = var(C_Mg_ha, na.rm = TRUE)) %>%
  group_by(upper_depth) %>%
  summarize(mean(within_treatment_variance, na.rm = TRUE))
site_var <- carbon_data %>% 
  group_by(upper_depth, site) %>%
  summarize(within_site_variance = var(C_Mg_ha, na.rm = TRUE)) %>%
  group_by(upper_depth) %>%
  summarize(mean(within_site_variance, na.rm = TRUE))
soil_type_var <- carbon_data %>%
  group_by(upper_depth, soil_type) %>%
  summarize(within_soil_type_variance = var(C_Mg_ha, na.rm = TRUE)) %>%
  group_by(upper_depth) %>%
  summarize(mean(within_soil_type_variance, na.rm = TRUE))



########## Carbon analysis #########
#I've spot checked that this aligns with soil_means$profile_carbon, as it should
wp_carbon_stocks <- carbon_means_long %>%
  group_by(site_name, treatment, soil) %>%
  summarize(wp_stock = sum(C_stock)) %>%
  ungroup()
  
#mean difference in stocks
mean_stock_difference <- wp_carbon_stocks %>% 
  pivot_wider(names_from = treatment, values_from = wp_stock, names_prefix = "wp_stock_") %>%
  summarize(mean_diff = mean(wp_stock_H - wp_stock_R, na.rm = T), sd_diff = sd(wp_stock_H - wp_stock_R, na.rm = T), std_error_diff = sd(wp_stock_H - wp_stock_R, na.rm = T) / sqrt(n())) 

carbon_stock_change_boxplot <- ggplot(data = wp_carbon_stocks %>% mutate(Treatment = ifelse(treatment == "H", "Hedgerow", "Row")), aes(x = soil, y = wp_stock, fill = Treatment)) +
  geom_boxplot() +
  xlab("Depth (cm)") +
  ylab("SOC Stock (Mg / ha)") +
  theme_bw() +
  theme(text = element_text(size = 16))
  

############ PCA and correlation analysis ##########
#PCAs of topsoil and subsoil data
topsoil_matrix <- topsoil_data %>%
  select(all_of(soil_health_vars), macroagg, microagg, infiltration_dry, infiltration_wet, surface_hardness, GWC, pH, sand, clay) %>%
  mutate_all(mean_impute) %>%
  as.matrix() 
subsoil_matrix <- subsoil_data %>%
  select(all_of(soil_health_vars), subsurface_hardness, GWC, pH, sand, clay) %>%
  mutate_all(mean_impute) %>%
  as.matrix() 

organic_sites <- c("Fong E", "Fong W", "CLBL", "FBF", "Gilmer", "JRR")
pca_labels <- gsub("_", " ", colnames(topsoil_matrix))
pca_labels[pca_labels == "MBC"] <- ""
top_pca <- prcomp(topsoil_matrix, scale = TRUE)
sub_pca <- prcomp(subsoil_matrix, scale = TRUE)
#variance explained by topsoil PCs
prop_explained <- top_pca$sdev^2 / sum(top_pca$sdev^2)
plot(prop_explained, type = 'b')

#topsoil PCA plot
plot(top_pca$x, xlim = c(-11,6),ylim=c(-6,6), col = brewer.pal(4, 'Dark2')[topsoil_data$soil_type], pch = ifelse(topsoil_data$treatment == "H", 19, 17), lwd = 1.5, cex = 1.5)
arrows(x0 = 0, y0 = 0, x1 = top_pca$rotation[,1]*10, y1 = top_pca$rotation[,2]*10, lwd = 1.5, length = .1)
text(x = top_pca$rotation[,1]*12, y = top_pca$rotation[,2]*11, labels = pca_labels, lwd = 1.5)
legend(x = 2, y = 6, col = brewer.pal(4, 'Dark2'), pch = 19, legend = c("Yolo", "Brentwood", "Capay", "Corning"))
legend(x = -10, y = 6, pch = c(19, 17), legend = c("Hedgerow","Crop"))

#subsoil PCA plot
plot(sub_pca$x, xlim = c(-6,6),ylim=c(-6,6), col = brewer.pal(4, 'Dark2')[subsoil_data$soil_type], pch = ifelse(subsoil_data$treatment == "H", 19, 17), lwd = 1.5, cex = 1.5)
arrows(x0 = 0, y0 = 0, x1 = sub_pca$rotation[,1]*10, y1 = sub_pca$rotation[,2]*10, lwd = 1.5, length = .1)
text(x = sub_pca$rotation[,1]*12, y = sub_pca$rotation[,2]*11, labels = pca_labels, lwd = 1.5)



#create matrix of soil health variables from means within plots
topmeans_matrix <- topsoil_means %>%
  select(all_of(c(soil_health_vars)), macroagg, microagg, infiltration_dry, infiltration_wet, surface_hardness, -profile_carbon) %>%
  as.matrix()
submeans_matrix <- subsoil_means %>%
  select(all_of(soil_health_vars), subsurface_hardness, -profile_carbon) %>%
  as.matrix()
carbon_matrix <- carbon_means %>%
  select(C_stock_0, C_stock_10, C_stock_20, C_stock_50, C_stock_75) %>%
  as.matrix()

top_corr_matrix <- cor(topmeans_matrix, use = "complete.obs", method = 'spearman')
sub_corr_matrix <- cor(submeans_matrix, use = "complete.obs", method = 'spearman')
total_corr_matrix <- cor(cbind(topmeans_matrix, submeans_matrix, carbon_matrix), use = "complete.obs", method = "spearman")
varnames_top <- c("Percent C", "Percent N", "POXc", "EOC", "MBC", "MBN", "Glucaminidase", "Glucosidase", "Cellulase", "Bulk Density", "Macroaggregates", "Microaggregates", "Dry Infiltration", "Wet Infiltration", "Surface Hardness")
varnames_sub <- c("Percent C", "Percent N", "POXc", "EOC", "MBC", "MBN", "Glucaminidase", "Glucosidase", "Cellulase", "Bulk Density", " Surface Hardness")
varnames_carbon <- c("C Stock 0-10", "C Stock 10-20", "C Stock 20-50", "C Stock 50-75", "C Stock 75-100")
colnames(top_corr_matrix) <- rownames(top_corr_matrix) <- varnames_top
colnames(sub_corr_matrix) <- rownames(sub_corr_matrix) <- varnames_sub
colnames(total_corr_matrix) <- rownames(total_corr_matrix) <- c(paste("Top", varnames_top), paste("Sub",varnames_sub), varnames_carbon)

corrplot(top_corr_matrix, method = "square")
corrplot(sub_corr_matrix, method = "square")
corrplot(total_corr_matrix, method = "square")


############ Pairwise comparison of soil types to each other on a few variables ########
#pairwise comparison of carbon across soil types stratified by row/hedgerow and combined
whole_profile_stock <- rowSums(carbon_matrix)

#helper function to run a pairwise comparison of the hypothesis that variable x does not differ between soil types for both rows and hedgerows
run_pairwise_comparison <- function(x){
  comparison_matrix <- matrix(NA, nrow = 4, ncol = 4)
  colnames(comparison_matrix) <- c("Yolo", "Brentwood", "Capay", "Corning")
  land_use <- carbon_means$treatment
  soil_types <- carbon_means$soil_type
  for(i in 1:4){
    if(i > 1){
      for(j in 1:(i-1)){
        row_data_0 <- x[land_use == "R"][soil_types[land_use == "R"] == i]
        row_data_1 <- x[land_use == "R"][soil_types[land_use == "R"] == j]
        diff_mean_rows <- mean(row_data_1) - mean(row_data_0) 
        hedgerow_data_0 <- x[land_use == "H"][soil_types[land_use == "H"] == i]
        hedgerow_data_1 <- x[land_use == "H"][soil_types[land_use == "H"] == j]
        diff_mean_hedgerows <- mean(hedgerow_data_1) - mean(hedgerow_data_0) 
        perm_dist <- lockstep_two_sample(x_matrix = cbind(row_data_0, hedgerow_data_0), y_matrix = cbind(row_data_1, hedgerow_data_1), exact = TRUE)
        combined_p_value <- npc(statistics = c(diff_mean_rows, diff_mean_hedgerows), distr = perm_dist, combine = "fisher", alternatives = "two-sided")
        comparison_matrix[i,j] <- combined_p_value
      }
    }
  }
  comparison_matrix
}

# pairwise_carbon_wp <- run_pairwise_comparison(whole_profile_stock)
# pairwise_MBC_top <- run_pairwise_comparison(topmeans_matrix[,'MBC'])
# pairwise_MBN_top <- run_pairwise_comparison(topmeans_matrix[,'MBN'])
# pairwise_POXc_top <- run_pairwise_comparison(topmeans_matrix[,'POXc'])
# pairwise_cellulase_top <-  run_pairwise_comparison(topmeans_matrix[,'cell'])
# pairwise_macroagg_top <-  run_pairwise_comparison(topmeans_matrix[,'macroagg'])
# pairwise_microagg_top <-  run_pairwise_comparison(topmeans_matrix[,'microagg'])


############# MANOVAs for management, soil type, and interaction ###########
#non-parametric paired one-way manova
#topsoil differences
topsoil_means_differences <- topsoil_means %>%
  select(site_name, soil_type, treatment, all_of(soil_health_vars), macroagg, microagg, infiltration_dry, infiltration_wet, surface_hardness) %>%
  pivot_wider(names_from = treatment, values_from = all_of(c(soil_health_vars, "macroagg", "microagg", "infiltration_dry", "infiltration_wet", "surface_hardness"))) %>%
  mutate(
    diff_per_N = per_N_H - per_N_R,
    diff_per_C = per_C_H - per_C_R,
    diff_BD = BD_H - BD_R,
    diff_POXc = POXc_H - POXc_R,
    diff_EOC = EOC_H - EOC_R,
    diff_EON = EON_H - EON_R,
    diff_MBC = MBC_H - MBC_R,
    diff_MBN = MBN_H - MBN_R,
    diff_glucam = glucam_H - glucam_R,
    diff_glucos = glucos_H - glucos_R,
    diff_cell = cell_H - cell_R,
    diff_macroagg = macroagg_H - macroagg_R,
    diff_microagg = microagg_H - microagg_R,
    diff_infiltration_dry = infiltration_dry_H - infiltration_dry_R,
    diff_infiltration_wet = infiltration_wet_H - infiltration_wet_R,
    diff_surface_hardness = surface_hardness_H - surface_hardness_R) %>%
  arrange(site_name) %>%
  select(site_name, soil_type, starts_with("diff_"))
#subsoil differences
subsoil_means_differences <- subsoil_means %>%
  select(site_name, soil_type, treatment, all_of(soil_health_vars), subsurface_hardness) %>%
  pivot_wider(names_from = treatment, values_from = all_of(c(soil_health_vars, "subsurface_hardness"))) %>%
  mutate(
    diff_per_N = per_N_H - per_N_R,
    diff_per_C = per_C_H - per_C_R,
    diff_BD = BD_H - BD_R,
    diff_POXc = POXc_H - POXc_R,
    diff_EOC = EOC_H - EOC_R,
    diff_EON = EON_H - EON_R,
    diff_MBC = MBC_H - MBC_R,
    diff_MBN = MBN_H - MBN_R,
    diff_glucam = glucam_H - glucam_R,
    diff_glucos = glucos_H - glucos_R,
    diff_cell = cell_H - cell_R,
    diff_subsurface_hardness = subsurface_hardness_H - subsurface_hardness_R) %>%
  arrange(site_name) %>%
  select(starts_with("diff_"))
#carbon differences
carbon_means_differences <- carbon_means %>%
  select(-starts_with("per_C_")) %>%
  pivot_wider(names_from = treatment, values_from = starts_with("C_stock_")) %>%
  mutate(diff_stock_0 = C_stock_0_H - C_stock_0_R) %>%
  mutate(diff_stock_10 = C_stock_10_H - C_stock_10_R) %>%
  mutate(diff_stock_20 = C_stock_20_H - C_stock_20_R) %>%
  mutate(diff_stock_50 = C_stock_50_H - C_stock_50_R) %>%
  mutate(diff_stock_75 = C_stock_75_H - C_stock_75_R) %>%
  arrange(site_name) %>%
  select(starts_with("diff_"))

soil_type <- as.numeric(topsoil_means_differences$soil_type)

diff_matrix_top <- topsoil_means_differences %>%
  select(-site_name, -soil_type) %>%
  as.matrix()
colnames(diff_matrix_top) <- paste("top_", gsub("diff_", "", colnames(diff_matrix_top)), sep = "")
diff_matrix_sub <- subsoil_means_differences %>%
  as.matrix()
colnames(diff_matrix_sub) <- paste("sub_", gsub("diff_", "", colnames(diff_matrix_sub)), sep = "")
diff_matrix_carbon <- carbon_means_differences %>%
  as.matrix()



#carbon stock is not included in this analysis
#differences are used for management and interaction analysis
diff_matrix <- cbind(diff_matrix_top, diff_matrix_sub)
#original means are used for soil types
mean_matrix <- cbind(topmeans_matrix, submeans_matrix)
colnames(mean_matrix) <- c(paste("top", colnames(topmeans_matrix), sep = "_"), paste("sub", colnames(submeans_matrix), sep = "_"))

#compute original test statistics
original_diff_means <- apply(diff_matrix, 2, mean)
original_soil_type_ANOVAs <- apply(mean_matrix, 2, get_ANOVA, group = topsoil_means$soil_type)
original_interaction_ANOVAs <- apply(diff_matrix, 2, get_ANOVA, group = soil_type)


#simulate from permutation distributions
paired_perm_dist <- lockstep_one_sample(diff_matrix, reps = 1e5)
soil_type_perm_dist <- lockstep_ANOVA(mean_matrix, group = topsoil_means$soil_type, reps = 10000)
interaction_perm_dist <- lockstep_ANOVA(diff_matrix, group = soil_type, reps = 1000)

#compute p-values from original test statistics and simulated permutation distributions
paired_p_values <- rep(1, length(original_diff_means))
soil_type_p_values <- rep(1, length(original_diff_means))
interaction_p_values <- rep(1, length(original_diff_means))
for(i in 1:length(paired_p_values)){
  paired_p_values[i] <- t2p(original_diff_means[i], paired_perm_dist[,i], alternative = "two-sided")
  soil_type_p_values[i] <- t2p(original_soil_type_ANOVAs[i], soil_type_perm_dist[,i], alternative = "greater")
  interaction_p_values[i] <- t2p(original_interaction_ANOVAs[i], interaction_perm_dist[,i], alternative = "greater")
}

#p-value for the intersection null, that there is no effect on any SH variable
#nonparametric one-way MANOVA
combined_paired_p_value <- npc(original_diff_means, distr = paired_perm_dist, combine = "fisher", alternative = "two-sided")
combined_soil_type_p_value <- npc(original_soil_type_ANOVAs, distr = soil_type_perm_dist, combine = "fisher", alternative = "two-sided")
combined_interaction_p_value <- npc(original_interaction_ANOVAs, distr = interaction_perm_dist, combine = "fisher", alternatives = "two-sided")

#p-values for partial tests, adjusted for multiplicity by closed testing.
closed_paired_p_values <- fwe_minp(paired_p_values, paired_perm_dist)


p_value_frame <- data.frame("variable" = colnames(diff_matrix), "difference_in_means" = original_diff_means, "p_value" = paired_p_values, "adjusted_p_value" = closed_paired_p_values)
rownames(p_value_frame) <- NULL


#Negative controls
#check that there aren't differences in texture, pH, and GWC
control_matrix <- topsoil_means %>%
  select(site_name, treatment, sand, clay, GWC, pH) %>%
  pivot_wider(names_from = treatment, values_from = all_of(c("sand","clay","GWC","pH"))) %>%
  mutate(
    diff_sand = sand_H - sand_R,
    diff_clay = clay_H - clay_R,
    diff_GWC = GWC_H - GWC_R,
    diff_pH = pH_H - pH_R
  ) %>%
  select(starts_with("diff")) %>%
  as.matrix()

control_diff_means <- apply(control_matrix, 2, mean)
control_perm_dist <- apply(control_matrix, 2, one_sample, reps = 10000)
control_p_values <- rep(1, length(control_diff_means))
for(i in 1:length(control_p_values)){
  control_p_values[i] <- t2p(control_diff_means[i], control_perm_dist[,i], alternative = "two-sided")
}
control_p_value_frame <- data.frame(variable = colnames(control_matrix), diff_mean = control_diff_means, p_value = control_p_values)
rownames(control_p_value_frame) <- NULL
