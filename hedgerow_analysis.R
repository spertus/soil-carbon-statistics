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

#this has only topsoil data, but includes physical indicators of soil health
#Site FBF is missing all biological data, we drop it.
hr_data <- read_csv("../Data/HR_Soil_Health_Updated.csv") %>%
  rename(site_name = `Site Name`, site = Site, hr_age = HR_age, compost = Compost, crops = Crops, cover_crop = Cover_Crop, fallow = Fallow, soil_type = Soil_Type, treatment = Treatment,  location = Location, upper_depth = Upper_Depth, profile_carbon = profileC, macroagg = Macroagg, microagg = Microagg, total_agg = Total_Agg, infiltration_dry = Infiltration_Dry, infiltration_wet = Infiltration_Wet, surface_hardness = SH_15) %>% 
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


#texture and pH should not differ between treatments but should differ between soil types
#GWC should not predict changes in soil health, microbial biomass should be controlled for GWC. Note that we do not see a significant difference in GWC between H and R
#these are names of the variables that are soil health indicators:
soil_health_vars <- c("per_C", "per_N", "profile_carbon", "POXc", "EOC", "MBC", "MBN", "glucam", "glucos", "cell", "BD", "macroagg", "microagg", "infiltration_dry", "infiltration_wet", "surface_hardness")


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
    pH = mean(pH), #this might require a modified average
    per_N = mean(per_N),
    per_C = mean(per_C),
    POXc = mean(POXc),
    GWC = mean(GWC),
    EOC = mean(EOC),
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
    surface_hardness = mean(surface_hardness, na.rm = T)
  ) %>%
  ungroup() %>%
  group_by(upper_depth) %>%
  mutate(infiltration_dry = mean_impute(infiltration_dry)) %>%
  mutate(infiltration_wet = mean_impute(infiltration_wet)) %>%
  mutate(surface_hardness = mean_impute(surface_hardness)) %>%
  ungroup() %>%
  mutate(soil_type = as_factor(soil_type)) 

topsoil_means <- soil_means %>% filter(upper_depth == 0)
subsoil_means <- soil_means %>% filter(upper_depth == 10)


#plot profile carbon
plot_profile_carbon_change <- ggplot(topsoil_means, aes(x = profile_carbon, y = site_name, color = treatment)) +
  geom_line(aes(group = site_name), color = "black") +
  geom_point(size = 3) +
  facet_grid(soil_type ~ ., scales = "free_y")

#PCAs of topsoil and subsoil data
topsoil_matrix <- topsoil_data %>%
  select(all_of(soil_health_vars), GWC, pH, sand, clay) %>%
  mutate_all(mean_impute) %>%
  as.matrix() 
subsoil_matrix <- subsoil_data %>%
  select(all_of(soil_health_vars), GWC, pH, sand, clay) %>%
  select(-macroagg, -microagg) %>%
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
plot(top_pca$x, xlim = c(-6,6),ylim=c(-6,6), col = brewer.pal(4, 'Dark2')[topsoil_data$soil_type], pch = ifelse(topsoil_data$treatment == "H", 19, 17), lwd = 1.5, cex = 1.5)
arrows(x0 = 0, y0 = 0, x1 = top_pca$rotation[,1]*10, y1 = top_pca$rotation[,2]*10, lwd = 1.5, length = .1)
text(x = top_pca$rotation[,1]*12, y = top_pca$rotation[,2]*11, labels = pca_labels, lwd = 1.5)

#subsoil PCA plot
plot(sub_pca$x, xlim = c(-6,6),ylim=c(-6,6), col = brewer.pal(4, 'Dark2')[subsoil_data$soil_type], pch = ifelse(subsoil_data$treatment == "H", 19, 17), lwd = 1.5, cex = 1.5)
arrows(x0 = 0, y0 = 0, x1 = sub_pca$rotation[,1]*10, y1 = sub_pca$rotation[,2]*10, lwd = 1.5, length = .1)
text(x = sub_pca$rotation[,1]*12, y = sub_pca$rotation[,2]*11, labels = pca_labels, lwd = 1.5)


#parametric manovas based on Gaussian assumptions
#create matrix of soil health variables from means within plots
topmeans_matrix <- topsoil_means %>%
  select(all_of(c(soil_health_vars))) %>%
  as.matrix()
submeans_matrix <- subsoil_means %>%
  select(all_of(soil_health_vars)) %>%
  select(-microagg, -macroagg) %>%
  as.matrix()

top_corr_matrix <- cor(topmeans_matrix, use = "complete.obs", method = 'spearman')
sub_corr_matrix <- cor(submeans_matrix, use = "complete.obs", method = 'spearman')
varnames_top <- c("Percent C", "Percent N", "Whole Profile SOC", "POXc", "EOC", "MBC", "MBN", "Glucaminidase", "Glucosidase", "Cellulase", "Bulk Density", "Macroaggregates", "Microaggregates", "Dry Infiltration", "Wet Infiltration", "Surface Hardness")
varnames_sub <- c("Percent C", "Percent N", "Whole Profile SOC", "POXc", "EOC", "MBC", "MBN", "Glucaminidase", "Glucosidase", "Cellulase", "Bulk Density", "Dry Infiltration", "Wet Infiltration", "Surface Hardness")
colnames(top_corr_matrix) <- rownames(top_corr_matrix) <- varnames_top
colnames(sub_corr_matrix) <- rownames(sub_corr_matrix) <- varnames_sub

corrplot(top_corr_matrix, method = "square")
corrplot(sub_corr_matrix, method = "square")

#whole profile carbon model
carbon_model <- lm(formula(paste("profile_carbon ~ ", paste(soil_health_vars[-c(3)], collapse = "+"))), data = topsoil_means)

#2-way
#we might want to take out soil texture
twoway_manova_model <- lm(topmeans_matrix ~ treatment * soil_type, data = topsoil_means)
summary(manova(twoway_manova_model))



#non-parametric paired one-way manova
topsoil_means_differences <- topsoil_means %>%
  select(site_name, soil_type, treatment, all_of(soil_health_vars)) %>%
  pivot_wider(names_from = treatment, values_from = all_of(soil_health_vars)) %>%
  mutate(
    diff_per_N = per_N_H - per_N_R,
    diff_per_C = per_C_H - per_C_R,
    diff_BD = BD_H - BD_R,
    diff_POXc = POXc_H - POXc_R,
    diff_EOC = EOC_H - EOC_R,
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
  select(-starts_with(soil_health_vars))
soil_type <- as.numeric(topsoil_means_differences$soil_type)

diff_matrix <- topsoil_means_differences %>%
  select(-site_name, -soil_type) %>%
  as.matrix()
diff_means <- apply(diff_matrix, 2, mean)
#ANOVA test statistic 
#for an example of this see the permuter github page https://github.com/statlab/permuter/blob/master/vignettes/examples_chapters1_4.Rmd
get_ANOVA_soiltype <- function(dependent_variable){
  block <- soil_type
  group_means <- tapply(dependent_variable, block, mean)
  group_sizes <- as.numeric(table(block))
  sum(group_sizes * group_means^2)
}
original_ANOVAs <- apply(diff_matrix, 2, get_ANOVA_soiltype)

#test for main effect of treatment and interaction with soil health 
paired_perm_dist <- apply(diff_matrix, 2, one_sample, reps = 10000)
interaction_perm_dist <- apply(diff_matrix, 2, function(x){k_sample(x = x, group = as.numeric(topsoil_means_differences$soil_type), reps = 10000)})
paired_p_values <- rep(1, length(diff_means))
interaction_p_values <- rep(1, length(diff_means))
for(i in 1:length(paired_p_values)){
  paired_p_values[i] <- t2p(diff_means[i], paired_perm_dist[,i], alternative = "two-sided")
  interaction_p_values[i] <- t2p(original_ANOVAs[i], interaction_perm_dist[,i], alternative = "two-sided")
}
#p-value for the intersection null, that there is no effect on any SH variable
#nonparametric one-way MANOVA
combined_paired_p_value <- npc(diff_means, distr = paired_perm_dist, combine = "fisher", alternatives = "two-sided")
combined_interaction_p_value <- npc(original_ANOVAs, distr = interaction_perm_dist, combine = "fisher", alternatives = "two-sided")
#p-values for partial tests, adjusted for multiplicity by closed testing.
closed_paired_p_values <- fwe_minp(paired_p_values, paired_perm_dist)


p_value_frame <- data.frame("variable" = colnames(diff_matrix), "difference_in_means" = diff_means, "p_value" = paired_p_values, "adjusted_p_value" = closed_paired_p_values)
rownames(p_value_frame) <- NULL


#TODO: nonparametric two-way MANOVA with soil type
#THIS NEEDS WORK!
shuffle <- function(x){sample(x, size = length(x), replace = TRUE)}
permutation_twoway_ANOVA <- function(y, treatment, block, B = 1000){
  n_treatment <- as.numeric(table(treatment))
  n_block <- as.numeric(table(block))
  grand_mean <- mean(y)
  original_treatment_means <- tapply(y, treatment, mean)
  original_block_means <- tapply(y, block, mean)
  original_test_stat <- sum(n_treatment * (original_treatment_means - grand_mean)^2) / sum(y - original_treatment_means - original_block_means + grand_mean)
  perm_dist <- rep(NA, B)
  for(i in 1:B){
    block_shuffled_y <- tapply(y, block, shuffle) %>% reduce(c)
    shuffled_treatment_mean <- tapply(block_shuffled_y, treatment, mean)
    perm_dist[i] <- sum(n_treatment * (shuffled_treatment_mean - grand_mean)^2)
  }
}


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


