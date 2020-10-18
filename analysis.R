source("functions.R")
library(tidyverse)
library(xtable)
################ reading data in ###############
#rs9 data is from LBL
rs9_data <- read_csv("../Data/RS92961_RS00681_JulyOctober2019surveys_drive.csv") %>%
  slice(-1) %>% #1st row is details/notes, not values 
  slice(-(310:n())) #there are extra rows from the conversion to csv
#Marin data is from silver lab
marin_data <- read_csv("../Data/mcp_duplicates.csv") %>%
  filter(!(substr(sample_name, 1, 4) %in% c('atro','blan','Blan','dumm','Dumm'))) #remove blanks, dummies, standards

############ estimating Marin plots parameters ############
#Marin plot parameters
marin_summary <- marin_data %>%
  mutate(depth = str_extract(sample_name, "0-10|10-30|30-50|50-100")) %>%
  mutate(plot = str_split(sample_name, "-", simplify = T)[,1]) %>%
  mutate(plot = str_replace(plot, " ", "")) %>%
  mutate(plot = str_to_upper(plot)) %>%
  group_by(sample_name) %>%
  summarize(carbon = mean(carbon), plot = first(plot), depth = first(depth)) %>% #mean over replicate measurements
  filter(plot != "TG2") %>% #plot TG2 has too few samples to estimate sigma_p
  group_by(plot, depth) %>%
  summarize(mu = mean(carbon), sigma_p = sd(carbon)) %>% #mean over cores
  ungroup() %>% 
  group_by(depth) %>%
  summarize(median_mu = median(mu), median_sigma_p = median(sigma_p))
  

################### estimating DC-EA error ################
#DC-EA measurement error based on Marin Data
marin_DCEA_error <- marin_data %>% 
  group_by(sample_name) %>%
  filter(!any(carbon == 0)) %>% #drop samples where any carbon measurements were 0
  filter(n() > 1) %>% #drop unduplicated samples
  summarize(variance = var(carbon), mean = mean(carbon), R_i = n()) %>%
  mutate(measurement_variance_estimate = variance / (mean^2 - variance / R_i)) %>%
  summarize(avg_measurement_variance = mean(measurement_variance_estimate, na.rm = T), median_measurement_variance = median(measurement_variance_estimate, na.rm = T))
#average sigma_delta_DCEA estimate is 0.144, median is .022
sigma_delta_dcea <- sqrt(marin_DCEA_error$avg_measurement_variance)


####################### estimating LOI error #######################
measurement_data <- rs9_data %>%
  select(LOI = 27, DC = 29, H2O = 26, sand = 30, silt = 31, clay = 32, pH = 33) %>%
  mutate_all(as.numeric)

LOI_DC_plot <- ggplot(data = measurement_data, aes(x = LOI, y = DC)) +
  geom_point() +
  geom_smooth(method = "lm")

#drop high carbon outlier
measurement_data <- measurement_data %>% filter(DC < 15)

#predicting DC just from LOI
LOI_DC_model <- lm(DC ~ LOI, data = measurement_data)
#parameter estimates: intercept -0.33, slope 0.62, R-squared = 0.964, RSE = 0.32


#what is a usual SOC concentration in Marin?
marin_SOC_bounds <- quantile(marin_data$carbon, c(.025,.975))

#95% of SOC measurments in marin are between 0.14% and 6%
#get index of DC measurements from LBL that are in range of marin data
LBL_index <- which(measurement_data$DC > marin_SOC_bounds[1] & measurement_data$DC < marin_SOC_bounds[2])
LBL_residuals <- residuals(LOI_DC_model)[LBL_index]
LBL_fitted <- fitted(LOI_DC_model)[LBL_index]

#note that everything below is in percents
hist(LBL_residuals, breaks = 30)
RMSE_loi <- sd(LBL_residuals) #.31

sigma_delta_loi <- RMSE_loi * sigma_delta_dcea / marin_summary$median_mu + sigma_delta_dcea + RMSE_loi / marin_summary$median_mu

################## estimating MIRS error #################

#from England and Viscarra-Rossel 2018
RMSE_mirs <- 0.11
sigma_delta_mirs <- RMSE_mirs * sigma_delta_dcea / marin_summary$median_mu + sigma_delta_dcea + RMSE_mirs / marin_summary$median_mu



################ generating figures ############
#average concentration and heterogeneity in topsoil (0-10 cm)
mu_top <- marin_summary$median_mu[1]
sigma_p_top <- marin_summary$median_sigma_p[1]
#average concentration and heterogeneity in deep soil (50-100 cm)
mu_deep <- marin_summary$median_mu[4]
sigma_p_deep <- marin_summary$median_sigma_p[4]
#measurement errors
sigma_delta_dcea <- 0.14
#constant additive error leads to better multiplicative error in topsoil then deep soil because low avg measurements are pertrubed more on a multiplicative scale 
sigma_delta_loi_top <- sigma_delta_loi[1]
sigma_delta_mirs_top <- sigma_delta_mirs[1]
sigma_delta_loi_deep <- sigma_delta_loi[4]
sigma_delta_mirs_deep <- sigma_delta_mirs[4]
#costs
C_0 <- 200
cost_c_low <- 5.00
cost_c_medium <- 20.00
cost_c_high <- 40.00
cost_P_loi <- 8.00
cost_P_dcea <- 11.00
cost_P_mirs <- 9.00
cost_M_dcea <- 15.00
cost_M_loi <- 1.25
cost_M_mirs <- 1.30


B_grid <- expand.grid(B = 270:5000, cost_c = c(cost_c_low, cost_c_medium, cost_c_high))


#from Vos 2005: LOI is estimated to be 1/12 the price of DC-EA.
#Vos estimates .20$ per LOI measurement and 2.35 for DC-EA
#root-mean squared predictive error (basically the measurement variance) is ~0.68% for LOI (when predicting DC-EA)


############# optimal composite sizes (topsoil) ############
optimal_composite_dcea_top <- get_optimal_composite_size(sigma_p = sigma_p_top, mu = mu_top, sigma_delta = sigma_delta_dcea, cost_c = c(cost_c_low, cost_c_medium, cost_c_high), cost_P = cost_P_dcea, cost_M = cost_M_dcea)
optimal_composite_loi_top <- get_optimal_composite_size(sigma_p = sigma_p_top, mu = mu_top, sigma_delta = sigma_delta_loi_top, cost_c = c(cost_c_low, cost_c_medium, cost_c_high), cost_P = cost_P_loi, cost_M = cost_M_loi)
optimal_composite_mirs_top <- get_optimal_composite_size(sigma_p = sigma_p_top, mu = mu_top, sigma_delta = sigma_delta_mirs_top, cost_c = c(cost_c_low, cost_c_medium, cost_c_high), cost_P = cost_P_mirs, cost_M = cost_M_mirs)

optimal_composite_dcea_deep <- get_optimal_composite_size(sigma_p = sigma_p_deep, mu = mu_deep, sigma_delta = sigma_delta_dcea, cost_c = c(cost_c_low, cost_c_medium, cost_c_high), cost_P = cost_P_dcea, cost_M = cost_M_dcea)
optimal_composite_loi_deep <- get_optimal_composite_size(sigma_p = sigma_p_deep, mu = mu_deep, sigma_delta = sigma_delta_loi_deep, cost_c = c(cost_c_low, cost_c_medium, cost_c_high), cost_P = cost_P_loi, cost_M = cost_M_loi)
optimal_composite_mirs_deep <- get_optimal_composite_size(sigma_p = sigma_p_deep, mu = mu_deep, sigma_delta = sigma_delta_mirs_deep, cost_c = c(cost_c_low, cost_c_medium, cost_c_high), cost_P = cost_P_mirs, cost_M = cost_M_mirs)



################# optimal variance under range of budgets (topsoil) ##########
optima_variance_dcea_top <- get_minimum_error(sigma_p = sigma_p_top, sigma_delta = sigma_delta_dcea, mu = mu_top, C_0 = C_0, cost_c = B_grid$cost_c, cost_P = cost_P_dcea, cost_M = cost_M_dcea, B = B_grid$B) %>%
  bind_cols(B_grid) %>%
  mutate(M = "DC-EA at 26.00 USD") %>%
  mutate(depth = "Topsoil (0-10 cm)")

optima_variance_loi_top <- get_minimum_error(sigma_p = sigma_p_top, sigma_delta = sigma_delta_loi_top, mu = mu_top, C_0 = C_0, cost_c = B_grid$cost_c, cost_P = cost_P_loi, cost_M = cost_M_loi, B = B_grid$B) %>%
  bind_cols(B_grid) %>%
  mutate(M = "LOI at 9.25 USD") %>%
  mutate(depth = "Topsoil (0-10 cm)")

optima_variance_mirs_top <- get_minimum_error(sigma_p = sigma_p_top, sigma_delta = sigma_delta_mirs_top, mu = mu_top,C_0 = C_0, cost_c = B_grid$cost_c, cost_P = cost_P_mirs, cost_M = cost_M_mirs, B = B_grid$B) %>%
  bind_cols(B_grid) %>%
  mutate(M = "MIRS at 10.30 USD") %>%
  mutate(depth = "Topsoil (0-10 cm)")

optima_variance_no_error_top <- get_minimum_error(sigma_p = sigma_p_top, sigma_delta = 0, mu = mu_top, C_0 = C_0, cost_c = B_grid$cost_c, cost_P = 0, cost_M = 0, B = B_grid$B, measurement_error = FALSE) %>%
  bind_cols(B_grid) %>%
  mutate(M = "No Measurement Error") %>%
  mutate(depth = "Topsoil (0-10 cm)")

optima_variance_top <- bind_rows(optima_variance_dcea_top, optima_variance_loi_top, optima_variance_mirs_top, optima_variance_no_error_top ) %>%
  mutate(cv = sqrt(optimum_variance) / mu_top) #add coefficient of variation


############ optimal variance under range of budgets (deep soi) ##########
optima_variance_dcea_deep <- get_minimum_error(sigma_p = sigma_p_deep, sigma_delta = sigma_delta_dcea, mu = mu_deep, C_0 = C_0, cost_c = B_grid$cost_c, cost_P = cost_P_dcea, cost_M = cost_M_dcea, B = B_grid$B) %>%
  bind_cols(B_grid) %>%
  mutate(M = "DC-EA at 26.00 USD") %>%
  mutate(depth = "Deep Soil (50-100 cm)")

optima_variance_loi_deep <- get_minimum_error(sigma_p = sigma_p_deep, sigma_delta = sigma_delta_loi_deep, mu = mu_deep, C_0 = C_0, cost_c = B_grid$cost_c, cost_P = cost_P_loi, cost_M = cost_M_loi, B = B_grid$B) %>%
  bind_cols(B_grid) %>%
  mutate(M = "LOI at 9.25 USD") %>%
  mutate(depth = "Deep Soil (50-100 cm)")

optima_variance_mirs_deep <- get_minimum_error(sigma_p = sigma_p_deep, sigma_delta = sigma_delta_mirs_deep, mu = mu_deep, C_0 = C_0, cost_c = B_grid$cost_c, cost_P = cost_P_mirs, cost_M = cost_M_mirs, B = B_grid$B) %>%
  bind_cols(B_grid) %>%
  mutate(M = "MIRS at 10.30 USD") %>%
  mutate(depth = "Deep Soil (50-100 cm)")

optima_variance_no_error_deep <- get_minimum_error(sigma_p = sigma_p_deep, sigma_delta = 0, mu = mu_deep, C_0 = C_0, cost_c = B_grid$cost_c, cost_P = 0, cost_M = 0, B = B_grid$B, measurement_error = FALSE) %>%
  bind_cols(B_grid) %>%
  mutate(M = "No Measurement Error") %>%
  mutate(depth = "Deep Soil (50-100 cm)")

optima_variance_deep <- bind_rows(optima_variance_dcea_deep, optima_variance_loi_deep, optima_variance_mirs_deep, optima_variance_no_error_deep) %>%
  mutate(cv = sqrt(optimum_variance) / mu_deep) # add coefficient of variation

################## plot optimal variances ##############
optima_variance <- bind_rows(optima_variance_top, optima_variance_deep) %>%
  mutate(depth = factor(depth, levels = c("Topsoil (0-10 cm)", "Deep Soil (50-100 cm)"))) %>%
  mutate(cost_c = paste(cost_c, "USD per core")) %>%
  mutate(cost_c = factor(cost_c, levels = paste(c(5,20,40), "USD per core")))

optima_std_error_plot <- ggplot(data = optima_variance, mapping = aes(x = B, y = sqrt(optimum_variance), color = as_factor(M))) +
  geom_line(size = 1.5) +
  labs(x = "Budget (USD)", y = "Standard Error", color = "Measurement") +
  scale_x_continuous(breaks = seq(500,5000,by=500)) +
  coord_cartesian(xlim = c(400,3000), ylim = c(0,.5)) +
  theme_bw() +
  theme(text = element_text(size = 16)) + 
  facet_grid(cost_c ~ depth)

optima_cv_plot <- ggplot(data = optima_variance, mapping = aes(x = B, y = cv, color = as_factor(M))) +
  geom_line(size = 1.5) +
  labs(x = "Budget (USD)", y = "Coefficient of Variation", color = "Measurement") +
  scale_x_continuous(breaks = seq(500,5000,by=500)) +
  coord_cartesian(xlim = c(400,3000), ylim = c(0,.4)) +
  theme_bw() +
  theme(text = element_text(size = 16)) + 
  facet_grid(cost_c ~ depth)

################ relative efficiencies ################ 
#relative efficiencies in the top soil
RE_dcea_loi_top <- get_relative_efficiency(sigma_p = sigma_p_top, mu = mu_top, cost_c = c(cost_c_low, cost_c_medium, cost_c_high), sigma_delta_1 = sigma_delta_dcea, cost_P_1 = cost_P_dcea, cost_M_1 = cost_M_dcea, sigma_delta_2 = sigma_delta_loi_top, cost_P_2 = cost_P_loi, cost_M_2 = cost_M_loi)
RE_dcea_mirs_top <- get_relative_efficiency(sigma_p = sigma_p_top, mu = mu_top, cost_c = c(cost_c_low, cost_c_medium, cost_c_high), sigma_delta_1 = sigma_delta_dcea, cost_P_1 = cost_P_dcea, cost_M_1 = cost_M_dcea, sigma_delta_2 = sigma_delta_mirs_top, cost_P_2 = cost_P_mirs, cost_M_2 = cost_M_mirs)
RE_mirs_loi_top <- get_relative_efficiency(sigma_p = sigma_p_top, mu = mu_top, cost_c = c(cost_c_low, cost_c_medium, cost_c_high), sigma_delta_1 = sigma_delta_mirs_top, cost_P_1 = cost_P_mirs, cost_M_1 = cost_M_mirs, sigma_delta_2 = sigma_delta_loi_top, cost_P_2 = cost_P_loi, cost_M_2 = cost_M_loi)

#relative efficiencies in deep soil
RE_dcea_loi_deep <- get_relative_efficiency(sigma_p = sigma_p_deep, mu = mu_deep, cost_c = c(cost_c_low, cost_c_medium, cost_c_high), sigma_delta_1 = sigma_delta_dcea, cost_P_1 = cost_P_dcea, cost_M_1 = cost_M_dcea, sigma_delta_2 = sigma_delta_loi_deep, cost_P_2 = cost_P_loi, cost_M_2 = cost_M_loi)
RE_dcea_mirs_deep <- get_relative_efficiency(sigma_p = sigma_p_deep, mu = mu_deep, cost_c = c(cost_c_low, cost_c_medium, cost_c_high), sigma_delta_1 = sigma_delta_dcea, cost_P_1 = cost_P_dcea, cost_M_1 = cost_M_dcea, sigma_delta_2 = sigma_delta_mirs_deep, cost_P_2 = cost_P_mirs, cost_M_2 = cost_M_mirs)
RE_mirs_loi_deep <- get_relative_efficiency(sigma_p = sigma_p_deep, mu = mu_deep, cost_c = c(cost_c_low, cost_c_medium, cost_c_high), sigma_delta_1 = sigma_delta_mirs_deep, cost_P_1 = cost_P_mirs, cost_M_1 = cost_M_mirs, sigma_delta_2 = sigma_delta_loi_deep, cost_P_2 = cost_P_loi, cost_M_2 = cost_M_loi)

RE_matrix <- cbind(rbind(RE_dcea_loi_top, RE_dcea_mirs_top), rbind(RE_dcea_loi_deep, RE_dcea_mirs_deep))
colnames(RE_matrix) <- c("Topsoil, 5 USD", "Topsoil, 20 USD", "Topsoil, 40 USD", "Deep Soil, 5 USD", "Deep Soil, 20 USD", "Deep Soil, 40 USD")
RE_matrix <- t(RE_matrix)

RE_table <- RE_matrix %>%
  round(digits = 2) %>%
  xtable()


######## minimum cost for a given variance #######
V_grid <- expand.grid(V = seq(.05, 2, by = .01)^2, cost_c = c(cost_c_low, cost_c_medium, cost_c_high))

optimal_costs_dcea_top <- get_minimum_cost(sigma_p = sigma_p_top, sigma_delta = sigma_delta_dcea, mu = mu_top, C_0 = C_0, cost_c = V_grid$cost_c, cost_P = cost_P_dcea, cost_M = cost_M_dcea, V = V_grid$V) %>%
  bind_cols(V_grid) %>%
  mutate(M = "DC-EA at 24.00 USD")


optimal_costs_loi_top <- get_minimum_cost(sigma_p = sigma_p_top, sigma_delta = sigma_delta_loi_top, mu = mu_top, C_0 = C_0, cost_c = V_grid$cost_c, cost_P = cost_P_loi, cost_M = cost_M_loi, V = V_grid$V) %>%
  bind_cols(V_grid) %>%
  mutate(M = "LOI at 4.25 USD")
optimal_costs_mirs_top <-  get_minimum_cost(sigma_p = sigma_p_top, sigma_delta = sigma_delta_mirs_top, mu = mu_top, C_0 = C_0, cost_c = V_grid$cost_c, cost_P = cost_P_mirs, cost_M = cost_M_mirs, V = V_grid$V) %>%
  bind_cols(V_grid) %>%
  mutate(M = "MIRS at 8.30 USD")


optimal_costs <- bind_rows(optimal_costs_dcea_top, optimal_costs_loi_top, optimal_costs_mirs_top) %>%
  mutate(composite_size = n/k)

optimal_cost_plot <- ggplot(optimal_costs, aes(x = sqrt(V), y = minimum_cost, color = as_factor(M))) +
  facet_grid(~ cost_c) +
  geom_line(size = 1.5) +
  coord_cartesian(xlim = c(.1,1.5), ylim = c(0,5000)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 200, linetype = "dashed") +
  labs(x = "Standard Error", y = "Minimum Budget (USD)", color = "Sample Cost") +
  scale_y_continuous(breaks = c(0,200,seq(500, 5000, by = 500))) +
  scale_x_continuous(breaks = c(.1,seq(0.25, 2.0, by = 0.25))) +
  theme_bw() +
  theme(text = element_text(size = 16)) 
