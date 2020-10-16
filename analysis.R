#plots 
source("sim_functions.R")

mu <- 4.45
sigma_p <- 1.42
sigma_delta_dcea <- 0.14
sigma_delta_loi <- 0.45
C_0 <- 200
cost_c_low <- 5.00
cost_c_medium <- 15.00
cost_c_high <- 40.00
cost_P_loi <- 4.00
cost_P_dcea <- 6.00
cost_dcea <- 8.50
cost_loi_low <- 0.70
cost_loi_high <- 6.00

B_grid <- expand.grid(B = 200:2200, cost_c = c(cost_c_low, cost_c_medium, cost_c_high))


#from Vos 2005: LOI is estimated to be 1/12 the price of DC-EA.
#Vos estimates .20$ per LOI measurement and 2.35 for DC-EA
#root-mean squared predictive error (basically the measurement variance) is ~0.68% for LOI (when predicting DC-EA)

#optimal composite size for each method
optimal_composite_dcea <- get_optimal_composite_size(sigma_p = sigma_p, mu = mu, sigma_delta = sigma_delta_dcea, cost_c = c(cost_c_low, cost_c_medium, cost_c_high), cost_P = cost_P_dcea, cost_M = cost_dcea)
optimal_composite_loi_cheap <- get_optimal_composite_size(sigma_p = sigma_p, mu = mu, sigma_delta = sigma_delta_loi, cost_c = c(cost_c_low, cost_c_medium, cost_c_high), cost_P = cost_P_loi, cost_M = cost_loi_low)
optimal_composite_loi_expensive <- get_optimal_composite_size(sigma_p = sigma_p, mu = mu, sigma_delta = sigma_delta_loi, cost_c = c(cost_c_low, cost_c_medium, cost_c_high), cost_P = cost_P_loi, cost_M = cost_loi_high)


#what is the optimal variance that can be attained at different budgets?
optima_variance_dcea <- get_minimum_error(sigma_p = sigma_p, sigma_delta = sigma_delta_dcea, mu = mu, C_0 = C_0, cost_c = B_grid$cost_c, cost_P = cost_P_dcea, cost_M = cost_dcea, B = B_grid$B) %>%
  bind_cols(B_grid) %>%
  mutate(M = "DC-EA at 8.50 USD")

optima_variance_loi_cheap <- get_minimum_error(sigma_p = sigma_p, sigma_delta = sigma_delta_loi, mu = mu, C_0 = C_0, cost_c = B_grid$cost_c, cost_P = cost_P_loi, cost_M = cost_loi_low, B = B_grid$B) %>%
  bind_cols(B_grid) %>%
  mutate(M = "LOI at 0.70 USD")

optima_variance_loi_expensive <- get_minimum_error(sigma_p = sigma_p, sigma_delta = sigma_delta_loi, mu = mu,C_0 = C_0, cost_c = B_grid$cost_c, cost_P = cost_P_loi, cost_M = cost_loi_low, B = B_grid$B) %>%
  bind_cols(B_grid) %>%
  mutate(M = "LOI at 6.00 USD")

optima_variance <- bind_rows(optima_variance_dcea, optima_variance_loi_cheap) 

optima_variance_plot <- ggplot(data = optima_variance, mapping = aes(x = B, y = sqrt(optimum_variance), color = as_factor(cost_c))) +
  geom_line(size = 1.5) +
  labs(x = "Budget (USD)", y = "Standard Error", color = "Cost per Sample") +
  scale_x_continuous(breaks = seq(300,2000,by=200)) +
  coord_cartesian(xlim = c(300,2000), ylim = c(0,2)) +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  facet_grid(~ M) 

#plot n and k for DCEA
optima_variance_dcea_long <- optima_variance_dcea %>%
  select(n, k, B, optimum_variance, cost_c) %>%
  pivot_longer(c(n, k, optimum_variance), names_to = "quantity")
dcea_plot <- ggplot(data = optima_variance_dcea_long, aes(x = B, y = value, color = as_factor(cost_c))) +
  geom_line() +
  facet_wrap(~ quantity, scales = "free_y")

#add a line to compute relative efficiency of M_1 and P_1 vs M_2 and P_2





######## minimum cost for a given variance #######
V_grid <- expand.grid(V = seq(.01, 2, by = .01)^2, cost_c = c(cost_c_low, cost_c_medium, cost_c_high))

optimal_costs_dcea <- get_minimum_cost(sigma_p = sigma_p, sigma_delta = sigma_delta_dcea, mu = mu, C_0 = C_0, cost_c = V_grid$cost_c, cost_P = cost_P, cost_M = cost_dcea, V = V_grid$V) %>%
  bind_cols(V_grid) %>%
  mutate(M = "DC-EA")

optimal_costs_loi_cheap <- get_minimum_cost(sigma_p = sigma_p, sigma_delta = sigma_delta_loi, mu = mu, C_0 = C_0, cost_c = V_grid$cost_c, cost_P = cost_P, cost_M = cost_loi_low, V = V_grid$V) %>%
  bind_cols(V_grid) %>%
  mutate(M = "LOI")


optimal_costs <- optimal_costs_dcea %>% 
  bind_rows(optimal_costs_loi_cheap) %>%
  mutate(composite_size = n/k)

optimal_cost_plot <- ggplot(optimal_costs_dcea, aes(x = sqrt(V), y = minimum_cost, color = as_factor(cost_c))) +
  geom_line() +
  coord_cartesian(xlim = c(.25,2), ylim = c(0,2000)) +
  geom_hline(yintercept = 0) +
  labs(x = "Standard Error", y = "Minimum Cost", color = "Sample Cost") +
  theme_bw() 
