library(shiny)
source("../functions.R")

####user interface####
ui <- fluidPage(
  withMathJax(),
  navbarPage("Estimation and Inference for Soil Organic Carbon",
             navbarMenu("App",
                        tabPanel("Explanation and Notation", source("explanation_page.R")$value),
                        tabPanel("Compute minimum standard error for a fixed budget",
                                 source("minimum_error_page.R")$value),
                        tabPanel("Compute minimum budget for a fixed standard error",
                                 source("minimum_budget_page.R")$value),
                        tabPanel("Compute power for fixed sample sizes",
                                 source("compute_power_page.R")$value),
                        tabPanel("Graph costs and standard errors across a range of composite sizes",
                                 source("plot_composites_page.R")$value),
                        tabPanel("Graph sample size needed to achieve a given power",
                                 source("plot_samplesize_effectsize_page.R")$value),
                        tabPanel("Graph power across a range of effect sizes",
                                 source("plot_power_effectsize_page.R")$value)
                        )
             )
           )

#### server ####
server <- function(input,output){
  
  #Minimum Error#
  output$error <- renderText({
    paste("Standard error of estimation:", round(get_minimum_error(B = input$B, C_0 = input$C_0, cost_c = input$cost_c, cost_P = input$cost_P, cost_A = input$cost_A, sigma_p = input$sigma_p, mu = input$mu, sigma_delta = input$sigma_delta)$optimum_variance, 2))
  })
  
  output$n_error <- renderText({
    paste("Sample size (number of cores):", round(get_minimum_error(B = input$B, C_0 = input$C_0, cost_c = input$cost_c, cost_P = input$cost_P, cost_A = input$cost_A, sigma_p = input$sigma_p, mu = input$mu, sigma_delta = input$sigma_delta)$n))
  })
  
  output$k_error <- renderText({
    paste("Number of assays:", round(get_minimum_error(B = input$B, C_0 = input$C_0, cost_c = input$cost_c, cost_P = input$cost_P, cost_A = input$cost_A, sigma_p = input$sigma_p, mu = input$mu, sigma_delta = input$sigma_delta)$k))
  })
  
  # Minimum Budget #
  output$budget <- renderText({
    solution <- get_minimum_cost(V = input$SE^2, C_0 = input$C_0_budget, cost_c = input$cost_c_budget, cost_P = input$cost_P_budget, cost_A = input$cost_A_budget, sigma_p = input$sigma_p_budget, mu = input$mu_budget, sigma_delta = input$sigma_delta_budget)
    
    paste("Minimum Budget:", round(solution$minimum_cost))
  })
  
  output$n <- renderText({
    solution <- get_minimum_cost(V = input$SE^2, C_0 = input$C_0_budget, cost_c = input$cost_c_budget, cost_P = input$cost_P_budget, cost_A = input$cost_A_budget, sigma_p = input$sigma_p_budget, mu = input$mu_budget, sigma_delta = input$sigma_delta_budget)
    
    paste("Sample size (number of cores):", round(solution$n))
  })
  
  output$k <- renderText({
    solution <- get_minimum_cost(V = input$SE^2, C_0 = input$C_0_budget, cost_c = input$cost_c_budget, cost_P = input$cost_P_budget, cost_A = input$cost_A_budget, sigma_p = input$sigma_p_budget, mu = input$mu_budget, sigma_delta = input$sigma_delta_budget)
    
    paste("Number of assays:", solution$k)
  })
  
  # Compute power 
  output$power <- renderText({
    
    solution <- round(get_power_two_sample(n_1 = input$n_1_power, k_1 = input$k_1_power, n_2 = input$n_2_power, input$k_2_power, sigma_p_1 = input$sigma_p_1_power, mu_1 = input$mu_1_power, sigma_p_2 = input$sigma_p_2_power, mu_2 = input$mu_2_power, sigma_delta = input$sigma_delta_power), 2)
    
    paste("Power:", solution)
  })
  
  
  # Minimum budget plot 
  output$plot_budget <- renderPlot({
    V_grid <- seq(input$min_V, input$max_V, length.out = 10000)^2
    
    solution <- get_minimum_cost(V = V_grid, C_0 = input$C_0_plot, cost_c = input$cost_c_plot, cost_P = input$cost_P_plot, cost_A = input$cost_A_plot, sigma_p = input$sigma_p_plot, mu = input$mu_plot, sigma_delta = input$sigma_delta_plot)
    
    costs <- solution$minimum_cost

    n_pos <- seq(min(V_grid), max(V_grid), length.out = 5)
    n_each <- solution$n[which(V_grid == n_pos)]

    plot(x = sqrt(V_grid), y = costs, type = 'l', lwd = 1.5, xlab = "Standard Error", ylab = "Minimum Cost")
    # text(x = n_pos, y = max(costs), labels = n_each)
  })

  
  # Plot Sample Size at a given power #
  output$plot_n_Delta <- renderPlot({
    mu_2_grid <- seq(input$mu_2_n_Delta[1], input$mu_2_n_Delta[2], length.out = 5000)
    delta_grid <- mu_2_grid - input$mu_1_n_Delta
    
    solution <- get_power_two_sample(sigma_p_1 = input$sigma_p_1_n_Delta, mu_1 = input$mu_1_n_Delta, sigma_p_2 = input$sigma_p_2_n_Delta, mu_2 = mu_2_grid, sigma_delta = input$sigma_delta_n_Delta, beta = 1 - input$power_n_Delta)
    
    
    plot(x = delta_grid, y = solution, type = 'l', lwd = 1.5, xlab = "Effect size (difference in average %SOC in plot 1 versus plot 2)", ylab = "Sample sizes needed from each plot", xlim = c(min(delta_grid),max(delta_grid)))
  })
  # Plot Power Effect Size#
  output$plot_power_Delta <- renderPlot({
    mu_2_grid <- seq(input$mu_2_power_Delta[1], input$mu_2_power_Delta[2], length.out = 5000)
    delta_grid <- mu_2_grid - input$mu_1_power_Delta
    
    solution <- get_power_two_sample(n_1 = input$n1_power_Delta, 
                                     k_1 = input$k1_power_Delta, 
                                     n_2 = input$n2_power_Delta, 
                                     k_2 = input$k2_power_Delta, 
                                     sigma_p_1 = input$sigma_p_1_power_Delta, 
                                     mu_1 = input$mu_1_power_Delta, 
                                     sigma_p_2 = input$sigma_p_2_power_Delta, 
                                     mu_2 = mu_2_grid, 
                                     sigma_delta = input$sigma_delta_power_Delta)
    
    
    plot(x = delta_grid, y = solution, type = 'l', lwd = 1.5, xlab = "Effect size (difference in average %SOC in plot 1 versus plot 2)", ylab = "Power of two-sample t-test", ylim = c(0,1), xlim = c(min(delta_grid),max(delta_grid)))
  })
  
  
  # Plot composites #
  output$plot_composites <- renderPlot({
    
    composite_grid <- get_composite_error_grid(n = input$n_composites, sigma_p = input$sigma_p_composites, sigma_delta = input$sigma_delta_composites, mu = input$mu_composites, C_0 = input$C_0_composites, cost_c = input$cost_c_composites, cost_A = input$cost_A_composites, cost_P = input$cost_P_composites) %>%
      pivot_longer(cols = c("std_error", "cost"), names_to = "quantity") %>%
      mutate(quantity = recode(quantity, std_error = "Standard Error", cost = "Cost"))
    
    
    ggplot(data = composite_grid, aes(x = composite_size, y = value)) +
      geom_line(linetype = "dashed") +
      geom_point(size = 3) +
      facet_grid(quantity ~ ., scales = "free") +
      xlab("Composite Size") +
      ylab("")
  })
  
  # plot a simulated plot and sampled points
  # surface <- eventReactive(input$refresh_plot_sample, {
  #   simulate_truth(size = c(50,50), nugget = 0, sill = (input$sigma_p_sample/100), range = input$range_sample, intercept = input$mu_sample/100, y_trend = FALSE) %>%
  #     mutate(z = 100*z)
  # })
  # 
  # 
  # samples <- eventReactive(input$refresh_samples_sample, {
  #   collect_sample(surface = surface(), design = input$type_sample, n_samp = input$n_sample)
  # })
  # 
  # output$plot_samples <- renderPlot({
  #   plot_surface_samples(surface = surface(), samples = samples())
  # })
  # output$population_mean <- renderText({
  #   paste("True (population) mean: ", round(mean(surface()$z), digits = 2), sep = "")
  # })
  # output$population_sd <- renderText({
  #   paste("True heterogeneity: ", round(sd(surface()$z), digits = 2), sep = "")
  # })
  # output$sample_mean <- renderText({
  #   paste("Sample mean: ", round(mean(samples()$z), digits = 2), sep = "")
  # })
  # output$sample_sd <- renderText({
  #   paste("Sample standard deviation: ", round(sd(samples()$z), digits = 2), sep = "")
  # })
  
  ## End ##
}

####Call####
shinyApp(ui, server)
