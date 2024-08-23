#plot power effect size

basicPage(
  titlePanel("Plot power across a range of effect sizes"),
    h5("Graph the power (y-axis) of a two-sample t-test against the effect size (x-axis), i.e. the difference in average SOC concentration between the two plots. The null hypothesis is no difference between average plot concentrations. Results assume simple random sampling, and that SOC is (marginally) normally distributed within plots. If SOC is not normally distributed results are approximately correct for large sample sizes."),
  sidebarPanel(
            numericInput("alpha_power_Delta", "Input significance level (\\(\\alpha\\))", value = 0.05),
            sliderInput("mu_2_power_Delta", label = "Input a range of average %SOC in plot 2 (\\(\\mu_2\\))", min = 0, max = 10, value = c(1, 4), step = 0.25),
            numericInput("mu_1_power_Delta", "Input average %SOC in plot 1 (\\(\\mu_1\\))", value = 1),
             numericInput("n1_power_Delta", "Input the number of cores from plot 1 (\\(n_1\\))", value = 30),
             numericInput("k1_power_Delta", "Input the number of assays of cores from plot 1 (\\(k_1\\))", value = 30),
             numericInput("n2_power_Delta", "Input the number of cores from plot 2 (\\(n_2\\))", value = 30),
             numericInput("k2_power_Delta", "Input the number of assays of cores from plot 2 (\\(k_2\\))", value = 30),
             numericInput("sigma_p_1_power_Delta", "Input heterogeneity of %SOC in plot 1 (\\(\\sigma_{p1}\\))", value = 1.00),
             numericInput("sigma_p_2_power_Delta", "Input heterogeneity of %SOC in plot 2 (\\(\\sigma_{p2}\\))", value = 1.00),
             numericInput("sigma_delta_power_Delta", "Input assay uncertainty (\\(\\sigma_\\delta\\))", value = 0)
  ),
  mainPanel(
    plotOutput("plot_power_Delta")
  )
)

# server <- function(input, output) {
# 
#   output$plot_power_Delta <- renderPlot({
#     mu_2_grid <- seq(input$min_mu_2_power_Delta, input$max_mu_2_power_Delta, length.out = 5000)
#     delta_grid <- mu_2_grid - input$mu_1_power_Delta
# 
#     solution <- get_power_two_sample(n_1 = input$n1_power_Delta, 
#                                      k_1 = input$k1_power_Delta, 
#                                      n_2 = input$n2_power_Delta, 
#                                      k_2 = input$k2_power_Delta, 
#                                      sigma_p_1 = input$sigma_p_1_power_Delta, 
#                                      mu_1 = input$mu_1_power_Delta, 
#                                      sigma_p_2 = input$sigma_p_2_power_Delta, 
#                                      mu_2 = mu_2_grid_power_Delta, 
#                                      sigma_delta = input$sigma_delta_power_Delta)
# 
# 
#     plot(x = delta_grid, y = solution, type = 'l', lwd = 1.5, xlab = "Effect size (difference in average %SOC in plot 1 versus plot 2)", ylab = "Power of two-sample t-test", ylim = c(0,1), xlim = c(min(delta_grid),max(delta_grid)))
#   })
# }
# 
# shinyApp(ui, server)