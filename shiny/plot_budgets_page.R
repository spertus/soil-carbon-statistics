library(shiny)


basicPage(
  titlePanel("Plot minimum budgets over a range of standard errors"),
  h5("Graph miminum budgets necessary to detect an effect of the standard error (x-axis) at a given power. Based on two-sided, two-sample t-test, with no compositing of samples."),
  sidebarPanel(
             numericInput("min_V", "Minimum Estimation Standard Error", value = .1),
             numericInput("max_V", "Maximum Estimation Standard Error", value = .2),
             numericInput("C_0_plot", "Fixed Cost", value = 0),
             numericInput("cost_P_plot", "Cost of Sample Prep", value = 5),
             numericInput("cost_A_plot", "Cost of Assay", value = 5),
             numericInput("cost_c_plot", "Cost of Sample", value = 30),
             numericInput("sigma_p_plot", "Sigma_p (plot sd)", value = 1.00),
             numericInput("mu_plot", "mu (plot mean concentration)", value = 2.5),
             numericInput("sigma_delta_plot", "Sigma_delta (measurement error sd)", value = 0.25)
  ),
  mainPanel(
    plotOutput("plot_budget")
  )
)
# 
# server <- function(input, output) {
# 
#   output$plot_budget <- renderPlot({
#     V_grid <- seq(input$min_V, input$max_V, length.out = 10000)^2
# 
#     cost_grid <- get_minimum_cost(V = V_grid, C_0 = input$C_0, cost_c = input$cost_c, cost_P = input$cost_P, cost_A = input$cost_M, sigma_p = input$sigma_p, mu = input$mu, sigma_delta = input$sigma_delta)$minimum_cost
# 
#     plot(x = sqrt(V_grid), y = cost_grid, type = 'l', lwd = 1.5, xlab = "Standard Error", ylab = "Minimum Cost (USD)")
#   })
# }
# 
# shinyApp(ui, server)
