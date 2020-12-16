basicPage(
  titlePanel("Plot plots and samples"),
  h5("A fake field with samples."),
  sidebarPanel(
    h3("Plot parameters"),
    sliderInput("mu_sample", label = "Input an average %SOC for the plot (\\(\\mu\\))", min = 0, max = 10, value = 2, step = 0.25),
    numericInput("sigma_p_sample", "Input the maximum plot heterogeneity for plot 1 (\\(\\sigma_{p}\\))", value = 1.00),
    numericInput("range_sample", "Input the variogram range", value = 5),
    actionButton("refresh_plot_sample", "Refresh Plot"),
    h3("Sample parameters"),
    numericInput("n_sample", "Input a sample size (\\(n\\))", value = 30),
    selectInput("type_sample", "Input a sampling design", choices = c("Simple random sample" = "simple random sample", "Transect sample" = "transect", "Well-spread sample" = "well-spread"))
  ),
  mainPanel(
    plotOutput("plot_samples")
  )
)