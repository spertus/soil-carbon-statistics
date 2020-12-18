basicPage(
  titlePanel("Plot plots and samples"),
  h5("Generate a ficticious plot and draw samples. To intialize, hit the buttons to generate plot and generate sample. The slider will determine the population mean. The heterogeneity input determines the plot variance (\\(\\sigma_p\\)) but is not 1 to 1. Output graphs the plot and samples. Text output is the true (population) mean and standard deviation, along with the mean and standard deviation estimated from the samples."),
  sidebarPanel(
    h3("Plot parameters"),
    sliderInput("mu_sample", label = "Input an average %SOC for the plot (\\(\\mu\\))", min = 0, max = 10, value = 2, step = 0.25),
    numericInput("sigma_p_sample", "Input the plot heterogeneity for plot 1", value = 1.00),
    numericInput("range_sample", "Input the variogram range", value = 5),
    actionButton("refresh_plot_sample", "Generate Plot"),
    h3("Sample parameters"),
    numericInput("n_sample", "Input a sample size (\\(n\\))", value = 30),
    selectInput("type_sample", "Input a sampling design", choices = c("Simple random sample" = "simple random sample", "Transect sample" = "transect", "Well-spread sample" = "well-spread")),
    actionButton("refresh_samples_sample", "Generate Sample")
  ),
  mainPanel(
    plotOutput("plot_samples"),
    textOutput("population_mean"),
    textOutput("population_sd"),
    textOutput("sample_mean"),
    textOutput("sample_sd")
  )
)