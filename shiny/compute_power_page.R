
basicPage(
  titlePanel("Compute power for fixed sample sizes"),
  h5("Given plot parameters and sample sizes, this app computes the power of a two-sample t-test where the null is no difference between average plot concentrations. Assumes uniform independent random sampling from the plot. Note that the validity of the test hinges on SOC being (marginally) normally distributed within the plot. If sample sizes are large, the test should be approximately correct by the central limit theorem."),
  sidebarPanel(
             numericInput("mu_1_power", "Input true average %SOC in plot 1 (\\(\\mu_1\\))", value = 1),
             numericInput("mu_2_power", "Input true average %SOC in plot 2 (\\(\\mu_2\\))", value = 2),
             numericInput("n_1_power", "Input the number of cores from plot 1 (\\(n_1\\))", value = 30),
             numericInput("k_1_power",  "Input the number of assays of cores from plot 1 (\\(k_1\\))", value = 30),
             numericInput("n_2_power", "Input the number of cores from plot 2 (\\(n_2\\))", value = 30),
             numericInput("k_2_power",  "Input the number of assays of cores from plot 2 (\\(k_2\\))", value = 30),
             numericInput("sigma_p_1_power", "Input heterogeneity of %SOC in plot 1 (\\(\\sigma_{p1}\\))", value = 1.00),
             numericInput("sigma_p_2_power", "Input heterogeneity of %SOC in plot 2 (\\(\\sigma_{p2}\\))", value = 1.00),
             numericInput("sigma_delta_power", "Input assay uncertainty (\\(\\sigma_\\delta\\))", value = 0)
  ),
  mainPanel(
    h3("Results"),
    span(textOutput("power"), style="color:blue; font-size:20px")
  )
)

