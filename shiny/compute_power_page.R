
basicPage(
  titlePanel("Compute power for fixed sample sizes"),
  h5("Given means, spatial heterogeneities, and sample sizes for two plots (or sampling times), this app computes the power of a two-sample t-test where the null is no difference between average plot concentrations. Assumes independent simple random sampling from each plot. The validity of the test relies on SOC being approximately normally distributed in each plot."),
  sidebarPanel(
             numericInput("alpha_power", "Input the significance level of the test (\\(\\alpha\\))", value = 0.05),
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

