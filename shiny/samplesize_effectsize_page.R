#plot power effect size

basicPage(
  titlePanel("Return sample size needed to achieve a given power"),
  h5("This app computes the total number of samples necessary for a level-\\(\\alpha\\) test to detect a difference (at power \\(1 - \\beta\\)) between two plots with given means and standard deviations. The test is a two-sided, two-sample, unpaired t-test, with uniform independent random sampling, no compositing of cores. Option to include assay variability. Also assumes that %SOC is normally distributed or else that sample sizes are fairly large in order for the t-test to be approximately valid. Finally, this calculator assumes that the sample size is equally allocated to the two plots. In practice, sampling is more efficient when the sample size is allocated more to the plot with higher heterogeneity. Note that the sample size returned is the total for both sampling times."),
  sidebarPanel(
              numericInput("power_samplesize_power", "Input the desired power (\\(1 - \\beta\\))", value = 0.8),
              numericInput("alpha_samplesize_power", "Input the significance level of the test (\\(\\alpha\\))", value = 0.05),
             numericInput("mu_1_samplesize_power", "Input average %SOC in plot 1 (\\(\\mu_1\\))", value = 1),
             numericInput("mu_2_samplesize_power", "Input average %SOC in plot 2 (\\(\\mu_1\\))", value = 1.5),
             numericInput("sigma_p_1_samplesize_power", "Input spatial heterogeneity of %SOC in plot 1 (\\(\\sigma_{p1}\\))", value = 1.00),
             numericInput("sigma_p_2_samplesize_power", "Input spatial heterogeneity of %SOC in plot 2 (\\(\\sigma_{p2}\\))", value = 1.00),
             numericInput("sigma_delta_samplesize_power", "Input assay variability (\\(\\sigma_\\delta\\))", value = 0)
  ),
  mainPanel(
    h3("Results"),
    span(textOutput("n_effectsize"), style="color:blue; font-size:20px"),
    span(textOutput("n_samplesize"), style="color:blue; font-size:20px")
  )
)
