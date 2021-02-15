#plot power effect size

basicPage(
  titlePanel("Return sample size needed to achieve a given power"),
  h5("Return the number of samples necessary to detect a difference between two plots at a particular power. Assumes a two-sided, two-sample, unpaired t-test, with uniform independent random sampling, no compositing of cores, and a significance level of \\(\\alpha = .05\\). Also assumes that %SOC is normally distributed or else that sample sizes are fairly large in order for the t-test to be approximately valid."),
  sidebarPanel(
             numericInput("mu_1_samplesize_power", "Input average %SOC in plot 1 (\\(\\mu_1\\))", value = 1),
             numericInput("mu_2_samplesize_power", "Input average %SOC in plot 2 (\\(\\mu_1\\))", value = 1.5),
             numericInput("power_samplesize_power", "Input the desired power (\\(\\beta\\))", value = 0.8),
             numericInput("sigma_p_1_samplesize_power", "Input heterogeneity of %SOC in plot 1 (\\(\\sigma_{p1}\\))", value = 1.00),
             numericInput("sigma_p_2_samplesize_power", "Input heterogeneity of %SOC in plot 2 (\\(\\sigma_{p2}\\))", value = 1.00),
             numericInput("sigma_delta_samplesize_power", "Input assay uncertainty (\\(\\sigma_\\delta\\))", value = 0)
  ),
  mainPanel(
    h3("Results"),
    span(textOutput("n_effectsize"), style="color:blue; font-size:20px"),
    span(textOutput("n_samplesize"), style="color:blue; font-size:20px")
  )
)
