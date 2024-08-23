#plot power effect size

basicPage(
  titlePanel("Plot sample size needed to achieve a given power"),
  h5("Graph sample sizes necessary to detect an effect of a particular size (x-axis) at a given power. Based on two-sided, two-sample t-test, with uniform independent random sampling, no compositing of cores, and a significance level of \\(\\alpha = .05\\). Also assumes SOC is normally distributed or else that sample sizes are fairly large (results for small sample sizes are suspect)."),
  sidebarPanel(
             numericInput("alpha_n_Delta", "Input the significance level (\\(\\alpha\\))", value = 0.05),
             numericInput("beta_n_Delta", "Input the desired power (\\(1 - \\beta\\))", value = 0.8),
             sliderInput("mu_2_n_Delta", label = "Input a range of average %SOC in plot 2 (\\(\\mu_2\\))", min = 0, max = 10, value = c(2, 4), step = 0.25),
             numericInput("mu_1_n_Delta", "Input average %SOC in plot 1 (\\(\\mu_1\\))", value = 1),
             numericInput("sigma_p_1_n_Delta", "Input heterogeneity of %SOC in plot 1 (\\(\\sigma_{p1}\\))", value = 1.00),
             numericInput("sigma_p_2_n_Delta", "Input heterogeneity of %SOC in plot 2 (\\(\\sigma_{p2}\\))", value = 1.00),
             numericInput("sigma_delta_n_Delta", "Input assay uncertainty (\\(\\sigma_\\delta\\))", value = 0)
  ),
  mainPanel(
    plotOutput("plot_n_Delta")
  )
)
