

basicPage(
  titlePanel("Plot costs and standard errors across a range of composite sizes"),
  h5("Graph costs and standard errors in %SOC (on the y-axes) against size of composites (on the x-axis). Compositing allows investigators to conduct fewer assays of a given number of samples (\\(n\\)). Compositing thus lowers costs but raises error when there is assay error. These plots allow us to explore the tradeoff."),
  sidebarPanel(
             numericInput("n_composites", "Input the number of cores (\\(n\\))", value = 100),
             numericInput("C_0_composites", "Input total fixed costs (\\(\\mbox{cost}_0\\))", value = 0),
             numericInput("cost_c_composites", "Input the cost of taking a single core (\\(\\mbox{cost}_c\\))", value = 20),
             numericInput("cost_P_composites", "Input the cost of preparing a single composite sample (\\(\\mbox{cost}_P\\))", value = 5),
             numericInput("cost_A_composites", "Input the cost of assaying a single composite sample (\\(\\mbox{cost}_A\\))", value = 5),
             numericInput("sigma_p_composites", "Input heterogeneity of %SOC in plot 1 (\\(\\sigma_{p1}\\))", value = 1.00),
             numericInput("mu_composites", "Input the average plot concentration in %SOC (\\(\\mu\\))", value = 2.5),
             numericInput("sigma_delta_composites", "Input the assay uncertainty (\\(\\sigma_\\delta\\))", value = 0.25)
  ),
  mainPanel(
    plotOutput("plot_composites")
  )
)

