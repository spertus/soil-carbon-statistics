library(shiny)


basicPage(
  titlePanel("Compute minimum budget for a fixed standard error"),
  h5("This app takes a standard error (estimation variance), costs, and parameters as inputs and returns the minimum budget that needs to estimate average SOC concentration. The number of cores (samples) and the assays to take is also returned.The currency used for costs and budget does not matter (results are scale invariant)."),
  sidebarPanel(
              numericInput("SE", "Input the desired estimation standard error (\\(SE\\))", value = .1),
              numericInput("C_0_budget", "Input total fixed costs (\\(C_0\\))", value = 0),
              numericInput("cost_c_budget", "Input the cost of taking a single core (\\(\\mbox{cost}_c\\))", value = 30),
              numericInput("cost_P_budget", "Input the cost of preparing a single composite sample (\\(\\mbox{cost}_P\\))", value = 5),
              numericInput("cost_A_budget", "Input the cost of assaying a single composite sample (\\(\\mbox{cost}_A\\))", value = 5),
              numericInput("sigma_p_budget", "Input the heterogeneity of the plot in %SOC (\\(\\sigma_{p}\\))", value = 1.00),
              numericInput("mu_budget", "Input the average plot concentration in %SOC (\\(\\mu\\))", value = 2.5),
              numericInput("sigma_delta_budget", "Input the assay uncertainty (\\(\\sigma_\\delta\\))", value = 0.25)
  ),
  textOutput("budget"),
  textOutput("n"),
  textOutput("k")
)
