library(shiny)


basicPage(
  titlePanel("Compute minimum standard error for a fixed budget"),
  h6("This app takes a budget, costs, and parmeters as inputs and returns the smallest achievable standard error to estimate average SOC concentration. The optimization is over the size of composites. The number of cores (samples) and assays to take is also returned. The currency used for costs and budget does not matter. The plot heterogeneity and mean are in percent SOC (grams SOC per hectograms of soil)."),
  sidebarPanel(
      numericInput("B", "Budget (\\(B\\))", value = 200),
      numericInput("C_0", "Input total fixed costs (\\(\\mbox{cost}_0\\))", value = 0),
      numericInput("cost_c", "Input the cost of taking a single core (\\(\\mbox{cost}_c\\))", value = 30),
      numericInput("cost_P", "Input the cost of preparing a single composite sample (\\(\\mbox{cost}_P\\))", value = 5),
      numericInput("cost_A", "Input the cost of assaying a single composite sample (\\(\\mbox{cost}_A\\))", value = 5),
      numericInput("sigma_p", "Input the plot heterogeneity in %SOC (\\(\\sigma_{p}\\))", value = 1.00),
      numericInput("mu", "Input average plot concentration in %SOC (\\(\\mu\\))", value = 2.5),
      numericInput("sigma_delta", "Input the assay uncertainty (\\(\\sigma_\\delta\\))", value = 0)
    ),
    mainPanel(
      h3("Results"),
      span(textOutput("error"), style="color:blue; font-size:20px"),
      span(textOutput("n_error"), style="color:blue; font-size:20px"),
      span(textOutput("k_error"), style="color:blue; font-size:20px")
    )
)



