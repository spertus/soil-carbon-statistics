library(shiny)


basicPage(
  h3("Explanation"),
  h6("This collection of shiny apps is designed to assist researchers with soil sampling design.The apps can compute power for SOC change detection, recommend sample sizes, and maximize precision when compositing is used to reduce the number of assays. "),
  
  h3("Notation"),
  h6("The apps refer to the following notation for plot parameters:"),
  h6("\\(\\cdot~~ \\mu \\): The plot average is the average percent soil organic carbon (%SOC) over the plot."),
  h6("\\(\\cdot~~ \\sigma_{p}\\): The spatial heterogeneity is the standard deviation of %SOC over the plot."),
  h6("\\(\\cdot~~ \\sigma_{\\delta}\\): The assay error is the average deviation of the assay from the true SOC concentration. It is closely related to the average percent error in assay. For example, if \\(\\sigma_{\\delta} = 0.1\\) we would expect most assays of a single sample to be within about 10% of the true SOC concentration of that sample. If \\(\\sigma_{\\delta} = 0\\), there is no assay error."),
  h6("For more details refer to Spertus (2021)\\, available at:"), 
     tagList("https://www.scirp.org/journal/paperinformation?paperid=107467")
)



