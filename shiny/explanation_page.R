library(shiny)


basicPage(
  h3("Explanation"),
  h6("This collection of shiny apps is designed to assist researchers with soil sampling design. The first two apps return sample and assay sizes that will minimize costs while achieving a certain precision or, vice versa, minimize precision while staying under a certain budget. The third app returns the power of a two-sample t-test to detect differences between two plots given sample and assay sizes. The 4th app graphs cost and precision across a range of possible composite sizes. The 5th and 6th apps plot the sample size and power (respectively) of a two-sample t-test, against the effect size."),
  
  h3("Notation"),
  h6("The apps refer to the following notation for plot parameters. Here are some slightly longform explanations of that notation. Notation for costs, budgets, and sample sizes also appears throughout the apps, but is fairly self-explanatory. For full details see [Spertus 2021 ArXiV link]."),
  h6("\\(\\cdot~~ \\mu \\): The plot average is the average %SOC over the plot. Formally it is defined as the integral \\( \\mu = \\int_{\\mathcal{P}} c(x,y) d\\mathcal{P} \\), where \\( \\mathcal{P} \\subset \\mathbb{R}^2 \\) is the set of points corresponding to the plot under study, and \\(c(x,y)\\) is the concentration at point \\( (x,y) \\in \\mathcal{P} \\)."),
  h6("\\(\\cdot~~ \\sigma_{p}\\): The plot heterogeneity is the average deviation of %SOC from the plot mean. Formally \\(\\sigma_{p} \\equiv \\sqrt{\\int_{\\mathcal{P}} (c(x,y) - \\mu)^2 d\\mathcal{P}} \\) where \\(\\mathcal{P} \\subset \\mathbb{R}^2 \\) is the plot under study and \\(c(x,y)\\) is the %SOC at point \\( (x,y) \\)."),
  h6("\\(\\cdot~~ \\sigma_{\\delta}\\): The assay error is the average deviation of the assay from the true SOC concentration. It is closely related to the average percent error in assay. For example, if \\(\\sigma_{\\delta} = 0.1\\) we would expect most assays of a single sample to be within about 10% of the true SOC concentration of that sample. If \\(\\sigma_{\\delta} = 0\\), there is no assay error. ")
  
)



