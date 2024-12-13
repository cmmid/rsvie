% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_economics_methods.R
\name{plot_economics_fit}
\alias{plot_economics_fit}
\title{Plot Economic Fit}
\usage{
plot_economics_fit(object)
}
\arguments{
\item{object}{An S4 object containing economic data and model outcomes.}
}
\value{
None. The function saves the generated plots to the specified directory.
}
\description{
This function generates and saves plots comparing fitted distributions to reference data for risks, QALYs, and costs across different age groups.
}
\details{
The function performs the following steps:
\itemize{
\item Checks if the directory for saving plots exists, and creates it if it doesn't.
\item Extracts unique age groups from the economic data.
\item Defines a helper function to extract mean and 95\% uncertainty intervals for model outcomes.
\item Extracts fitted model outcomes for risks, QALYs, and costs.
\item Prepares data frames for risks, QALYs, and costs by merging reference data with fitted model outcomes.
\item Defines a helper function to convert data frames to long format for plotting.
\item Converts the prepared data frames to long format.
\item Generates and saves plots comparing fitted model outcomes to reference data for risks, QALYs, and costs.
}
}
\examples{
\dontrun{
  plot_economics_fit(my_economic_model)
}
}
