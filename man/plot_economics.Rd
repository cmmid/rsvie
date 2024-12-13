% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_economics_methods.R
\name{plot_economics}
\alias{plot_economics}
\title{Plot Economics Data}
\usage{
plot_economics(object)
}
\arguments{
\item{object}{An object containing economic data frames (\code{econ_df}, \code{risks_df}, \code{risks_vhr_df}) and an economic name (\code{econ_name}).}
}
\description{
This function generates and saves plots for economic data, risks data, and very high-risk data from the given object.
}
\details{
The function performs the following steps:
\itemize{
\item Checks if the directory for saving plots exists, and creates it if it does not.
\item Processes the economic data frame (\code{econ_df}) to create a combined data frame for QALY loss and Costs.
\item Generates and saves a plot for the economic data with uncertainty ribbons and lines.
\item Generates and saves a plot for the risks data with uncertainty ribbons and lines.
\item Generates and saves a plot for the very high-risk data with lines.
}
}
