% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/programmeClass.R
\name{run_state,RSVProgramme-method}
\alias{run_state,RSVProgramme-method}
\title{A setMethod for the \code{run_state} function}
\usage{
\S4method{run_state}{RSVProgramme}(object, S_i = 1, direct = FALSE)
}
\arguments{
\item{object}{An object on which the state is to be run.}

\item{S_i}{An integer representing the state index. Default is 1.}

\item{direct}{A logical value indicating whether to run the state directly. Default is FALSE.}
}
\value{
a matrix with all the values of the state variables in the dynamic transmission model
}
\description{
A setMethod for the \code{run_state} function
}
