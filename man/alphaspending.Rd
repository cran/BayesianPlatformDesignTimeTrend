% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plot_listoffunctionforplotgeneration.R
\name{alphaspending}
\alias{alphaspending}
\title{alphaspending}
\usage{
alphaspending(res)
}
\arguments{
\item{res}{A list of output matrix of a number of trial replicates}
}
\value{
The error rate at each interim analysis
}
\description{
This function estimates the mean error rate spent at each interim analysis for a trial
Example usage:
\enumerate{
\item sapply(res = result, fun = alphaspending) will generate list of the proportion of trial replicates are stopped at each stage for all scenarios in result where result is a list containing output data for different scenario
\item sapply(sapply(result,FUN = alphaspending),sum) will generate the type I error rate or power for all scenario on the result list
\item alpha(result) generate the proportion of trial replicates are stopped at each stage where result is the output data for one specific scenario
\item sum(alpha(result)) will generate the type I error rate or power for a specific scenario
}
}
\examples{
\dontrun{alphaspending(res)}
}
\author{
Ziyan Wang
}
