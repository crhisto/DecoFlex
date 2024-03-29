% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_guided_deconvolution.R
\name{performance_metrics_general}
\alias{performance_metrics_general}
\title{Calculate and Print Performance Metrics}
\usage{
performance_metrics_general(
  true_prop,
  calc_prop,
  title = NULL,
  verbose = FALSE
)
}
\arguments{
\item{true_prop}{A matrix or dataframe of the true proportions. Rows
correspond to cell types and columns to samples.}

\item{calc_prop}{A matrix or dataframe of the calculated proportions. It
should have the same structure as true_prop.}

\item{title}{Optional title for the print output.}

\item{verbose}{If TRUE, the performance metrics will be printed.
Default is FALSE.}
}
\value{
A list with three elements: 'pm_clusters' contains the performance
metrics by clusters, 'pm_samples' contains the performance metrics by
samples, true_prop' contains the true proportions, and 'calc_prop' contains
the calculated proportions.
}
\description{
This function computes various performance metrics to assess the accuracy of
calculated cell type proportions relative to the true proportions. The
performance metrics are computed both by clusters and by samples.
The results can be optionally printed.
}
\details{
The function computes performance metrics, such as mean absolute error,
root mean squared error, and correlation coefficient, using the
'calculate_performance_metrics' function. These metrics are calculated for
each cell type (i.e., by clusters) and for each sample.
If verbose is set to TRUE, these metrics will be printed with the
'print_performance_metrics' function.
}
