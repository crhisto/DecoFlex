% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_guided_deconvolution.R
\name{calculate_min_correlation_incremental_markers}
\alias{calculate_min_correlation_incremental_markers}
\title{Calculate Minimum Correlation with Incremental Markers}
\usage{
calculate_min_correlation_incremental_markers(
  list_markers,
  reference,
  minimum_markers = 4,
  min_delta_cor_threshold = 0.05,
  verbose = TRUE,
  verbose_detailed = TRUE
)
}
\arguments{
\item{list_markers}{A list of marker genes for each group, ordered by log2
fold change in descending order.}

\item{reference}{The reference dataset to which the correlations will be
calculated.}

\item{minimum_markers}{The minimum number of marker genes to be selected,
default is set to 4. This can be adjusted based on the analysis
requirements.}

\item{min_delta_cor_threshold}{Minimum difference between the correlation of
t and t-1 steps. The idea is to skip genes that have small changes in
the correlation. For default it is a 5\%.}

\item{verbose}{If TRUE, prints detailed information during the function
execution. Default is TRUE.}

\item{verbose_detailed}{If TRUE, prints even more detailed information
during the function execution. Default is TRUE.}
}
\value{
A data frame with the minimum number of markers required to achieve
a satisfactory level of correlation.
}
\description{
This function processes a list of marker genes for each group, sorted in
descending order by log2 fold change. The objective is to determine the
number of genes that must be selected to achieve the minimum correlation
between samples in the reference dataset. The function is designed to
incrementally add markers, starting with those that display the greatest
change (as indicated by the log2 fold change). For each new set of markers,
the function calculates the correlation between samples in the reference
dataset.
}
\details{
The function begins by including the most significant markers (those with
the highest log2 fold change) in the analysis and calculates the correlation
between samples in the reference dataset. This process is repeated, each
time adding more markers, until a satisfactory level of correlation is
reached. This way, the function attempts to balance the inclusion of
informative markers against the potential noise introduced by less
significant markers.
}
