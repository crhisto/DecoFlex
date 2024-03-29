% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{calculate_p_values_wt_vs_ko}
\alias{calculate_p_values_wt_vs_ko}
\title{Calculate p-values for WT vs KO}
\usage{
calculate_p_values_wt_vs_ko(results.prop)
}
\arguments{
\item{results.prop}{DataFrame, the results that will be used to perform the
t-tests.}
}
\value{
No return value. The function prints the result of the t-test for each cell
type.
}
\description{
This function calculates p-values for wild type (WT) vs knockout (KO) based
on given results. For more information, check:
https://www.cyclismo.org/tutorial/R/pValues.html
}
\details{
The function takes a DataFrame of results and transposes it. Then for each
cell type present in the data, it extracts the WT and KO results, performs a
t-test, and prints the result.
}
