% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/external_resources.R
\name{getESET}
\alias{getESET}
\title{Create an ExpressionSet object}
\usage{
getESET(exprs, fdata, pdata)
}
\arguments{
\item{exprs}{A matrix or data frame of expression values. Rows correspond to
features (e.g., genes) and columns correspond to samples.}

\item{fdata}{A data frame or matrix of feature data. Each row corresponds to
feature in the expression set and columns correspond to feature variables or
a annotations. The row names should match the row names of the \code{exprs}
parameter.}

\item{pdata}{A data frame or matrix of phenotype data. Each row corresponds
to a sample in the expression set and columns correspond to phenotype
variables. The row names should match the column names of the \code{exprs}
parameter.}
}
\value{
An ExpressionSet object which includes the expression data (\code{exprs}),
phenotype data (\code{pdata}), and feature data (\code{fdata}).
}
\description{
This function creates an ExpressionSet object which is a container for
storing gene expression data along with related experimental data.
An ExpressionSet object includes an expression matrix, feature data (fdata)
and phenotype data (pdata).
}
