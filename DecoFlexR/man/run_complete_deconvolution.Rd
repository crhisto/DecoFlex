% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NMMFlex_wrapper.R
\name{run_complete_deconvolution}
\alias{run_complete_deconvolution}
\title{title Run Complete Deconvolution}
\usage{
run_complete_deconvolution(
  x_matrix,
  y_matrix = NULL,
  z_matrix = NULL,
  k,
  gamma = 1,
  alpha = 0,
  beta = 0,
  delta_threshold = 1e-10,
  max_iterations = 200,
  print_limit = 100,
  proportion_constraint_h = TRUE,
  regularize_w = NULL,
  alpha_regularizer_w = 0,
  fixed_w = NULL,
  fixed_h = NULL,
  fixed_a = NULL,
  fixed_b = NULL,
  initialized_w = NULL,
  initialized_h = NULL,
  initialized_a = NULL,
  initialized_b = NULL,
  init_method_w = "random_based.uniform",
  init_method_h = "random_based.uniform",
  init_method_a = "random_based.uniform",
  init_method_b = "random_based.uniform",
  partial_w_fixed = NULL,
  partial_h_fixed = NULL,
  w_mask_fixed = NULL,
  h_mask_fixed = NULL,
  batches_partial_fixed = 1,
  verbose = TRUE
)
}
\arguments{
\item{x_matrix}{A matrix, Input matrix corresponding to a dataframe.}

\item{y_matrix}{A matrix or NULL, Input matrix.}

\item{z_matrix}{A matrix or NULL, Input matrix.}

\item{k}{An integer, The rank used for deconvolution.}

\item{gamma}{A numeric, The gamma parameter value. Default is 1.}

\item{alpha}{A numeric, The alpha parameter value. Default is 0.0.}

\item{beta}{A numeric, The beta parameter value. Default is 0.0.}

\item{delta_threshold}{A numeric, The convergence threshold for stopping the
deconvolution iterations. Default is 1e-10.}

\item{max_iterations}{An integer, The maximum number of iterations to perform
during deconvolution. Default is 200.}

\item{print_limit}{An integer, The iteration interval at which to print
progress messages during deconvolution. Default is 100.}

\item{proportion_constraint_h}{A logical, Whether to apply a proportion
constraint to matrix H. Default is TRUE.}

\item{regularize_w}{A matrix or NULL, The regularization matrix for W.
Default is NULL.}

\item{alpha_regularizer_w}{A numeric, The alpha regularization parameter for
W. Default is 0.}

\item{fixed_w}{A matrix or NULL, The fixed matrix for W. Default is NULL.}

\item{fixed_h}{A matrix or NULL, The fixed matrix for H. Default is NULL.}

\item{fixed_a}{A matrix or NULL, The fixed matrix for A. Default is NULL.}

\item{fixed_b}{A matrix or NULL, The fixed matrix for B. Default is NULL.}

\item{initialized_w}{A matrix or NULL, The initial value for W. Default is
NULL.}

\item{initialized_h}{A matrix or NULL, The initial value for H. Default is
NULL.}

\item{initialized_a}{A matrix or NULL, The initial value for A. Default is
NULL.}

\item{initialized_b}{A matrix or NULL, The initial value for B. Default is
NULL.}

\item{init_method_w}{A string, The initialization method for W. Default is
'random_based.uniform'.}

\item{init_method_h}{A string, The initialization method for H. Default is
'random_based.uniform'.}

\item{init_method_a}{A string, The initialization method for A. Default is
'random_based.uniform'.}

\item{init_method_b}{A string, The initialization method for B. Default is
'random_based.uniform'.}

\item{partial_w_fixed}{A matrix or NULL, The partial fixed matrix for W.
Default is NULL.}

\item{partial_h_fixed}{A matrix or NULL, The partial fixed matrix for H.
Default is NULL.}

\item{w_mask_fixed}{A matrix or NULL, A mask matrix for W defining the
positions that are fixed with TRUE values. Default is NULL.}

\item{h_mask_fixed}{A matrix or NULL, A mask matrix for H defining the
positions that are fixed with TRUE values. Default is NULL.}

\item{batches_partial_fixed}{An integer, Defines the number of batches
where the procedure will inject again the the partially fixed matrices.
Default is 1.}

\item{verbose}{A logical, Whether to display progress messages during
deconvolution. Default is TRUE.}
}
\value{
The deconvoluted matrix.
}
\description{
Function that runs the src with three matrices as input. You
can define if the model is complete or partial in terms of the input and
therefore the output. For example, if some matrices are fixed (e.g., WH),
and Y is provided, then the variable Gamma = 0, resulting in Divergence(X|WG)
being zero, while the rest can still be calculated.

Deconvolution using three input matrices.
}
\details{
Raises an error if any of the input matrices.
Raises an error if the shape of the input matrices and the rank k is not
compatible.
Raises an error if any of the parameter values are invalid or not in the
expected range.
}
