% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{choose_model_from_grid_search_object}
\alias{choose_model_from_grid_search_object}
\title{Choose Model from Grid Search Result}
\usage{
choose_model_from_grid_search_object(
  grid_search_result,
  alpha = NULL,
  beta = NULL
)
}
\arguments{
\item{grid_search_result}{A list containing the results of a grid search.
Each element of the list should represent one model and include the alpha
and beta parameters for that model.}

\item{alpha}{Numeric, the alpha parameter of the model to be returned.}

\item{beta}{Numeric, the beta parameter of the model to be returned.}
}
\value{
The model from 'grid_search_result' that matches the provided alpha and beta
parameters, or NULL if no such model exists.
}
\description{
This function retrieves a specific model from a grid search result based on
provided alpha and beta parameters.
}
\details{
The function iterates through the 'grid_search_result' list until it finds a
model that matches the provided alpha and beta parameters.
If a matching model is found, the function returns the model; otherwise, it
returns NULL.
}
