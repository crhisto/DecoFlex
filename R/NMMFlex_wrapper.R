####################### DecoFlex Package: NMMFlex Module #######################
##
## This file is a module of the R package Decoflex. It provides a wrapper for
## the NMMFlex algorithm, which is used for flexible non-negative multiple
## matrix factorization NMMFlex.
##
## Included Files:
## 1. environment_creation.R: This file contains functions related to creating
##    the environment for NMMFlex.
## 2. run_standard_deconvolution.R: This file contains functions for running the
##    standard deconvolution using NMMFlex.
##
## Usage:
## - The environment_creation.R file is used to set up the necessary environment
##   for NMMFlex.
## - The run_standard_deconvolution.R file provides functions to perform
##   standard deconvolution using NMMFlex.
##
## Note:
## Please ensure that all required dependencies for Decoflex and NMMFlex are
## properly installed before using this module.
##

#' Set up environment for NMMFlex
#'
#' The environment_creation function installs and verifies the NMMFlex Python
#' package to prepare the required environment for using NMMFlex in R.
#'
#' @import reticulate
#'
#' @return This function does not return any value.
#'
#' @export
environment_creation <- function(){
  reticulate::py_install(
    packages="/Users/crhisto/Documents/GitHub/NMMFlex/NMMFlexPy",
    pip=TRUE)
  NMMFlex <- reticulate::import('NMMFlex')

  if(is.null(NMMFlex)){
    stop('The NMMFlex package is not loaded and perhaps not installed.')
  }else{
    message('The external Python package NMMFlex has been imported and therefore
            is installed')
  }
}

#' @title Run standard deconvolution
#'
#' @description This function is a wrapper for the MMNF library in python that
#' is used for running a standard deconvolution process. This involves
#' estimating the source distribution in a mixture model given some observed
#' data.
#'
#' @name run_standard_deconvolution
#'
#' @param bulk_data_x A numeric vector, matrix, or data frame representing the
#'  bulk data for the deconvolution.
#' @param references_w A numeric vector, matrix, or data frame representing the
#'  references for the deconvolution.
#' @param markers A numeric vector or NULL. Represents the markers to be used
#'  in the deconvolution. By default, this is set to NULL.
#' @param max_iterations A positive integer. Represents the maximum number of
#'  iterations for the deconvolution process. By default, this is set to 10000.
#' @param verbose A boolean. If TRUE, the function will print additional details
#'  during the execution. By default, this is set to FALSE.
#'
#' @return A list containing the results of the deconvolution. This might
#'  include the estimated sources and the final weights depending on the
#'  implementation of MMNF library in python.
#'
#' @import reticulate
#'
#' @export
run_standard_deconvolution <- function(bulk_data_x,
                                       references_w,
                                       markers = NULL,
                                       max_iterations=10000,
                                       verbose = FALSE){

  # 1. Calculation of the k equal to the number of cell references in the
  # reference matrix.
  k = ncol(references_w)

  # 2. I need to conciliate the number of genes in the bulk data and the
  # reference.
  genes_intersection <- intersect(rownames(bulk_data_x),
                                  rownames(references_w))
  bulk_data_x <- subset(bulk_data_x,
                        rownames(bulk_data_x) %in% genes_intersection)
  references_w <- references_w[genes_intersection, ]

  # 3. In case that I have marker genes.
  if(!is.null(markers)){
    #I take the total of genes shared between the samples
    bulk_data_x <- subset(bulk_data_x, rownames(bulk_data_x) %in% markers)
    references_w <- references_w[rownames(references_w) %in% markers, ]
  }

  # 4. Print the parameters
  if(verbose){
    message('NMMF parameters: ')
    message(paste0('markers: ', length(markers)))
    message(paste0('max_iterations: ', max_iterations))
    message(paste0('k: ', k))
    message(paste0('genes_intersection: ', length(genes_intersection)))
  }

  start_time <- Sys.time()

  # 5. Let's try the deconvolution method library. Import the NMMFlex python
  # library
  NMMFlex <- reticulate::import('NMMFlex')

  # For the reticulate library I have to specify the parameter type.
  deco_result <- NMMFlex$factorization$run_deconvolution_multiple(
    x_matrix = data.frame(bulk_data_x),
    y_matrix = NULL,
    z_matrix = NULL,
    k = as.integer(k), alpha = as.double(0.0), beta = as.double(0.0),
    delta_threshold = 1e-10,
    max_iterations = as.integer(max_iterations),
    proportion_constraint_h = as.logical(TRUE),
    fixed_w = data.frame(references_w),
    verbose = verbose)

  end_time <- Sys.time()
  end_time - start_time

  # 6. We return the result of the deconvolution.
  return(deco_result)
}
