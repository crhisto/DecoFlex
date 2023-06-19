####################### DecoFlex Package: NMMFlex Module #######################
##
## This file is a module of the R package Decoflex. It provides a wrapper for
## the NMMFlex algorithm, which is used for flexible non-negative multiple
## matrix factorization NMMFlex.
##
## Included Files:
## 1. environment_creation: This contains procedures related to creating
##    the environment for NMMFlex.
## 2. run_standard_deconvolution: This contains the wrapper for running the
##    standard deconvolution using NMMFlex.
##
## Usage:
## - The environment_creation function is used to set up the necessary
##   environment
##   for NMMFlex.
## - The run_standard_deconvolution function provides the wrapper to perform
##   standard deconvolution using NMMFlex.
##
## Note:
## Please ensure that all required dependencies for Decoflex and NMMFlex are
## properly installed before using this module.
##


#' Load Required Libraries
#'
#' This function is designed to load the NMMFlex library, a Python package.
#' The function uses the `reticulate` package to facilitate interoperability
#' between R and Python, allowing the use of Python scripts within the R
#' environment. NMMFlex needs to be located at the specific path on your local
#' machine for this function to work.
#'
#' @return A message indicating whether or not the NMMFlex library was
#' successfully loaded.
#'
#' @details
#' This function first attempts to install the Python package from a local
#' directory using `reticulate::py_install`. It then tries to import the
#' package. If the import is successful, it confirms that the package is
#' installed and loaded. If not, it throws an error message. Make sure to
#' adjust the path to the location of the NMMFlex package on your machine.
#'
#' @note Ensure that the Python package is located at the provided path on your
#' machine before using this function.
#'
#' @import reticulate
#'
#' @seealso \code{\link{reticulate}}
#'
#' @export
load_libraries <- function(){
  # Download the NMMFlexPy code from www.github.com/crhisto/NMMFlex/NMMFlexPy
  reticulate::py_install(
    packages="[git_hub_folder]/NMMFlex/NMMFlexPy",
    pip=TRUE)
  NMMFlex <- reticulate::import('NMMFlex')

  if(is.null(NMMFlex)){
    stop('The NMMFlex package is not loaded and perhaps not installed.')
  }else{
    message('The external Python package NMMFlex has been imported and ',
            'therefore is installed')
  }
}

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
  # Download the NMMFlexPy code from www.github.com/crhisto/NMMFlex/NMMFlexPy
  reticulate::py_install(
    packages="[git_hub_folder]/NMMFlex/NMMFlexPy",
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
    message('NMMFlex parameters: ')
    message(paste0('markers: ', length(markers)))
    message(paste0('max_iterations: ', max_iterations))
    message(paste0('k: ', k))
    message(paste0('genes_intersection: ', length(genes_intersection)))
  }

  start_time <- Sys.time()

  # 5. Let's try the deconvolution method library. Import the NMMFlex python
  # library
  NMMFlex <- reticulate::import('NMMFlex')
  NMMFlex_factorization <- NMMFlex$factorization()

  # For the reticulate library I have to specify the parameter type.
  deco_result <- NMMFlex_factorization$run_deconvolution_multiple(
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
  return(trans_mult_deco_to_R(deco_result))
}

#' @title title Run Complete Deconvolution
#'
#' @description Function that runs the src with three matrices as input. You
#' can define if the model is complete or partial in terms of the input and
#' therefore the output. For example, if some matrices are fixed (e.g., WH),
#' and Y is provided, then the variable Gamma = 0, resulting in Divergence(X|WG)
#'  being zero, while the rest can still be calculated.
#'
#' @title Run Complete Deconvolution
#' @description Deconvolution using three input matrices.
#'
#' @param x_matrix A matrix, Input matrix corresponding to a dataframe.
#' @param y_matrix A matrix or NULL, Input matrix.
#' @param z_matrix A matrix or NULL, Input matrix.
#' @param k An integer, The rank used for deconvolution.
#' @param gamma A numeric, The gamma parameter value. Default is 1.
#' @param alpha A numeric, The alpha parameter value. Default is 0.0.
#' @param beta A numeric, The beta parameter value. Default is 0.0.
#' @param delta_threshold A numeric, The convergence threshold for stopping the
#'  deconvolution iterations. Default is 1e-10.
#' @param max_iterations An integer, The maximum number of iterations to perform
#'  during deconvolution. Default is 200.
#' @param print_limit An integer, The iteration interval at which to print
#'  progress messages during deconvolution. Default is 100.
#' @param proportion_constraint_h A logical, Whether to apply a proportion
#'  constraint to matrix H. Default is TRUE.
#' @param regularize_w A matrix or NULL, The regularization matrix for W.
#'  Default is NULL.
#' @param alpha_regularizer_w A numeric, The alpha regularization parameter for
#'  W. Default is 0.
#' @param fixed_w A matrix or NULL, The fixed matrix for W. Default is NULL.
#' @param fixed_h A matrix or NULL, The fixed matrix for H. Default is NULL.
#' @param fixed_a A matrix or NULL, The fixed matrix for A. Default is NULL.
#' @param fixed_b A matrix or NULL, The fixed matrix for B. Default is NULL.
#' @param initialized_w A matrix or NULL, The initial value for W. Default is
#'  NULL.
#' @param initialized_h A matrix or NULL, The initial value for H. Default is
#'  NULL.
#' @param initialized_a A matrix or NULL, The initial value for A. Default is
#'  NULL.
#' @param initialized_b A matrix or NULL, The initial value for B. Default is
#' NULL.
#' @param init_method_w A string, The initialization method for W. Default is
#'  'random_based.uniform'.
#' @param init_method_h A string, The initialization method for H. Default is
#'  'random_based.uniform'.
#' @param init_method_a A string, The initialization method for A. Default is
#'  'random_based.uniform'.
#' @param init_method_b A string, The initialization method for B. Default is
#'  'random_based.uniform'.
#' @param verbose A logical, Whether to display progress messages during
#'  deconvolution. Default is TRUE.
#'
#' @return The deconvoluted matrix.
#'
#' @details
#' Raises an error if any of the input matrices.
#' Raises an error if the shape of the input matrices and the rank k is not
#' compatible.
#' Raises an error if any of the parameter values are invalid or not in the
#' expected range.
#'
#' @export
run_complete_deconvolution <- function(x_matrix,
                                       y_matrix=NULL,
                                       z_matrix=NULL,
                                       k,
                                       gamma=1, alpha=0.0, beta=0.0,
                                       delta_threshold=1e-10,
                                       max_iterations=200,
                                       print_limit=100,
                                       proportion_constraint_h=TRUE,
                                       regularize_w=NULL,
                                       alpha_regularizer_w=0,
                                       fixed_w=NULL, fixed_h=NULL, fixed_a=NULL,
                                       fixed_b=NULL,
                                       initialized_w=NULL,
                                       initialized_h=NULL,
                                       initialized_a=NULL,
                                       initialized_b=NULL,
                                       init_method_w='random_based.uniform',
                                       init_method_h='random_based.uniform',
                                       init_method_a='random_based.uniform',
                                       init_method_b='random_based.uniform',
                                       verbose=TRUE){

  # 1. Let's try the deconvolution method library. Import the NMMFlex python
  # library
  NMMFlex <- reticulate::import('NMMFlex')
  NMMFlex_factorization <- NMMFlex$factorization()

  # 2. Call the function the runs the multiple deconvolution
  deco_result <- NMMFlex_factorization$run_deconvolution_multiple(
    x_matrix=x_matrix,
    y_matrix=y_matrix,
    z_matrix=z_matrix,
    k=as.integer(k),
    gamma=as.double(gamma),
    alpha=as.double(alpha),
    beta=as.double(beta),
    delta_threshold=as.double(delta_threshold),
    max_iterations=as.integer(max_iterations),
    print_limit=as.integer(print_limit),
    proportion_constraint_h=as.logical(proportion_constraint_h),
    regularize_w=regularize_w,
    alpha_regularizer_w=alpha_regularizer_w,
    fixed_w=fixed_w,
    fixed_h=fixed_h,
    fixed_a=fixed_a,
    fixed_b=fixed_b,
    initialized_w=initialized_w,
    initialized_h=initialized_h,
    initialized_a=initialized_a,
    initialized_b=initialized_b,
    init_method_w=init_method_w,
    init_method_h=init_method_h,
    init_method_a=init_method_a,
    init_method_b=init_method_b,
    verbose=as.logical(verbose))


  # 3. We return the result of the deconvolution.
  return(trans_mult_deco_to_R(deco_result))
}

#' @title Run Grid Search
#'
#' @description Performs a grid search in a parallelized manner over different
#' values of alpha and beta in the Non-negative matrix factorization. The search
#' is conducted over the combinations of given alpha and beta values.
#'
#' @param bulk_data_methylation A data.frame, Bulk data methylation matrix.
#' @param bulk_data_expression A data.frame, Bulk data expression matrix.
#' @param data_expression_auxiliary A data.frame, Auxiliary expression data
#'  matrix.
#' @param k An integer, Number of clusters.
#' @param alpha_list A list of numerics, List of alpha values to be considered
#'  in the grid search. If NULL, no grid search is performed over alpha values.
#' @param beta_list A list of numerics, List of beta values to be considered in
#'  the grid search. If NULL, no grid search is performed over beta values.
#' @param delta_threshold A numeric, The threshold value for convergence.
#'  Default is 1e-20.
#' @param max_iterations An integer, Maximum number of iterations for
#'  convergence. Default is 200.
#' @param print_limit An integer, Limit for print statements. Default is 100.
#' @param threads An integer, Number of CPU threads to be used. If 0, then it
#'  uses the total number of CPUs minus one. Default is 0.
#' @param proportion_constraint_h A logical, Whether to apply proportion
#'  constraint on H matrix. Default is TRUE.
#' @param regularize_w A logical, Whether to apply regularization on W matrix.
#'  If NULL, no regularization is applied.
#' @param alpha_regularizer_w_list A list of numerics, List of alpha regularizer
#'  values to be considered in the grid search. If NULL, no grid search is
#'  performed over alpha regularizer values.
#' @param fixed_w A matrix or NULL, Fixed matrix W. If NULL, it implies that W
#'  is not fixed.
#' @param fixed_h A matrix or NULL, Fixed matrix H. If NULL, it implies that H
#'  is not fixed.
#' @param fixed_a A matrix or NULL, Fixed matrix A. If NULL, it implies that A
#'  is not fixed.
#' @param fixed_b A matrix or NULL, Fixed matrix B. If NULL, it implies that B
#'  is not fixed.
#' @param model_type If the model is core expression or methylation, therefore
#'  Y and X assignation are going to change.
#'  Values: c('core_expression', 'core_methylation')
#'
#' @return A list of results from the grid search.
#'
#' @details
#' Raises an error if gamma, alpha, and beta are all 0, indicating that the
#' model cannot be executed.
#' Suggests the user to switch to the more direct function
#' run_deconvolution_multiple() if any two of gamma, alpha, beta are zero.
#'
#' @export
run_grid_search <- function(bulk_data_methylation,
                            bulk_data_expression=NULL,
                            data_expression_auxiliary=NULL,
                            k,
                            alpha_list=NULL,
                            beta_list=NULL,
                            delta_threshold=1e-20,
                            max_iterations=200,
                            print_limit=100,
                            threads=0,
                            proportion_constraint_h=TRUE,
                            regularize_w=NULL,
                            alpha_regularizer_w_list=NULL,
                            fixed_w=NULL, fixed_h=NULL,
                            fixed_a=NULL, fixed_b=NULL,
                            model_type){

  # 1. Let's try the deconvolution method library. Import the NMMFlex python
  # library
  NMMFlex <- import('NMMFlex')
  NMMFlex_grid_search <- NMMFlex$grid_search()

  # 2. Initialization of
  x_value <-  NULL
  y_value <- NULL

  # 3. Depend on the core model, I send the methylation or the expression as
  # core or as a secondary matrix.
  if(model_type == 'core_methylation'){
    x_value <-  bulk_data_methylation
    y_value <- bulk_data_expression
  }else if(model_type == 'core_expression'){
    x_value <-  bulk_data_expression
    y_value <- bulk_data_methylation
  }else{
    # To avoid problems in old code without the new parameter
    stop(
      "The 'model_type' variable holds an incorrect value. Please ensure that",
      " it corresponds to either 'core_methylation' or 'core_expression' to ",
      "signify the model type, which should be either methylation core or ",
      "expression core. The current value is: ", model_type)
  }

  # 4. Call the function the runs the multiple deconvolution
  deco_grid_search_result <-
    NMMFlex_grid_search$grid_search_parallelized_alpha_beta(
      bulk_data_methylation=x_value,
      bulk_data_expression=y_value,
      data_expression_auxiliary=data_expression_auxiliary,
      k=as.integer(k),
      alpha_list=alpha_list,
      beta_list=beta_list,
      delta_threshold=as.double(delta_threshold),
      max_iterations=as.integer(max_iterations),
      print_limit=as.integer(print_limit),
      threads=as.integer(threads),
      proportion_constraint_h=as.logical(proportion_constraint_h),
      regularize_w=regularize_w,
      alpha_regularizer_w_list=alpha_regularizer_w_list,
      fixed_w=fixed_w,
      fixed_h=fixed_h,
      fixed_a=fixed_a,
      fixed_b=fixed_b)

  # 5. We return the result of the deconvolution.
  return(trans_grid_search_to_R(deco_grid_search_result))
}

#' @title Translate Multiple Deconvolution Results to R object
#'
#' @description This function takes the output of a multiple deconvolution
#' process and translates it into an R-friendly format, returning a list
#' of all relevant results. The idea is to preserve the python data in the R
#' allocated memory, otherwise under certain circumstance the object will be
#' null.
#'
#' @param multiple_deconvolution_results A list, The result object returned 
#' from a multiple deconvolution function.
#'
#' @return A list, with the following elements:
#'   "x", "y", "z": The original input matrices.
#'   "x_hat", "y_hat", "z_hat": The estimated matrices after 
#'    deconvolution.
#'   "w", "h": The matrices of factor loadings and factor scores.
#'   "a", "b": Additional matrices related to the deconvolution model.
#'   "is_x_sparse", "is_y_sparse", "is_z_sparse", "is_model_sparse": 
#'    Boolean values indicating if the corresponding matrices are sparse.
#'   "initialized_w", "initialized_h", "initialized_a", "initialized_b": 
#'    The initial values used for the corresponding matrices in the 
#'    deconvolution process.
#'   "iterations": The number of iterations performed during the 
#'    deconvolution process.
#'   "divergence_value", "delta_divergence_value": The final divergence 
#'    value and the change in divergence in the last iteration.
#'   "running_info": Additional information about the running process.
#'   "alpha", "beta", "alpha_regularizer_w": The parameters used in the 
#'    deconvolution model.
#'
#' @details This function is designed to simplify the handling of multiple 
#' deconvolution results by translating the returned object from Python into 
#' an R-friendly format.
#'
#' @export
trans_mult_deco_to_R <- function(multiple_deconvolution_results){
  
  return_list <- list(
    'x' = multiple_deconvolution_results$x,
    'y' = multiple_deconvolution_results$y,
    'z' = multiple_deconvolution_results$z, 
    
    'x_hat' = multiple_deconvolution_results$x_hat,
    'y_hat' = multiple_deconvolution_results$y_hat,
    'z_hat' = multiple_deconvolution_results$z_hat,
    
    'w' = multiple_deconvolution_results$w,
    'h' = multiple_deconvolution_results$h,
    
    'a' = multiple_deconvolution_results$a,
    'b' = multiple_deconvolution_results$b,
    
    'is_x_sparse' = multiple_deconvolution_results$is_x_sparse,
    'is_y_sparse' = multiple_deconvolution_results$is_y_sparse,
    'is_z_sparse' = multiple_deconvolution_results$is_z_sparse,
    'is_model_sparse' = multiple_deconvolution_results$is_model_sparse,
    
    'initialized_w' = multiple_deconvolution_results$initialized_w,
    'initialized_h' = multiple_deconvolution_results$initialized_h,
    'initialized_a' = multiple_deconvolution_results$initialized_a,
    'initialized_b' = multiple_deconvolution_results$initialized_b,
    
    'iterations' = multiple_deconvolution_results$iterations,
    'divergence_value' = multiple_deconvolution_results$divergence_value,
    'delta_divergence_value' = 
      multiple_deconvolution_results$delta_divergence_value,
    'running_info' = multiple_deconvolution_results$running_info,
    
    'alpha' = multiple_deconvolution_results$alpha,
    'beta' = multiple_deconvolution_results$beta,
    'alpha_regularizer_w' = multiple_deconvolution_results$alpha_regularizer_w
  )
  
  return(return_list)
}

#' @title Translate Grid Search Results to R
#'
#' @description This function takes the output of a grid search and translates 
#' it into an R-friendly format, returning a list
#' of all relevant results. It works by iterating over the grid search results,
#' converting those into native R objects.The idea is to preserve the python 
#' data in the R allocated memory, otherwise under certain circumstance the 
#' object will be null.
#'
#' @param grid_search_result A list. This is the result object returned from a 
#'  grid search process.
#' @param verbose A logical. If TRUE, the function will print information about
#'  the progress of the transformation for each model in the console. Default 
#'  is FALSE.
#'
#' @return A list of all the elements contained in grid_search_result, 
#'  transformed into an R-friendly format. Each element in the list corresponds
#'  to a model from the grid search and is named "alpha_beta_" followed by the 
#'  alpha and beta parameters of the model.
#'
#' @details This function is designed to simplify the handling of grid search 
#' results by translating the returned objects from Python into an R-friendly 
#' format. The function uses \code{trans_mult_deco_to_R()} to transform each 
#' individual deconvolution result.
#'
#' @export
trans_grid_search_to_R <- function(grid_search_result,
                                   verbose = FALSE){
  
  #Iteration over the results, converting those to native R objects
  return_list <- list()
  for(i in 1:length(grid_search_result)){
    if(verbose){
      print(paste0('Model: ',i, '') )
    }
    
    # Getting each model
    object_deco <- grid_search_result[[i]]$get()
    
    # Transforming the object of the specific deconvolution.
    object_deco_transformed <- trans_mult_deco_to_R(object_deco)
    
    # Adding it to the list.
    return_list[[paste0('alpha_beta_',
                        object_deco$alpha,
                        '_',
                        object_deco$beta)]] <- object_deco_transformed
  }
  
  return(return_list)
}