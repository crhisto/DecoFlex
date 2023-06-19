#################  Utility Functions for the DecoFlex Package #################
##
## This script contains utility functions used throughout the DecoFlex package.
## The functions included are as follows:
##
## - `create_dataframe_results_multiple_deconvolutions()`: Creates a dataframe
##    with the summary of the gridsearch procedure.
##
## - `performance_plot_multiple_deconvolutions()`: Plots the performance of
##    different deconvolutions based on given parameters.
##
## - `plot_heatmap_alpha_beta()`: Plots the divergence value for each alpha,
##    beta value to identify the optimal parameters.
##
## - `choose_model_from_grid_search_object()`: Returns a specific model from
##    the gridSearch object based on the alpha and beta values.
##
## - `save_plot_as_pdf()`: Saves a given plot as a PDF in a specified path.
##
## - `plot_performance_for_each_model()`: Plots the performance for each model
##    generated from the grid search.
##
## - `plot_deconvolution_performance_x_y_z_total()`: Detailed plots for each
##    model showing the divergence of x, y, z and total.
##
## - `calculate_sparseness_matrix()`: Calculates the sparseness of a given
##    matrix in percentages.
##
## - `order_result_generic()`: Orders the result matrix based on specified
##    starting column and range.
##
## - `order_simple_vector()`: Orders a simple vector based on specified starting
##    index and range.
##
## - `calculate_p_values_wt_vs_ko()`: Calculates the p-values for WT vs KO from
##    the provided result prop.
##
## @note All the functions implemented in this script assume that the user
##       provides appropriate and correctly formatted inputs.
##       The error handling within each function is specific to the expected
##       input types, and may not function as expected for all incorrect inputs.
##


#' Plot Deconvolution Performance Metrics
#'
#' This function plots various performance metrics for a deconvolution process
#' across multiple iterations. Function for gridSearch analysis: General plots.
#'
#' @param running_info A data frame with iteration number and the corresponding
#'  values for divergence and delta divergence.
#' @param title_addition An optional string that is added to the plot title.
#'
#' @details
#' This function creates three plots to visualize the performance of a
#' deconvolution process:
#' 1. The natural logarithm of the divergence value over iterations.
#' 2. The natural logarithm of the delta divergence value over iterations.
#' 3. Both the above metrics in the same plot.
#' The divergence value provides a measure of how much the deconvolution results
#' deviate from the ground truth, while the delta divergence gives a measure of
#' how much the divergence changes between successive iterations.
#'
#' @return
#' The function returns a list of ggplot objects, each representing a different
#' visualization of the deconvolution performance.
#'
#' @import ggplot2
#' @importFrom ggplot2 .data
#'
#'
# @export
plot_deconvolution_performance <- function(running_info, title_addition = ''){
  colnames(running_info) <- c('iteration', 'divergence_value_x',
                              'divergence_value_y', 'divergence_value_z',
                              'divergence_value', 'delta_divergence_value')
  rownames(running_info) <- seq(1:nrow(running_info))
  running_info.df <- data.frame(running_info)

  # Solution for the error 'no visible binding for global variable' posted in:
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  # and solved in https://stackoverflow.com/questions/9439256 using:
  # @importFrom ggplot2 .data, and the .data$ prefix.

  # Iterations vs divergence value
  ggplot2::ggplot(running_info.df,
                  ggplot2::aes(.data$iteration, log(.data$divergence_value))) +
    ggplot2::geom_line() + ggplot2::theme_bw() +
    ggplot2::ggtitle("Iterations vs log(divergence)")

  # Iterations vs delta
  ggplot2::ggplot(running_info.df,
                  ggplot2::aes(
                    .data$iteration, log(.data$delta_divergence_value))) +
    ggplot2::geom_line() + ggplot2::theme_bw() +
    ggplot2::ggtitle("Iterations vs log(delta divergence)")

  colors_by_lines <- c("divergence_value"="black",
                       "delta_divergence_value"="red")

  # Plot with both metric included
  ggplot2::ggplot(running_info.df, ggplot2::aes(x=.data$iteration)) +
    ggplot2::geom_line(ggplot2::aes(y = log(.data$divergence_value),
                  color = 'divergence_value')) +
    ggplot2::geom_line(ggplot2::aes(y = log(.data$delta_divergence_value),
                  color = 'delta_divergence_value')) +
    ggplot2::scale_color_manual(name = "Metric",
                       values=colors_by_lines) +
    ggplot2::ylab('log(metric)') +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(paste0("Iterations vs log(divergence)",
                            " and log(delta divergence)",
                            ' - ', title_addition)) +
    ggplot2::theme(plot.title = ggplot2::element_text(size=10))
}

#' Create Data Frame of Multiple Deconvolution Results
#'
#' This function creates a summary data frame from the results of multiple
#' deconvolutions, typically obtained from a grid search procedure.
#'
#' @param grid_info A list of deconvolution results obtained from a grid search.
#' @param verbose Shows details about the added deconvolutions.
#'
#' @details
#' The function iterates through all the deconvolution results in the input list.
#' For each result, it extracts relevant performance metrics such as alpha, beta,
#' divergence value, and delta divergence value. These metrics are then combined
#' into a single data frame. Additionally, it computes the logarithm of the
#' divergence and delta divergence values, and includes them in the data frame
#' as well.
#'
#' @return
#' A data frame that contains performance metrics for each deconvolution result,
#' including alpha, beta, divergence value, delta divergence value, and their
#' logarithmic forms.
#'
#' @export
create_dataframe_results_multiple_deconvolutions <- function(grid_info,
                                                             verbose=TRUE){

  performance_df <- NULL

  # Extracting information from the models.
  for(i in 1:length(grid_info)){
    if(verbose){
      print(paste0('Model: ', i))
    }

    object_deco <- grid_info[[i]]
    performance_df <- rbind(performance_df,
                            c(object_deco$alpha,
                              object_deco$beta,
                              object_deco$divergence_value,
                              object_deco$delta_divergence_value))

    if(verbose){
      print("Running information:")
    }
  }

  colnames(performance_df) <- c('alpha', 'beta', 'divergence_value',
                                'delta_divergence_value')
  rownames(performance_df) <- performance_df[,1]

  performance_df <- as.data.frame(performance_df)

  # Create log forms
  performance_df$log_divergence_value <- log(performance_df$divergence_value)
  performance_df$log_delta_divergence_value <-
    log(performance_df$delta_divergence_value)

  performance_df
}

#' Plot Performance of Multiple Deconvolutions
#'
#' This function generates a plot showing the performance of multiple
#' deconvolutions. The plot visualizes the divergence and delta divergence
#' metrics against a main parameter of interest. Besides it is used for plotting
#' the results for the gridSearch in a general level.
#'
#' @param performance_df A data frame that contains performance metrics for each
#'  deconvolution result. The data frame should include the main parameter,
#'  divergence, and delta divergence values.
#' @param main_parameter A string that specifies the name of the main parameter
#'  to be plotted on the x-axis.This parameter should be a column in the
#'  'performance_df' data frame.
#'
#' @details
#' The function uses ggplot2 to generate a line plot where the y-axis represents
#'  the logarithmic values of the divergence and delta divergence metrics,
#'  and the x-axis represents the main parameter. The divergence and delta
#'  divergence metrics are differentiated by color.
#'
#' @return
#' A ggplot object showing the performance of multiple deconvolutions against
#' the specified main parameter.
#'
#' @import ggplot2
#' @import_from ggplot2 .data
#'
#'
#' @export
performance_plot_multiple_deconvolutions <- function(performance_df,
                                                     main_parameter = 'alpha'){

  # Solution for the error 'no visible binding for global variable' posted in:
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  # and solved in https://stackoverflow.com/questions/9439256 using:
  # @importFrom ggplot2 .data, and the .data$ prefix.

  # Plot with both metric included
  ggplot2::ggplot(performance_df, ggplot2::aes_string(x=main_parameter)) +
    ggplot2::geom_line(ggplot2::aes(y = log(.data$divergence_value),
                  color = 'divergence')) +
    ggplot2::geom_line(ggplot2::aes(y = log(.data$delta_divergence_value),
                  color = 'delta_divergence')) +
    ggplot2::scale_color_manual(name = "Metric",
                       values=c("blue", "green"),
                       labels = c("divergence", "delta_divergence")) +
    ggplot2::ylab('log(metric)') +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(paste0(main_parameter,
                            " vs (divergence and delta divergence)"))
}

#' Plot Heatmap of Alpha and Beta Values
#'
#' This function generates a heatmap showing the divergence value for each
#' combination of alpha and beta values. This visualization aids in identifying
#' the optimal alpha and beta parameters.
#'
#' @param performance_df A data frame that contains performance metrics for each
#'  alpha, beta combination. The data frame should include columns for alpha,
#'  beta, and divergence value.
#' @param fill_value A string indicating the column to be used for filling the
#'  heatmap. Default is "divergence_value".
#' @param title A string specifying the title of the plot. Default is
#'  'Divergence performance: alpha vs beta'.
#'
#' @details
#' The function uses ggplot2 to generate a heatmap where the x-axis represents
#' alpha values and the y-axis represents beta values. The color fill of each
#' tile is determined by the 'fill_value' parameter, which is typically set to
#' 'divergence_value'.
#'
#' @return
#' A ggplot object showing a heatmap of divergence values for different
#' combinations of alpha and beta values.
#'
#' @import ggplot2
#' @import grDevices
#' @import RColorBrewer
#'
#' @export
plot_heatmap_alpha_beta <- function(
    performance_df,
    fill_value="divergence_value",
    title = 'Divergence performance: alpha vs beta'){

  performance_df.plot <- performance_df

  performance_df.plot$alpha <- factor(performance_df$alpha)
  performance_df.plot$beta <- factor(performance_df$beta)

  color_combination_gradient <-
    grDevices::colorRampPalette(RColorBrewer::brewer.pal(10, "RdYlBu"))(256)

  ggplot2::ggplot(performance_df.plot,
         ggplot2::aes_string(x='alpha', y='beta', fill=fill_value)) +
    ggplot2::ggtitle(title) +
    ggplot2::geom_tile(color = "black") +
    ggplot2::theme_light() +
    ggplot2::scale_fill_gradientn(colours = color_combination_gradient) +
    ggplot2::coord_fixed() +
    ggplot2::guides(fill = ggplot2::guide_colourbar(barwidth = 0.5,
                                  barheight = 10))
}

#' Choose Model from Grid Search Result
#'
#' This function retrieves a specific model from a grid search result based on
#' provided alpha and beta parameters.
#'
#' @param grid_search_result A list containing the results of a grid search.
#'  Each element of the list should represent one model and include the alpha
#'  and beta parameters for that model.
#' @param alpha Numeric, the alpha parameter of the model to be returned.
#' @param beta Numeric, the beta parameter of the model to be returned.
#'
#' @details
#' The function iterates through the 'grid_search_result' list until it finds a
#' model that matches the provided alpha and beta parameters.
#' If a matching model is found, the function returns the model; otherwise, it
#' returns NULL.
#'
#' @return
#' The model from 'grid_search_result' that matches the provided alpha and beta
#' parameters, or NULL if no such model exists.
#'
#' @export
choose_model_from_grid_search_object <- function(grid_search_result,
                                                 alpha=NULL, beta=NULL){
  return_object = NULL

  # let's plot model 0.001 and 1 beta and alpha 0.
  for(i in 1:length(grid_search_result)){
    print(paste0('Model: ',i, '') )
    object_deco <- grid_search_result[[i]]

    if(object_deco$alpha==alpha && object_deco$beta==beta){
      return_object <- object_deco
      break
    }
  }
  return_object
}

#' Save Plot as PDF
#'
#' This function saves a given plot object as a PDF file to a specific path.
#'
#' @param plot_object The plot object to be saved.
#' @param file_name String, the name of the file to save the plot as. It should
#'  not include the path or file extension.
#' @param width.parameter Numeric, the width of the output PDF in inches.
#'  Optional, default is NULL.
#' @param height.parameter Numeric, the height of the output PDF in inches.
#'  Optional, default is NULL.
#' @param path String, the directory where the PDF should be saved. Optional,
#'  default is '/mnt_volumen/GIT_REPOSITORIES/MTF/results/'.
#'
#' @details
#' The function first checks whether a path was provided. If not, it uses a
#' default path. Then, it opens a new PDF device, prints the 'plot_object' to
#' this device, and then closes the device, which saves the plot to the
#' specified file. The 'width.parameter' and 'height.parameter' arguments allow
#' for control over the dimensions of the output PDF. If these are not provided,
#' the PDF is created with the default dimensions.
#'
#' @return
#' The function is used for its side effects of creating a PDF file and does not
#' return anything.
#'
#' @import grDevices
#'
#' @export
save_plot_as_pdf <- function(plot_object, file_name, width.parameter = NULL,
                             height.parameter = NULL, path = NULL){

  # Assign the default path if is null.
  if(is.null(path)){
    path <- '/mnt_volumen/GIT_REPOSITORIES/MTF/results/'
  }

  if(is.null(width.parameter) && is.null(height.parameter)){
    grDevices::pdf(file = paste0(path, file_name))
  }else{
    grDevices::pdf(file = paste0(path, file_name),
        width = width.parameter,
        height = height.parameter)
  }

  print(plot_object)
  grDevices::dev.off()
}

#' Plot Performance for Each Model
#'
#' This function plots the performance of each model from a grid search result
#' based on specified alpha and beta values.
#'
#' @param grid_search_result A list object that contains the results of a grid
#'  search.
#' @param alpha Numeric, the alpha value to filter the models to plot. Optional,
#'  default is NULL.
#' @param beta Numeric, the beta value to filter the models to plot. Optional,
#'  default is NULL.
#' @param plot_type Character, specifies the type of plot to produce. Options
#'  are 'summary' (default) or other types supported by the underlying plotting
#'  function.
#' @param scale_parameters Logical, whether to scale the parameters in the plot.
#'  Optional, default is FALSE.
#' @param save.pdf Logical, whether to save the plot as a PDF. Optional, default
#'  is FALSE.
#' @param path.pdf String, the path where to save the PDF file if 'save.pdf' is
#'  TRUE. Optional, default is NULL.
#'
#' @details
#' The function iterates over each model in the grid search result. For each
#' model, it produces a plot based on the 'plot_type' argument.
#' It filters the models to plot based on the 'alpha' and 'beta' arguments.
#'
#' If 'save.pdf' is TRUE, it saves the plot as a PDF file to the specified path
#' or the current working directory if no path is specified.
#'
#' @return
#' The function is used for its side effects of creating and optionally saving
#' plots. It does not return anything.
#'
#'
#' @export
plot_performance_for_each_model <- function(grid_search_result, alpha=NULL,
                                            beta=NULL, plot_type = 'summary',
                                            scale_parameters = FALSE,
                                            save.pdf = FALSE, path.pdf = NULL){
  for(i in 1:length(grid_search_result)){
    object_deco <- grid_search_result[[i]]

    if(plot_type == 'summary'){
      plot_deco <- plot_deconvolution_performance(
        object_deco$running_info,
        title_addition = paste0('alpha=',
                                object_deco$alpha,
                                ', beta=',
                                object_deco$beta))
    }else{
      plot_deco <- plot_deconvolution_performance_x_y_z_total(
        object_deco$running_info,
        title_addition = paste0('alpha=',
                                object_deco$alpha,
                                ', beta=',
                                object_deco$beta),
        alpha = object_deco$alpha,
        beta = object_deco$beta,
        scale_parameters = scale_parameters)
    }

    if(is.null(alpha) && is.null(beta)){
      print(plot_deco)
    }else if(object_deco$alpha==alpha && object_deco$beta==beta){
      print(plot_deco)
    }

    if(save.pdf){
      scaled <- '_scaled'
      if(!scale_parameters)
        scaled <- 'no_scaled'

      save_plot_as_pdf(plot_deco, paste0(plot_type,
                                         scaled,
                                         '_divergence_performance_',
                                         'alpha_', object_deco$alpha, '_beta_',
                                         object_deco$beta, '.pdf'),
                       10, 10, path = path.pdf)
    }

  }
}

#' Plot Detailed Deconvolution Performance
#'
#' This function plots detailed performance of the deconvolution with respect
#' to x, y, z, and total divergence. If specified, the y and z divergence values
#' can be scaled by the given alpha and beta parameters. Also you can see 
#' and additional zoom plot in a specific area.
#'
#' @param running_info Data frame, the running information of the deconvolution.
#' @param title_addition Character, an optional string to be added to the plot
#'  title. Optional, default is ''.
#' @param alpha Numeric, the alpha value to scale the y divergence. Optional,
#'  default is NULL.
#' @param beta Numeric, the beta value to scale the z divergence. Optional,
#'  default is NULL.
#' @param scale_parameters Logical, whether to scale the y and z divergence
#'  values by the given alpha and beta. Optional, default is FALSE.
#' @param include_divergences_x_y_z Logical. Indicates if the x, y and z 
#'  divergences are showed in the plot. Optional, default is TRUE
#' @param include_divergences Logical. Indicates if the divergence is showed 
#'  in the plot. Optional, default is TRUE.
#' @param include_zoom Logical. Indicates if the zoom component is included. 
#'  Optional, default is FALSE.
#' @param zoom_xlim_limits List with both limits to show the zoom: c(100, 1000)
#' @param include_title Logical. Indicates if the title is added to the plot. 
#'  Optional, default is TRUE.
#'
#' @details
#' The function generates a line plot showing the log of the divergence values
#' (for x, y, z, total, and delta) against iterations.
#' The line color is specified according to the type of divergence.
#' The divergence values for y and z can be scaled by the specified alpha and
#' beta parameters if 'scale_parameters' is TRUE. For more information about
#' the zoom functionality: 
#' https://datavizpyr.com/how-to-zoom-in-on-a-plot-in-r/
#'
#' @return
#' A ggplot2 object representing the plot.
#'
#' @import ggplot2 tidyverse ggforce
#' @importFrom ggplot2 .data
#'
#'
#' @export
plot_deconvolution_performance_x_y_z_total <- function(
    running_info,
    title_addition = '',
    alpha=NULL, beta=NULL,
    scale_parameters = FALSE, 
    include_divergences_x_y_z = TRUE,
    include_divergences = TRUE, 
    include_zoom = FALSE,
    zoom_xlim_limits = c(0, 10),
    include_title = TRUE){
  
  running_info.df <- data.frame(running_info)
  
  if(!is.null(alpha) && !is.null(beta) && scale_parameters){
    alpha_scale = alpha
    beta_scale = beta
    title_addition <- paste0('Scaled - ', title_addition)
  }else{
    alpha_scale = 1
    beta_scale = 1
  }
  
  title <- "Iterations vs "
  
  colors_by_lines <- c("delta_divergence_value"="red")
  
  # Solution for the error 'no visible binding for global variable' posted in:
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  # and solved in https://stackoverflow.com/questions/9439256 using:
  # @importFrom ggplot2 .data, and the .data$ prefix.
  
  # Plot with both metric included
  plot <- ggplot2::ggplot(running_info.df, ggplot2::aes(x=.data$iteration)) +
    ggplot2::geom_line(ggplot2::aes(
      y = log(.data$delta_divergence_value),
      color = 'delta_divergence_value'),
      linetype = 'solid', size = 0.3) +
    ggplot2::ylab('log(metric)') +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(size=10))
  
  if(include_divergences_x_y_z){
    plot <- plot + 
      ggplot2::geom_line(ggplot2::aes(
        y = log(.data$divergence_value_x),
        color = 'divergence_value_x'),
        linetype = 'solid', size = 1) +
      ggplot2::geom_line(ggplot2::aes(
        y = log(.data$divergence_value_y * alpha_scale),
        color = 'divergence_value_y'),
        linetype = 'solid', size = 1) +
      ggplot2::geom_line(ggplot2::aes(
        y = log(.data$divergence_value_z * beta_scale),
        color = 'divergence_value_z'),
        linetype = 'solid', size = 1)
    
    title <- paste0(title, 'log(divergence)')
    
    colors_by_lines <- c(colors_by_lines, 
                         "divergence_value_x"="#377EB8",
                         "divergence_value_y"="#4DAF4A",
                         "divergence_value_z"="#984EA3")
  }
  
  if(include_divergences){
    plot <- plot + 
      ggplot2::geom_line(ggplot2::aes(
        y = log(.data$divergence_value),
        color = 'divergence_value'),
        linetype = 'solid',  size = 0.3)
    
    colors_by_lines <- c(colors_by_lines, 
                         "divergence_value"="black")
  }
  
  if(include_zoom){
    plot <- plot + ggforce::facet_zoom(xlim = zoom_xlim_limits)
  }
  
  title <- paste0(title, ' and log(delta divergence)')
  
  # Adding the title
  if(include_title){
    plot <- plot + ggplot2::ggtitle(paste0(title, ' - ', title_addition))
  }
  
  #Adding color lines
  plot <- plot +     
    ggplot2::scale_color_manual(name = "Metric", values=colors_by_lines)

  print(plot)
}

#' Calculate Sparseness of a Matrix
#'
#' This function calculates the sparseness of a given matrix, i.e., the
#' percentage of zero and NA values in the matrix.
#'
#' @param matrix_object Matrix, the matrix for which to calculate sparseness.
#'
#' @details
#' The sparseness of a matrix is defined as the percentage of zero and NA values
#' out of the total number of entries in the matrix. This function counts the
#' number of non-zero and non-NA values in the matrix, and then calculates
#' sparseness as (1 - (non-zero/total)) * 100.
#'
#' @return
#' Numeric, the sparseness of the matrix in percentage.
#'
#' @import Matrix
#'
#'
#' @export
calculate_sparseness_matrix <- function(matrix_object){
  # Entries in the data matrix
  entries.matrix <- nrow(matrix_object) * ncol(matrix_object)

  # Number of nulls, zero and NA values.
  non.zero.na.values.matrix <- Matrix::nnzero(matrix_object, na.counted = FALSE)

  # % sparseness
  sparseness.matrix <- (1-(non.zero.na.values.matrix/entries.matrix)) * 100
  sparseness.matrix
}

#' Order Columns in a Matrix
#'
#' This function orders columns in a matrix based on a given starting position.
#'
#' @param matrix Matrix, the matrix to be ordered.
#' @param start Integer, the column from which to start ordering (default = 7).
#' @param verbose Logical, if TRUE, displays detailed messages during execution
#'  (default = FALSE).
#'
#' @details
#' The function orders the columns in the matrix starting from the position
#' indicated by the 'start' parameter. It creates a new matrix by concatenating
#' the ordered columns and returns it.
#'
#' @return
#' Matrix, the ordered matrix.
#'
#'
#' @export
order_result_generic <- function(matrix, start = 7, verbose = FALSE){

  matrix_return <- NULL
  index <- NULL
  to <- ncol(matrix)

  for(counter in 1:to){
    index <- counter+start-1

    if(verbose){
      print(paste0(counter, ", ", index))
    }

    matrix_return <- cbind(matrix_return, matrix[,paste0('cluster_',index)])
  }
  colnames(matrix_return) <- paste0('cluster_',seq(start,index))

  matrix_return
}

#' Order Elements in a Vector
#'
#' This function orders elements in a vector based on a given starting position.
#'
#' @param vector Vector, the vector to be ordered.
#' @param start Integer, the position from which to start ordering
#'  (default = 5).
#' @param verbose Logical, if TRUE, displays detailed messages during execution
#'  (default = FALSE).
#'
#' @details
#' The function orders the elements in the vector starting from the position
#' indicated by the 'start' parameter. It creates a new vector by concatenating
#' the ordered elements and returns it.
#'
#' @return
#' Vector, the ordered vector.
#'
#'
#' @export
order_simple_vector <- function(vector, start = 5, verbose = FALSE){
  vector_return <- NULL
  index <- NULL
  to <- length(vector)

  for(counter in 1:to){
    index <- counter+start-1

    if(verbose){
      print(paste0(counter, ", ", index))
    }

    vector_return <- cbind(vector_return, vector[paste0('cluster_',index)])
  }
  colnames(vector_return) <- paste0('cluster_',seq(start,index))
  rownames(vector_return) <- 'real'

  vector_return
}

#' Calculate p-values for WT vs KO
#'
#' This function calculates p-values for wild type (WT) vs knockout (KO) based
#' on given results. For more information, check:
#' https://www.cyclismo.org/tutorial/R/pValues.html
#'
#' @param results.prop DataFrame, the results that will be used to perform the
#'  t-tests.
#'
#' @details
#' The function takes a DataFrame of results and transposes it. Then for each
#' cell type present in the data, it extracts the WT and KO results, performs a
#' t-test, and prints the result.
#'
#' @return
#' No return value. The function prints the result of the t-test for each cell
#' type.
#'
#' @import stats
#'
#'
#' @export
calculate_p_values_wt_vs_ko <- function(results.prop){
  results.prop_test <- t(results.prop)
  for (celltype in colnames(results.prop_test)) {
    print(celltype)

    x <- results.prop_test[c('wt_1', 'wt_2', 'wt_3'), celltype]
    y <- results.prop_test[c('ko_1', 'ko_2', 'ko_3'), celltype]
    test.result <- stats::t.test(as.numeric(x),as.numeric(y))
    print(test.result)
  }
}
