#' Tree Guided Deconvolution Functions
#'
#' This script, tree_guided_deconvolution.R, is part of the Decoflex package.
#' It includes the following functions:
#'
#' \describe{
#'   \item{add_grouping_single_cell()}{A function for adding single cell
#'   groupings.}
#'   \item{aggregate_true_top_clustering_prop()}{A function to aggregate true
#'   top clustering proportions.}
#'   \item{calculate_joint_cor()}{A function to calculate joint correlation.}
#'   \item{calculate_markers()}{A function to calculate markers.}
#'   \item{calculate_markers_level_intersection()}{A function to calculate
#'   markers level intersection.}
#'   \item{calculate_min_correlation_incremental_markers()}{A function to
#'   calculate minimum correlation incremental markers.}
#'   \item{calculate_performance_metrics()}{A function to calculate performance
#'   metrics.}
#'   \item{calculate_prop_back_propagation()}{A function to calculate proportion
#'   back propagation.}
#'   \item{FindMarkers.proxy()}{A proxy function to find markers.}
#'   \item{join_back_propagation_proportions_leaves()}{A function to join back
#'   propagation proportions leaves.}
#'   \item{performance_metrics_general()}{A function to calculate general
#'   performance metrics.}
#'   \item{print_performance_metrics()}{A function to print performance
#'   metrics.}
#'   \item{recalculate_proportions()}{A function to recalculate proportions.}
#'   \item{run_deconvolution_simulation_generic()}{A function to run generic
#'   deconvolution simulations.}
#'   \item{run_deconvolution_tree_guided()}{A function to run tree-guided
#'   deconvolution.}
#'   \item{run_deconvolution_tree_guided_recursive()}{A recursive function to
#'   run tree-guided deconvolution.}
#' }
#'


#' @title Deconvolution Simulation
#'
#' @description This function performs a detailed simulation of the
#' deconvolution process on single cell expression data. It allows for the
#' inclusion of multiple cluster levels (both top level and sub level clusters),
#' and customization of various parameters such as marker strategy, log fold
#' change threshold, and adjusted p-value threshold. The function also supports
#' a minimum correlation strategy for deconvolution, providing flexibility in
#' the simulation approach. Verbose mode can be enabled for additional execution
#' details.
#' Note: If you encounter an error related to the function 'checkHT', you can
#' load it from: https://svn.r-project.org/R/trunk/src/library/utils/R/head.R
#'
#' @name run_deconvolution_simulation
#'
#' @param single_cell_data_exp A data frame representing the single cell data
#'  for the deconvolution simulation.
#' @param top_clusters_var A vector indicating the top clusters variables.
#' @param top_clusters_list A list containing the top clusters.
#' @param sub_clusters_var A vector indicating the sub clusters variables.
#' @param sub_clusters_list A list containing the sub clusters.
#' @param sample A numeric or integer value indicating the sample size for the
#'  deconvolution simulation.
#' @param use_min_cor_strategy A boolean. If TRUE, the function uses minimum
#'  correlation strategy. Default is FALSE.
#' @param delete_shared_level_markers A boolean. If TRUE, shared level markers
#'  will be deleted. Default is TRUE.
#' @param delete_shared_internal_markers A boolean. If TRUE, shared internal
#'  markers will be deleted. Default is TRUE.
#' @param filter_markers A vector of markers to filter, or NULL. Default is NULL
#' @param param.logfc.threshold A numeric value indicating the log fold change
#'  threshold. Default is 2.0.
#' @param param.p_val_adj A numeric value indicating the adjusted p-value
#'  threshold. Default is 0.05.
#' @param marker_strategy A character string specifying the marker strategy.
#'  Default is 'keep_default_values'.
#' @param verbose A boolean. If TRUE, the function will print additional details
#'  during the execution. Default is FALSE.
#'
#' @return A list containing the results of the deconvolution simulation. This
#'  might include the estimated sources, the final weights, and any other
#'  relevant output depending on the actual implementation of the function.
#'
#' @export
run_deconvolution_simulation <- function(single_cell_data_exp,
                                         top_clusters_var, top_clusters_list,
                                         sub_clusters_var, sub_clusters_list,
                                         sample,
                                         use_min_cor_strategy = FALSE,
                                         delete_shared_level_markers = TRUE,
                                         delete_shared_internal_markers = TRUE,
                                         filter_markers = NULL,
                                         param.logfc.threshold = 2.0,
                                         param.p_val_adj = 0.05,
                                         marker_strategy ='keep_default_values',
                                         verbose = FALSE){

  # In case I want I clean running.
  if(verbose){
    result_deco <- run_deconvolution_simulation_generic(
      single_cell_data_exp = single_cell_data_exp,
      top_clusters_var = top_clusters_var,
      top_clusters_list = top_clusters_list,
      sub_clusters_var = sub_clusters_var,
      sub_clusters_list = sub_clusters_list,
      sample = sample,
      use_min_cor_strategy = use_min_cor_strategy,
      delete_shared_level_markers = delete_shared_level_markers,
      delete_shared_internal_markers = delete_shared_internal_markers,
      filter_markers = filter_markers,
      param.logfc.threshold = param.logfc.threshold,
      param.p_val_adj = param.p_val_adj,
      marker_strategy = marker_strategy,
      verbose = verbose)

  }else{
    result_deco <- (run_deconvolution_simulation_generic(
      single_cell_data_exp = single_cell_data_exp,
      top_clusters_var = top_clusters_var,
      top_clusters_list = top_clusters_list,
      sub_clusters_var = sub_clusters_var,
      sub_clusters_list = sub_clusters_list,
      sample = sample,
      use_min_cor_strategy = use_min_cor_strategy,
      delete_shared_level_markers = delete_shared_level_markers,
      delete_shared_internal_markers = delete_shared_internal_markers,
      filter_markers = filter_markers,
      param.logfc.threshold = param.logfc.threshold,
      param.p_val_adj = param.p_val_adj,
      marker_strategy = marker_strategy,
      verbose = verbose))
  }

  return(result_deco)
}

#' @title Check Acceptability of 'n' for 'head()' and 'tail()' Methods
#'
#' @description This function checks the acceptability of 'n' for 'head()' and
#' 'tail()' methods. If you encounter an error related to this function, you
#' can load it from:
#' https://svn.r-project.org/R/trunk/src/library/utils/R/head.R
#'
#' @name checkHT
#'
#' @param n The number to be checked for its acceptability.
#' @param d Data input for the function.
#'
#' @return Returns TRUE if 'n' is acceptable, and FALSE otherwise.
#'
#' @export
checkHT <- function(n, d) {
  len <- length(n)
  msg <- if(len == 0 || all(is.na(n)))
    gettext(
      "invalid 'n' - must contain at least one non-missing element, got none.")
  else if(!(is.numeric(n) || is.logical(n)))
    gettext("invalid 'n' - must be numeric, possibly NA.")
  else if(is.null(d) && len > 1L)
    gettextf("invalid 'n' - must have length one when dim(x) is NULL, got %d",
             len)
  else if(!is.null(d) && len > length(d))
    gettextf("invalid 'n' - length(n) must be <= length(dim(x)), got %d > %d",
             len,
             length(d))
  else return(invisible())
  stop(msg, domain = NA)
}

#' @title Generic Deconvolution Simulation
#'
#' @description This function is a generalized form of the deconvolution
#' simulation process applicable to single cell data. It allows users to modify
#' various parameters like the marker strategy, log fold change threshold, and
#' adjusted p-value threshold. By setting use_min_cor_strategy to TRUE, the
#' function will apply a minimum correlation strategy during the deconvolution.
#'
#' @name run_deconvolution_simulation_generic
#'
#' @param single_cell_data_exp A data frame representing the single cell data
#'  for the deconvolution simulation.
#' @param top_clusters_var A vector indicating the top clusters variables.
#' @param top_clusters_list A list containing the top clusters.
#' @param sub_clusters_var A vector indicating the sub clusters variables.
#' @param sub_clusters_list A list containing the sub clusters.
#' @param sample A numeric or integer value indicating the sample size for the
#'  deconvolution simulation.
#' @param use_min_cor_strategy A boolean. If TRUE, the function uses minimum
#'  correlation strategy. Default is FALSE.
#' @param delete_shared_level_markers A boolean. If TRUE, shared level markers
#'  will be deleted. Default is TRUE.
#' @param delete_shared_internal_markers A boolean. If TRUE, shared internal
#'  markers will be deleted. Default is TRUE.
#' @param filter_markers A vector of markers to filter, or NULL. Default is NULL
#' @param param.logfc.threshold A numeric value indicating the log fold change
#'  threshold. Default is 2.0.
#' @param param.p_val_adj A numeric value indicating the adjusted p-value
#'  threshold. Default is 0.05.
#' @param marker_strategy A character string specifying the marker strategy, or
#'  NULL. Default is NULL.
#' @param verbose A boolean. If TRUE, the function will print additional details
#'  during the execution. Default is FALSE.
#'
#' @return A list containing the results of the deconvolution simulation. This
#'  might include the estimated sources,
#'  the final weights, and any other relevant output depending on the actual
#'  implementation of the function.
#'
#'
#' @export
run_deconvolution_simulation_generic <- function(
    single_cell_data_exp,
    top_clusters_var,
    top_clusters_list,
    sub_clusters_var,
    sub_clusters_list,
    sample,
    use_min_cor_strategy = FALSE,
    delete_shared_level_markers = TRUE,
    delete_shared_internal_markers = TRUE,
    filter_markers = NULL,
    param.logfc.threshold = 2.0,
    param.p_val_adj = 0.05,
    marker_strategy = 'keep_default_values',
    verbose = FALSE){

  # 1. calculate the simulated bulk data.
  pseudo.bulk.data <- generateBulk_allcells(single_cell_data_exp,
                                            ct.varname = sub_clusters_var,
                                            sample = "sample",
                                            ct.sub = NULL,
                                            verbose = verbose)

  # 2. Call the more generic function
  result_deco <- run_deconvolution_tree_guided(
    bulk_data = data.frame(pseudo.bulk.data$pseudo_eset@assayData$exprs),
    single_cell_data_exp = single_cell_data_exp,
    top_clusters_var = top_clusters_var,
    top_clusters_list = top_clusters_list,
    sub_clusters_var = sub_clusters_var,
    sub_clusters_list = sub_clusters_list,
    sample = sample,
    use_min_cor_strategy = use_min_cor_strategy,
    true_proportions = as.data.frame.matrix(pseudo.bulk.data$truep),
    delete_shared_level_markers = delete_shared_level_markers,
    delete_shared_internal_markers = delete_shared_internal_markers,
    filter_markers = filter_markers,
    param.logfc.threshold = param.logfc.threshold,
    param.p_val_adj = param.p_val_adj,
    marker_strategy = marker_strategy,
    verbose = verbose)

  # 3. Return the results.
  return(list(result_deco = result_deco, pseudo.bulk.data = pseudo.bulk.data))
}




#' run_deconvolution_simulation_generic_recursive
#'
#' @description Perform a recursive deconvolution simulation on given single
#' cell data. This function handles the first top node of the hierarchy and then
#' send the rest to the recursive algorithm
#'
#' @param single_cell_data_exp A data frame that contains single cell expression
#'  data.
#' @param hierarchy A list that outlines the hierarchy of cell types.
#' @param sub_clusters_var A variable that represents the subclusters in the
#'  data.
#' @param sample A variable that represents the samples in the data.
#' @param use_min_cor_strategy A logical parameter indicating whether to use
#'  minimum correlation strategy. Defaults to TRUE.
#' @param delete_shared_level_markers A logical parameter indicating whether to
#'  delete shared level markers. Defaults to FALSE.
#' @param delete_shared_internal_markers A logical parameter indicating whether
#'  to delete shared internal markers. Defaults to FALSE.
#' @param deconvolute_just_top A logical parameter indicating whether to only
#'  deconvolute top. Defaults to FALSE.
#' @param filter_markers A parameter for specifying markers to filter, defaults
#'  to NULL.
#' @param param.logfc.threshold A numeric value specifying the threshold for
#'  log fold change, defaults to 2.0.
#' @param param.p_val_adj A numeric value specifying the threshold for adjusted
#'  p-values, defaults to 0.05.
#' @param marker_strategy Strategy for selecting markers, defaults to NULL.
#' @param verbose A logical parameter indicating whether to print detailed
#'  messages during the execution of the function. Defaults to FALSE.
#'
#' @return A list containing the deconvolution result and the pseudo bulk data.
#'
#' @export
run_deconvolution_simulation_generic_recursive <- function(
    single_cell_data_exp,
    hierarchy,
    sub_clusters_var,
    sample,
    use_min_cor_strategy = TRUE,
    delete_shared_level_markers = FALSE,
    delete_shared_internal_markers = FALSE,
    deconvolute_just_top = FALSE,
    filter_markers = NULL,
    param.logfc.threshold = 2.0,
    param.p_val_adj = 0.05,
    marker_strategy = 'keep_default_values',
    verbose = FALSE){

  # TODO: 0. Check the hierarchy object

  # 1. calculate the simulated bulk data. I CAN'T add the root cell type list
  # meaning (ct.sub = hierarchy$`1`$celltype_list) because we need a simulation
  # of all cell-types; this is due to the fact that the method should recover
  # the proportions of a mixture of cell-types, event the unknow ones. Given
  # this, ct.sub must be equal to NULL.
  pseudo.bulk.data <- generateBulk_allcells(single_cell_data_exp,
                                            ct.varname = sub_clusters_var,
                                            sample = sample,
                                            ct.sub = NULL,
                                            verbose = verbose)

  # First let's asign the column names (samples) again, just in case
  # those start with a number; in that case R will add an X and everything
  # will be wrong.
  bulk_data.df <- data.frame(pseudo.bulk.data$pseudo_eset@assayData$exprs)
  colnames(bulk_data.df) <- colnames(pseudo.bulk.data$pseudo_eset@assayData$exprs)

  # 2.Filter the original single cell data to have just the celltypes that are
  # needed for the process.
  single_cell_data_exp <-
    single_cell_data_exp[,
                         single_cell_data_exp[[sub_clusters_var]] %in%
                           hierarchy$`1`$celltype_list]

  # 2.1. In case that I have just one sample, I need to be sure that the data
  # frame is built.
  true_proportions_parameter <- NULL
  if(is.null(nrow(pseudo.bulk.data$truep))){
    true_proportions_parameter <- data.frame(t(pseudo.bulk.data$truep))
    rownames(true_proportions_parameter) <-
      colnames(pseudo.bulk.data$pseudo_eset)
  }else{
    true_proportions_parameter <- as.data.frame.matrix(pseudo.bulk.data$truep)
  }

  # 3. Call the more generic function. I send the next level of clustering
  # object to initiate the recursive process. :).
  # Also the top_proportion is null since it is the first level of the tree.
  result_deco <- run_deconvolution_tree_guided_recursive(
    result_deco_top = NULL,
    bulk_data = bulk_data.df,
    true_proportions = true_proportions_parameter[,hierarchy$`1`$celltype_list],
    single_cell_data_exp = single_cell_data_exp,
    sub_clusters_var = sub_clusters_var,
    hierarchy = hierarchy$`1`,
    sample = sample,
    use_min_cor_strategy = use_min_cor_strategy,
    delete_shared_level_markers = delete_shared_level_markers,
    delete_shared_internal_markers = delete_shared_internal_markers,
    deconvolute_just_top = deconvolute_just_top,
    param.logfc.threshold = param.logfc.threshold,
    param.p_val_adj = param.p_val_adj,
    marker_strategy = marker_strategy,
    verbose = verbose)

  # 4. Return the results.
  return(list(result_deco = result_deco, pseudo.bulk.data = pseudo.bulk.data))
}


#' Run Deconvolution Tree Guided Recursively
#'
#' This function uses a tree-guided method to perform a recursive deconvolution
#' operation on bulk and single-cell data. The recursive method helps in
#' analyzing hierarchically structured data, usually encountered in biological
#' systems where a set of cell types are hierarchically organized based on
#' their molecular similarity. The method starts at the top level,
#' identifies clusters, calculates markers, performs deconvolution, and
#' propagates the information back to guide the deconvolution at the next level.
#'
#' @param result_deco_top A list representing the result of deconvolution at
#'  the top level of the hierarchy. If NULL (default), the function
#'  starts at the root of the tree.
#' @param bulk_data The bulk data on which the deconvolution needs to be
#'  performed.
#' @param true_proportions An optional matrix representing the true proportions
#'  of cell types. This can be used to calculate the accuracy of
#'  the deconvolution. If NULL (default), the accuracy calculation is skipped.
#' @param single_cell_data_exp A data frame representing single cell expression
#'  data, used to guide the deconvolution process.
#' @param sub_clusters_var A string representing the variable in
#'  single_cell_data_exp which indicates sub-clusters.
#' @param hierarchy A list structure describing the hierarchy of cell type
#'  clustering.
#' @param sample A string representing the sample name.
#' @param use_min_cor_strategy A boolean indicating whether to use the minimum
#'  correlation strategy during the marker calculation phase. Default is TRUE.
#' @param delete_shared_level_markers A boolean indicating whether to delete
#'  shared level markers during the marker calculation phase.
#' Default is FALSE.
#' @param delete_shared_internal_markers A boolean indicating whether to delete
#'  shared internal markers during the marker calculation phase.
#'  Default is FALSE.
#' @param deconvolute_just_top A boolean indicating whether to perform
#'  deconvolution only at the top level. If TRUE (default is FALSE), the
#'  function stops after deconvoluting the top level and returns the result.
#' @param deconvolute_top_hierarchy_limit An integer indicating the top
#'  hierarchy limit for deconvolution. Default is 1.
#' @param filter_markers An optional list of markers to be used for filtering
#'  during the marker calculation phase. If NULL (default), no filtering is
#'  performed.
#' @param param.logfc.threshold A numeric value representing the log-fold-change
#'  threshold for marker calculation. Default is 2.0.
#' @param param.p_val_adj A numeric value representing the adjusted p-value
#'  threshold for marker calculation. Default is 0.05.
#' @param test.use.value A string indicating the statistical test to use for
#'  marker calculation. Default is 'wilcox'.
#' @param marker_strategy An optional string representing the marker selection
#'  strategy. If NULL (default), a default strategy is used.
#' @param minimum_markers The minimum number of marker genes to be selected,
#'  default is set to 4. This can be adjusted based on the analysis
#'  requirements.
#' @param min_delta_cor_threshold Minimum difference between the correlation of
#'  t and t-1 steps. The idea is to skip genes that have small changes in
#'  the correlation. For default it is a 5%.
#' @param verbose A boolean indicating whether to print detailed messages during
#'  the execution of the function. Default is FALSE.
#'
#' @return A list containing the result of the top level deconvolution, a list
#'  of deconvolution results at each level, and the back-propagated proportions
#'  at the top level.
#'
#' @import stats
#'
#' @export
run_deconvolution_tree_guided_recursive <- function(
    result_deco_top =  NULL,
    bulk_data,
    true_proportions = NULL,
    single_cell_data_exp,
    sub_clusters_var,
    hierarchy,
    sample,
    use_min_cor_strategy = TRUE,
    delete_shared_level_markers = FALSE,
    delete_shared_internal_markers = FALSE,
    deconvolute_just_top = FALSE,
    deconvolute_top_hierarchy_limit = 1,
    filter_markers = NULL,
    param.logfc.threshold = 2.0,
    param.p_val_adj = 0.05,
    test.use.value = 'wilcox',
    marker_strategy = 'keep_default_values',
    ordering_strategy = 'pvalue_foldchange',
    minimum_markers = 4,
    percentile_markers = NULL,
    min_delta_cor_threshold = 0.05,
    percentile_markers.min_corr = NULL,
    max_iterations = 10000,
    delta_threshold = 1e-10,
    markers.clusters.parameter = NULL,
    verbose = FALSE){

  message(paste0('run_deconvolution_tree_guided_recursive - Verbose: ',
                 verbose))

  # Let's begin with a basic validation about the samples names staring with #
  current_columns <- unique(single_cell_data_exp[[sample]])
  column_validation <- all(
    unlist(lapply(current_columns,
                  function(text) !grepl("^[[:digit:]]+", text))))
  if(!column_validation){
    stop(paste0('There are samples in the single cell data that start with a',
    ' number. Please make sure that all samples start with a letter. Some R ',
    'functionalities do not work with this kind of naming conventions.'))
  }

  # Creation of the level name that I will use for the top deconvolution
  top_clusters_var <- paste0('metacluster_level_', hierarchy$tree_level)

  if(verbose){
    message('Recursive version tree-guided method: level ',
            hierarchy$tree_level, ', leaf ', hierarchy$leaf_number,
            ', strategy: ', use_min_cor_strategy)
  }

  #-1. Checking the parameters in the object to locally use it
  param.logfc.threshold.use <- param.logfc.threshold
  if(!is.null(hierarchy$parameters$param.logfc.threshold)){
    param.logfc.threshold.use <- hierarchy$parameters$param.logfc.threshold
    message('Using param.logfc.threshold from object hierarchy: ',
            param.logfc.threshold)
  }

  param.p_val_adj.use <- param.p_val_adj
  if(!is.null(hierarchy$parameters$param.p_val_adj)){
    param.p_val_adj.use <- hierarchy$parameters$param.p_val_adj
    message('Using param.p_val_adj from object hierarchy: ',
            param.p_val_adj.use)
  }

  minimum_markers.use <- minimum_markers
  if(!is.null(hierarchy$parameters$minimum_markers)){
    minimum_markers.use <- hierarchy$parameters$minimum_markers
    message('Using minimum_markers from object hierarchy: ',
            minimum_markers)
  }

  min_delta_cor_threshold.use <- min_delta_cor_threshold
  if(!is.null(hierarchy$parameters$min_delta_cor_threshold)){
    min_delta_cor_threshold.use <- hierarchy$parameters$min_delta_cor_threshold
    message('Using min_delta_cor_threshold from object hierarchy: ',
            min_delta_cor_threshold)
  }

  markers.clusters.parameter.use <- markers.clusters.parameter
  if(!is.null(hierarchy$parameters$markers.clusters.parameter)){
    markers.clusters.parameter.use <- hierarchy$parameters$markers.clusters.parameter
    message('Using markers.clusters.parameter from object hierarchy: ',
            paste0(length(markers.clusters.parameter.use$total_markers), ' markers.'))
  }


  # 0. Create the top clustering on the dataset in a given variable or in the
  # final configuration if next_level_clustering is NULL
  if(!is.null(hierarchy$next_level_clustering)){
    single_cell_data_exp <- add_grouping_single_cell(
      single_cell_data_exp = single_cell_data_exp,
      sub_clusters_var = sub_clusters_var,
      level_group_name = top_clusters_var,
      hierarchy_level_groups = hierarchy$next_level_clustering)

    # I got the list of top clusters
    top_clusters_list <-
      unique(stats::na.omit(single_cell_data_exp[[top_clusters_var]]))
  }else{
    # Because this is the end of the tree, I will use the original cell type
    # variable and the original names as well
    top_clusters_list <- unique(stats::na.omit(hierarchy$celltype_list))
    top_clusters_var <- sub_clusters_var
  }

  # 1. Create the reference for the top clustering configuration
  reference_w.top.cluster <- decoflex_build_cell_reference(
    x = single_cell_data_exp,
    ct.sub = top_clusters_list,
    ct.varname = top_clusters_var,
    sample = sample,
    verbose = verbose)

  # 2. Calculate the Marker genes for the top clusters
  # In this case, I can pass the object as parameter
  markers.top.clusters.object <- NULL
  if(is.null(markers.clusters.parameter.use)){
    markers.top.clusters.object <- calculate_markers(
      single_cell_data_exp = single_cell_data_exp,
      reference = reference_w.top.cluster$basis,
      group_clusters_var = top_clusters_var,
      use_min_cor_strategy = use_min_cor_strategy,
      delete_shared_internal_markers = delete_shared_internal_markers,
      filter_markers = filter_markers,
      param.logfc.threshold = param.logfc.threshold.use,
      param.p_val_adj = param.p_val_adj.use,
      test.use.value = test.use.value,
      marker_strategy = marker_strategy,
      ordering_strategy = ordering_strategy,
      minimum_markers = minimum_markers.use,
      percentile_markers = percentile_markers,
      min_delta_cor_threshold = min_delta_cor_threshold.use,
      percentile_markers.min_corr = percentile_markers.min_corr,
      verbose = verbose)
  }else{
    markers.top.clusters.object <- markers.clusters.parameter.use
  }

  # 2.1. Assign the final markers and the top marker list
  markers.top.clusters <- markers.top.clusters.object$total_markers
  list.top.marker <- markers.top.clusters.object$list_markers

  # 3. Run de deconvolution of the top clusters using the markers genes for the
  # top clusters.
  result_deco_top_cluster <- run_standard_deconvolution(
    bulk_data_x = bulk_data,
    references_w = reference_w.top.cluster$basis,
    markers = markers.top.clusters,
    max_iterations = max_iterations,
    delta_threshold = delta_threshold,
    verbose = verbose)

  # 3. Run the proportion back propagation with respect to the top proportions
  # and the current level-leaf If this is Null, is because the current process
  # is the root, meaning the first clustering.
  back_propagation_proportions_top = NULL
  if(!is.null(result_deco_top)){
    #3.1 calculate the proportions for the top model
    back_propagation_proportions_top <- calculate_prop_back_propagation(
      top_proportions = result_deco_top,
      level_proportions = result_deco_top_cluster$h,
      verbose = verbose)
  }else{

    if(verbose){
      message('Current process corresponding to the root tree. Not back ',
              'propagation activated. Therefore the workflow will use the top ',
              'proportions.')
    }

    back_propagation_proportions_top <- result_deco_top_cluster$h
  }

  # 4.If the current process is a simulation, I calculate the associated
  # metrics.
  if(!is.null(true_proportions)){
    # 3.a. Calculate the relative proportions with respect to the subset of
    # subcluster taken.
    true_prop_relative <- true_proportions[, hierarchy$celltype_list]

    # rescaling of the proportions relative JUST TO THE SUBCLUSTERS
    true_prop_relative <- true_prop_relative/rowSums(true_prop_relative)

    # 3.b calculate the proportions for the top model
    true_prop_top_clustering <- aggregate_true_top_clustering_prop(
      single_cell_data_exp = single_cell_data_exp,
      true_prop = true_prop_relative,
      top_clusters_var = top_clusters_var,
      top_clusters_list = top_clusters_list,
      sub_clusters_var = sub_clusters_var)

    # 4. Now I can calculate the metrics for both configurations: samples and
    # celltypes
    pm_top_clustering <- performance_metrics_general(
      true_prop = true_prop_top_clustering,
      calc_prop = result_deco_top_cluster$h,
      title = paste0('Top clustering: ', top_clusters_var),
      verbose = verbose)

    # adding the metric information to the result of the top model
    result_deco_top_cluster[['simulation_metrics']] <- pm_top_clustering
  }


  # 5. Check the marker genes for the level to have the shared level markers
  # If the parameter is false, I will avoid this calculation.
  # In addition, if I'm analizing a cluster with just one cell type, I will
  # avoid this calculation since it doesn't make sense to
  # calculate it.
  shared_level_markers = NULL
  if(delete_shared_level_markers){
    shared_level_markers <- calculate_markers_level_intersection(
      single_cell_data_exp = single_cell_data_exp,
      top_clusters_var = top_clusters_var,
      sub_clusters_var = sub_clusters_var,
      top_clusters_list = top_clusters_list,
      param.logfc.threshold = param.logfc.threshold.use,
      param.p_val_adj = param.p_val_adj.use,
      minimum_markers = minimum_markers.use,
      percentile_markers = percentile_markers,
      min_delta_cor_threshold = min_delta_cor_threshold.use,
      percentile_markers.min_corr = percentile_markers.min_corr,
      verbose = verbose)
  }

  # Create a variable to save the results of each deconvolution
  deco_results_list <- list()

  # If we want to deconvolute just the top because that is the initial result
  # that is important in a specific experiment.
  if(deconvolute_just_top){
    hierarchy$next_level_clustering = NULL

    if(verbose){
      message('The deconvolution has been executed just for the TOP LEVEL',
              'following the parameter: deconvolute_just_top = TRUE. ')
    }
  }else{
    if(verbose){
      message('The deconvolution will execute down this level: ',
              'deconvolute_just_top = FALSE ')
    }
  }

  # 6. For each top cluster create the deconvolution for the subclusters
  # inside. if the subclustering is the final level, this is going to
  # stop the entire recursive process.
  activate_back_propagation <- TRUE
  for(current_hierarchy in hierarchy$next_level_clustering) {

    top_cluster <- paste0('subcluster_level_', current_hierarchy$tree_level,
                          '_leaf_',  current_hierarchy$leaf_number)

    deconvolute_current_level <- current_hierarchy[['deconvolute']]

    # In case that I want to keep deconvoluting until the last leaf.
    if(is.null(deconvolute_current_level)){
      deconvolute_current_level <- TRUE
    }else{

      if(deconvolute_current_level){
        # Normal deconvolution
      }else{
        activate_back_propagation <- FALSE

        # nothing, now I can use the real value to keep deconvoluting or
        # stopping.
        message('The group is not going to be deconvoluted: ',
                'deconvolute_current_level=FALSE -> ', top_cluster)

        # continue with the loop
        next
      }
    }

    # 6.1. check the subcluster that are part of the top cluster analyzed.
    values_current_top_cluster <-
      single_cell_data_exp[[top_clusters_var]] == top_cluster
    values_subclusters <- single_cell_data_exp[[sub_clusters_var]]
    subclusters_list <- unique(values_subclusters[values_current_top_cluster])
    # This is the final list of celltypes that are related with this leaf
    subclusters_list <-
      levels(factor(subclusters_list[!is.na(subclusters_list)]))

    # 6.2. This is just a double check since both list must coincide.
    if(length(intersect(subclusters_list, current_hierarchy$celltype_list)) !=
       length(current_hierarchy$celltype_list)){
      stop(paste0('The hierarchy definition does not coincide with the single cell ',
                  'data. Check input parameters.'))
    }

    # 6.3. There are cases where the leaf is just one cell type, therefore, the
    # deconvolution is not needed and the proportions calculated on the top
    # level correspond to the leaf of one cell type. Low rank >1 [2,3,4]
    result_deco_subclusters = NULL
    if(length(subclusters_list) > 1){
      result_deco_subclusters <- run_deconvolution_tree_guided_recursive(
        #' filter the top proportions with just the one for the leaf
        result_deco_top = back_propagation_proportions_top[top_cluster,],
        bulk_data = bulk_data,
        single_cell_data_exp = single_cell_data_exp,
        sub_clusters_var = sub_clusters_var,
        hierarchy = current_hierarchy,
        sample = sample,
        use_min_cor_strategy = use_min_cor_strategy,
        markers.clusters.parameter = markers.clusters.parameter.use,
        true_proportions = true_proportions[, current_hierarchy$celltype_list],
        delete_shared_level_markers = delete_shared_level_markers,
        delete_shared_internal_markers = delete_shared_internal_markers,
        filter_markers = filter_markers,
        param.logfc.threshold = param.logfc.threshold.use,
        param.p_val_adj = param.p_val_adj.use,
        test.use.value = test.use.value,
        marker_strategy = marker_strategy,
        ordering_strategy = ordering_strategy,
        minimum_markers = minimum_markers.use,
        percentile_markers = percentile_markers,
        min_delta_cor_threshold = min_delta_cor_threshold.use,
        max_iterations = max_iterations,
        delta_threshold = delta_threshold,
        percentile_markers.min_corr = percentile_markers.min_corr,
        verbose = verbose)

    }else if(length(subclusters_list) == 1){

      if(verbose){
        message('Deconvolution with just one celltype (',
                top_cluster, '): ', subclusters_list)
      }

      # I create manually the object expected for a normal deconvolution
      # getting the same proportions of the top clustering rescaled but
      # just for the current clustering.
      proportions_one_celltype <- back_propagation_proportions_top[top_cluster,]

      # Since the top proportion has a length of 1, it also has the name of the
      # top cluster, therefore I need to rename it with a more deep cluster,
      # meaning the current one.
      rownames(proportions_one_celltype) <- subclusters_list

      # Since the name is corresponding to a actual cluster I change the name
      # for it but I keep the standard one to know that it is only one in the
      # calculation.
      markers.top.clusters.object$list_markers[[top_cluster]]$group_name <-
        paste0('markers_', subclusters_list)
      result_deco_subclusters <- list(
        result_deco_top_cluster = NULL,
        deco_results_list = NULL,
        back_propagation_proportions_top = proportions_one_celltype,
        back_propagation_proportions_top_detailed = proportions_one_celltype,
        markers.top.clusters.object = markers.top.clusters.object)
    }

    # 6.3.a. Add the results to the global results.
    deco_results_list[[top_cluster]] <- result_deco_subclusters

  }

  # 7. Joining the back_propagation proportions of each model if all the
  # submodels are calculated
  back_propagation_proportions_top_detailed = NULL
  if(!is.null(hierarchy$next_level_clustering) & activate_back_propagation){

    # The difference between back_propagation_proportions_top and
    # back_propagation_proportions_top_detailed is that the first one is the
    # conciliation of the current top deconvolution with the relative top
    # (before percentage of the leaf), but the second is the concatenation of
    # the proportions that have been already rescaled of the leaves
    # deconvolutions
    back_propagation_proportions_top_detailed <-
      join_back_propagation_proportions_leaves(
        deco_results_list = deco_results_list,
        verbose = verbose)

  }else{
    message('Current process corresponding to a final leaf. Not back ',
            'propagation detailed activated, therefore I will put the top one.')
    back_propagation_proportions_top_detailed = back_propagation_proportions_top
  }

  # 8. Return the final variables.
  return(list(
    result_deco_top_cluster = result_deco_top_cluster,
    deco_results_list = deco_results_list,
    back_propagation_proportions_top = back_propagation_proportions_top,
    back_propagation_proportions_top_detailed =
      back_propagation_proportions_top_detailed,
    markers.top.clusters.object = markers.top.clusters.object))
}

#' Run marker selection OMiC
#'
#' This function return multiple objects with the OMiC marker selection method,
#' result, such you can visualize how the marker genes for each celltype has
#' been selected. Its the same function implemented in the initial part of the
#' function: run_deconvolution_tree_guided_recursive.
#'
#' @param result_deco_top A list representing the result of deconvolution at
#'  the top level of the hierarchy. If NULL (default), the function
#'  starts at the root of the tree.
#' @param bulk_data The bulk data on which the deconvolution needs to be
#'  performed.
#' @param true_proportions An optional matrix representing the true proportions
#'  of cell types. This can be used to calculate the accuracy of
#'  the deconvolution. If NULL (default), the accuracy calculation is skipped.
#' @param single_cell_data_exp A data frame representing single cell expression
#'  data, used to guide the deconvolution process.
#' @param sub_clusters_var A string representing the variable in
#'  single_cell_data_exp which indicates sub-clusters.
#' @param hierarchy A list structure describing the hierarchy of cell type
#'  clustering.
#' @param sample A string representing the sample name.
#' @param use_min_cor_strategy A boolean indicating whether to use the minimum
#'  correlation strategy during the marker calculation phase. Default is TRUE.
#' @param delete_shared_level_markers A boolean indicating whether to delete
#'  shared level markers during the marker calculation phase.
#' Default is FALSE.
#' @param delete_shared_internal_markers A boolean indicating whether to delete
#'  shared internal markers during the marker calculation phase.
#'  Default is FALSE.
#' @param deconvolute_just_top A boolean indicating whether to perform
#'  deconvolution only at the top level. If TRUE (default is FALSE), the
#'  function stops after deconvoluting the top level and returns the result.
#' @param deconvolute_top_hierarchy_limit An integer indicating the top
#'  hierarchy limit for deconvolution. Default is 1.
#' @param filter_markers An optional list of markers to be used for filtering
#'  during the marker calculation phase. If NULL (default), no filtering is
#'  performed.
#' @param param.logfc.threshold A numeric value representing the log-fold-change
#'  threshold for marker calculation. Default is 2.0.
#' @param param.p_val_adj A numeric value representing the adjusted p-value
#'  threshold for marker calculation. Default is 0.05.
#' @param test.use.value A string indicating the statistical test to use for
#'  marker calculation. Default is 'wilcox'.
#' @param marker_strategy An optional string representing the marker selection
#'  strategy. If NULL (default), a default strategy is used.
#' @param minimum_markers The minimum number of marker genes to be selected,
#'  default is set to 4. This can be adjusted based on the analysis
#'  requirements.
#' @param min_delta_cor_threshold Minimum difference between the correlation of
#'  t and t-1 steps. The idea is to skip genes that have small changes in
#'  the correlation. For default it is a 5%.
#' @param verbose A boolean indicating whether to print detailed messages during
#'  the execution of the function. Default is FALSE.
#'
#' @return Object with the marker gene structure.
#'
#' @import stats
#'
#' @export
run_marker_selection_OMiC <- function(
    result_deco_top =  NULL,
    bulk_data,
    true_proportions = NULL,
    single_cell_data_exp,
    sub_clusters_var,
    hierarchy,
    sample,
    use_min_cor_strategy = TRUE,
    delete_shared_level_markers = FALSE,
    delete_shared_internal_markers = FALSE,
    filter_markers = NULL,
    param.logfc.threshold = 2.0,
    param.p_val_adj = 0.05,
    test.use.value = 'wilcox',
    marker_strategy = 'keep_default_values',
    ordering_strategy = 'pvalue_foldchange',
    minimum_markers = 4,
    percentile_markers = NULL,
    min_delta_cor_threshold = 0.05,
    percentile_markers.min_corr = NULL,
    verbose = FALSE){

  message(paste0('run_deconvolution_tree_guided_recursive - Verbose: ',
                 verbose))

  # Let's begin with a basic validation about the samples names staring with #
  current_columns <- unique(single_cell_data_exp[[sample]])
  column_validation <- all(
    unlist(lapply(current_columns,
                  function(text) !grepl("^[[:digit:]]+", text))))
  if(!column_validation){
    stop(paste0('There are samples in the single cell data that start with a',
                ' number. Please make sure that all samples start with a letter. Some R ',
                'functionalities do not work with this kind of naming conventions.'))
  }

  # Creation of the level name that I will use for the top deconvolution
  top_clusters_var <- paste0('metacluster_level_', hierarchy$tree_level)

  if(verbose){
    message('Recursive version tree-guided method: level ',
            hierarchy$tree_level, ', leaf ', hierarchy$leaf_number,
            ', strategy: ', use_min_cor_strategy)
  }

  #-1. Checking the parameters in the object to locally use it
  param.logfc.threshold.use <- param.logfc.threshold
  if(!is.null(hierarchy$parameters$param.logfc.threshold)){
    param.logfc.threshold.use <- hierarchy$parameters$param.logfc.threshold
    message('Using param.logfc.threshold from object hierarchy: ',
            param.logfc.threshold)
  }

  param.p_val_adj.use <- param.p_val_adj
  if(!is.null(hierarchy$parameters$param.p_val_adj)){
    param.p_val_adj.use <- hierarchy$parameters$param.p_val_adj
    message('Using param.p_val_adj from object hierarchy: ',
            param.p_val_adj.use)
  }

  minimum_markers.use <- minimum_markers
  if(!is.null(hierarchy$parameters$minimum_markers)){
    minimum_markers.use <- hierarchy$parameters$minimum_markers
    message('Using minimum_markers from object hierarchy: ',
            minimum_markers)
  }

  min_delta_cor_threshold.use <- min_delta_cor_threshold
  if(!is.null(hierarchy$parameters$min_delta_cor_threshold)){
    min_delta_cor_threshold.use <- hierarchy$parameters$min_delta_cor_threshold
    message('Using min_delta_cor_threshold from object hierarchy: ',
            min_delta_cor_threshold)
  }



  # 0. Create the top clustering on the dataset in a given variable or in the
  # final configuration if next_level_clustering is NULL
  if(!is.null(hierarchy$next_level_clustering)){
    single_cell_data_exp <- add_grouping_single_cell(
      single_cell_data_exp = single_cell_data_exp,
      sub_clusters_var = sub_clusters_var,
      level_group_name = top_clusters_var,
      hierarchy_level_groups = hierarchy$next_level_clustering)

    # I got the list of top clusters
    top_clusters_list <-
      unique(stats::na.omit(single_cell_data_exp[[top_clusters_var]]))
  }else{
    # Because this is the end of the tree, I will use the original cell type
    # variable and the original names as well
    top_clusters_list <- unique(stats::na.omit(hierarchy$celltype_list))
    top_clusters_var <- sub_clusters_var
  }

  # 1. Create the reference for the top clustering configuration
  reference_w.top.cluster <- decoflex_build_cell_reference(
    x = single_cell_data_exp,
    ct.sub = top_clusters_list,
    ct.varname = top_clusters_var,
    sample = sample,
    verbose = verbose)

  # 2. Calculate the Marker genes for the top clusters
  markers.top.clusters.object <- calculate_markers(
    single_cell_data_exp = single_cell_data_exp,
    reference = reference_w.top.cluster$basis,
    group_clusters_var = top_clusters_var,
    use_min_cor_strategy = use_min_cor_strategy,
    delete_shared_internal_markers = delete_shared_internal_markers,
    filter_markers = filter_markers,
    param.logfc.threshold = param.logfc.threshold.use,
    param.p_val_adj = param.p_val_adj.use,
    test.use.value = test.use.value,
    marker_strategy = marker_strategy,
    ordering_strategy = ordering_strategy,
    minimum_markers = minimum_markers.use,
    percentile_markers = percentile_markers,
    min_delta_cor_threshold = min_delta_cor_threshold.use,
    percentile_markers.min_corr = percentile_markers.min_corr,
    verbose = verbose)

  return(markers.top.clusters.object)
}

#' Add Grouping to Single Cell Data
#'
#' The `add_grouping_single_cell` function iterates over a list of cell types
#' and assigns labels to each single cell in a dataset based on its group at
#' each level. The labels are in the format 'subcluster_level_#_leaf_#'.
#'
#' @param single_cell_data_exp A list containing the single cell data. The
#'  function will add new labels to this list.
#' @param sub_clusters_var A string representing the variable used to indicate
#' the subclusters in `single_cell_data_exp`.
#' @param level_group_name A string representing the name of the level group
#'  variable in `single_cell_data_exp`. The function will add labels to this
#'  variable.
#' @param hierarchy_level_groups A list of cell types that define the hierarchy
#' levels. Each cell type should be a list itself with the elements
#' `celltype_list`, `tree_level`, and `leaf_number`.
#'
#' @return The function returns the `single_cell_data_exp` list with added
#' labels to the `level_group_name` variable.
#'
#'
#' @export
add_grouping_single_cell <- function(single_cell_data_exp,
                                     sub_clusters_var,
                                     level_group_name,
                                     hierarchy_level_groups){

  # Iteration on the cell type list.
  for (level in hierarchy_level_groups) {
    # For each group on the level, I add the labeling on the single cell with
    # a standard name: subcluster_level_#_leaf_#
    single_cell_data_exp[[level_group_name]][single_cell_data_exp[[sub_clusters_var]] %in%
                                               level$celltype_list] <-
      paste0('subcluster_level_',
             level$tree_level,
             '_leaf_',
             level$leaf_number)
  }

  # Returning the single cell data with the new labeling of clustering
  return(single_cell_data_exp)
}

#' Tree-guided Deconvolution of Bulk Data
#'
#' `run_deconvolution_tree_guided` is a generic function for deconvoluting bulk
#' data based on a reference single-cell dataset using a tree-guided method.
#'
#' @param bulk_data A numeric matrix representing the bulk data to be
#'  deconvoluted.
#' @param single_cell_data_exp A list or data.frame representing the single cell
#'  data.
#' @param top_clusters_var A string representing the variable in
#'  `single_cell_data_exp` used to indicate the top clusters.
#' @param top_clusters_list A list or vector indicating the top clusters to be
#'  considered.
#' @param sub_clusters_var A string representing the variable in
#'  `single_cell_data_exp` used to indicate the subclusters.
#' @param sub_clusters_list A list or vector indicating the subclusters to be
#'  considered.
#' @param sample A string or numeric value indicating the specific sample to
#'  use.
#' @param use_min_cor_strategy A logical value indicating whether to use the
#'  minimal correlation strategy. Default is FALSE.
#' @param true_proportions A numeric vector representing the true proportions of
#'  the clusters, if known. Default is NULL.
#' @param delete_shared_level_markers A logical value indicating whether to
#'  delete shared markers at the same level. Default is TRUE.
#' @param delete_shared_internal_markers A logical value indicating whether to
#'  delete shared markers within a group. Default is TRUE.
#' @param filter_markers A character string specifying the filtering strategy
#'  for markers. Default is NULL.
#' @param param.logfc.threshold A numeric value representing the log fold change
#'  threshold for marker genes. Default is 2.0.
#' @param param.p_val_adj A numeric value representing the adjusted p-value
#'  threshold for marker genes. Default is 0.05.
#' @param test.use.value A string specifying the test to use for marker gene
#'  identification. Default is 'wilcox'.
#' @param marker_strategy A character string specifying the overall strategy
#'  for marker gene handling. Default is NULL.
#' @param minimum_markers The minimum number of marker genes to be selected,
#'  default is set to 4. This can be adjusted based on the analysis
#'  requirements.
#' @param min_delta_cor_threshold Minimum difference between the correlation of
#'  t and t-1 steps. The idea is to skip genes that have small changes in
#'  the correlation. For default it is a 5%.
#' @param verbose A logical value indicating whether to print detailed output
#'  during execution. Default is FALSE.
#'
#' @return A list containing the deconvoluted bulk data.
#'
#' @export
run_deconvolution_tree_guided <- function(bulk_data, single_cell_data_exp,
                                          top_clusters_var, top_clusters_list,
                                          sub_clusters_var, sub_clusters_list,
                                          sample,
                                          use_min_cor_strategy = FALSE,
                                          true_proportions = NULL,
                                          delete_shared_level_markers = TRUE,
                                          delete_shared_internal_markers = TRUE,
                                          filter_markers = NULL,
                                          param.logfc.threshold = 2.0,
                                          param.p_val_adj = 0.05,
                                          test.use.value = 'wilcox',
                                          marker_strategy = 'keep_default_values',
                                          minimum_markers = 4,
                                          percentile_markers = NULL,
                                          min_delta_cor_threshold = 0.05,
                                          percentile_markers.min_corr = NULL,
                                          max_iterations = 10000,
                                          delta_threshold = 1e-10,
                                          verbose = FALSE){

  if(verbose){
    print(paste0('Use minimal correlation strategy: ', use_min_cor_strategy))
  }

  # 1. Create the reference for the top clustering configuration
  reference_w.top.cluster <- decoflex_build_cell_reference(
    x = single_cell_data_exp,
    ct.sub = top_clusters_list,
    ct.varname = top_clusters_var,
    sample = sample,
    verbose = verbose)

  # 3. Calculate the Marker genes for the top clusters
  markers.top.clusters.object <- calculate_markers(
    single_cell_data_exp = single_cell_data_exp,
    reference = reference_w.top.cluster$basis,
    group_clusters_var = top_clusters_var,
    use_min_cor_strategy = use_min_cor_strategy,
    delete_shared_internal_markers = delete_shared_internal_markers,
    filter_markers = filter_markers,
    param.logfc.threshold = param.logfc.threshold,
    param.p_val_adj = param.p_val_adj,
    marker_strategy = marker_strategy,
    minimum_markers = minimum_markers,
    percentile_markers = percentile_markers,
    min_delta_cor_threshold = min_delta_cor_threshold,
    percentile_markers.min_corr = percentile_markers.min_corr,
    verbose = verbose)

  # 3.1. Assing the final markers and the top marker list
  markers.top.clusters <- markers.top.clusters.object$total_markers
  list.top.marker <- markers.top.clusters.object$list_markers

  # 2. Run de deconvolution of the top clusters using the markers genes for the
  # top clusters.
  result_deco_top_cluster <- run_standard_deconvolution(
    bulk_data_x = bulk_data,
    references_w = reference_w.top.cluster$basis,
    markers = markers.top.clusters,
    max_iterations = max_iterations,
    delta_threshold = delta_threshold,
    verbose = verbose)

  # In case the process is a simulation
  if(!is.null(true_proportions)){
    # 3.a. Calculate the relative proportions with respect to the subset of
    # subcluster taken.
    true_prop_relative <-
      as.data.frame.matrix(true_proportions[, sub_clusters_list])
    # rescaling of the proportions relative JUST TO THE SUBCLUSTERS
    true_prop_relative <- true_prop_relative/rowSums(true_prop_relative)


    # 3.b calculate the proportions for the top model
    true_prop_top_clustering <- aggregate_true_top_clustering_prop(
      single_cell_data_exp = single_cell_data_exp,
      true_prop = true_prop_relative,
      top_clusters_var = top_clusters_var,
      top_clusters_list = top_clusters_list,
      sub_clusters_var = sub_clusters_var)

    # 4. Now I can calculate the metrics for both configurations: samples and
    # celltypes
    pm_top_clustering <- performance_metrics_general(
      true_prop = true_prop_top_clustering,
      calc_prop = result_deco_top_cluster$h,
      title = 'Top clustering',
      verbose = verbose)

    # Adding the metric information to the result of the top model
    result_deco_top_cluster[['simulation_metrics']] <- pm_top_clustering
  }


  # Create a variable to save the results of each deconvolution
  deco_results_list <- list()

  # Check the marker genes for the level to have the shared level markers
  # If the parameter is false, I will avoid this calculation.
  shared_level_markers = NULL
  if(delete_shared_level_markers){
    shared_level_markers <- calculate_markers_level_intersection(
      single_cell_data_exp = single_cell_data_exp,
      top_clusters_var = top_clusters_var,
      sub_clusters_var = sub_clusters_var,
      top_clusters_list = top_clusters_list,
      param.logfc.threshold = param.logfc.threshold,
      param.p_val_adj = param.p_val_adj,
      marker_strategy = marker_strategy,
      verbose = verbose)
  }

  # 5. For each top cluster create the deconvolution for the subclusters inside
  for (top_cluster in top_clusters_list) {

    # 5.1. check the subcluster that are part of the top cluster analized.
    values_current_top_cluster <-
      single_cell_data_exp[[top_clusters_var]] == top_cluster
    values_subclusters <- single_cell_data_exp[[sub_clusters_var]]
    subclusters_list <- unique(values_subclusters[values_current_top_cluster])
    subclusters_list <-
      levels(factor(subclusters_list[!is.na(subclusters_list)]))

    # 5.2. Create the reference for the subclusters
    reference_w_subclusters <- decoflex_build_cell_reference(
      x = single_cell_data_exp,
      ct.sub = subclusters_list,
      ct.varname = sub_clusters_var,
      sample = sample,
      verbose = verbose)

    # 5.3. Calculate the marker genes for the separated deconvolution
    markers_subclusters.object <- calculate_markers(
      single_cell_data_exp = single_cell_data_exp,
      reference = reference_w_subclusters$basis,
      group_clusters_var = sub_clusters_var,
      use_min_cor_strategy = use_min_cor_strategy,
      delete_shared_level_markers = delete_shared_level_markers,
      shared_level_markers = shared_level_markers,
      delete_shared_internal_markers = delete_shared_internal_markers,
      filter_markers = filter_markers,
      param.logfc.threshold = param.logfc.threshold,
      param.p_val_adj = param.p_val_adj,
      test.use.value = test.use.value,
      marker_strategy = marker_strategy,
      minimum_markers = minimum_markers,
      percentile_markers = percentile_markers,
      min_delta_cor_threshold = min_delta_cor_threshold,
      percentile_markers.min_corr = percentile_markers.min_corr,
      verbose = verbose)

    # 5.3.1. Assing the final markers and the top marker list
    markers_subclusters <- markers_subclusters.object$total_markers
    list.subcluster.marker <- markers_subclusters.object$list_markers

    # 5.4. Calculate the deconvolution with the marker genes.
    result_deco_subclusters <- run_standard_deconvolution(
      bulk_data_x = bulk_data,
      references_w = reference_w_subclusters$basis,
      markers = markers_subclusters,
      max_iterations = max_iterations,
      delta_threshold = delta_threshold,
      verbose = verbose)

    # 5.4.1. Adding detailed info of clusters
    result_deco_subclusters[['markers_subclusters.object']] <-
      list.subcluster.marker

    # In case that the current process is a simulation
    if(!is.null(true_proportions)){
      # 5.5. Recalculate the proportions for each final model.
      deco_results_current_list <- list()
      deco_results_current_list[[top_cluster]] <- result_deco_subclusters
      partial_prop_top_cluster <- recalculate_proportions(
        top_cluster_deco_list = result_deco_top_cluster,
        top_clusters_list = top_cluster,
        subcluster_deco_list = deco_results_current_list,
        verbose = verbose)

      # 5.6. Calculate metrics for each deconvolution based on the true
      # proportions.
      pm_subclustering <- performance_metrics_general(
        t(true_prop_relative[, subclusters_list]),
        partial_prop_top_cluster,
        title = paste0('Top cluster: ', top_cluster),
        verbose = verbose)

      # 5.7.a. adding the metric information to the result of the top model
      result_deco_subclusters[['simulation_metrics']] <- pm_subclustering
    }

    # 5.7.b. Add the results to the global results.
    deco_results_list[[top_cluster]] <- result_deco_subclusters

  }

  # 6. Recalculate proportions
  total_prop <- recalculate_proportions(
    top_cluster_deco_list = result_deco_top_cluster,
    top_clusters_list = top_clusters_list,
    subcluster_deco_list = deco_results_list,
    verbose = verbose)

  # In case the current process is a simulation.
  if(!is.null(true_proportions)){
    # 7. Calculate metrics of final proportions.
    # I need to convert this proportions to a normal dataset.
    true_prop_relative <- true_prop_relative[, sub_clusters_list]
    print('FINAL CALCULATION...')
    pm_final_model <- performance_metrics_general(
      true_prop = true_prop_relative,
      calc_prop = t(total_prop),
      title = paste0('Complete model'),
      verbose = verbose)
  }else
  {
    pm_final_model =  NULL
  }


  # 8. Return the final variables.
  return(list(proportions = total_prop,
              result_deco_top_cluster = result_deco_top_cluster,
              deco_results_list = deco_results_list,
              metric_final_model = pm_final_model,
              markers.top.clusters.object = markers.top.clusters.object))
}

#' Calculation of Marker Genes Intersection at Different Levels
#'
#' `calculate_markers_level_intersection` computes the intersection of marker
#' genes at various levels of clustering. This function is currently under
#' testing and development, intending to streamline the code and improve the
#' algorithm performance.
#'
#' @param single_cell_data_exp A list or data.frame representing the single cell
#'  data.
#' @param top_clusters_var A string representing the variable in
#'  `single_cell_data_exp` used to indicate the top clusters.
#' @param sub_clusters_var A string representing the variable in
#'  `single_cell_data_exp` used to indicate the subclusters.
#' @param top_clusters_list A list or vector indicating the top clusters to be
#'  considered.
#' @param param.logfc.threshold A numeric value representing the log fold change
#'  threshold for marker genes. Default is 2.0.
#' @param param.p_val_adj A numeric value representing the adjusted p-value
#'  threshold for marker genes. Default is 0.05.
#' @param marker_strategy A character string specifying the overall strategy for
#'  marker gene handling. Default is NULL.
#' @param minimum_markers The minimum number of marker genes to be selected,
#'  default is set to 4. This can be adjusted based on the analysis
#'  requirements.
#' @param min_delta_cor_threshold Minimum difference between the correlation of
#'  t and t-1 steps. The idea is to skip genes that have small changes in
#'  the correlation. For default it is a 5%.
#' @param verbose A logical value indicating whether to print detailed output
#'  during execution. Default is FALSE.
#'
#' @return A list containing the intersection of marker genes at different
#' levels.
#'
#' @export
calculate_markers_level_intersection <- function(single_cell_data_exp,
                                                 top_clusters_var,
                                                 sub_clusters_var,
                                                 top_clusters_list,
                                                 param.logfc.threshold = 2.0,
                                                 param.p_val_adj = 0.05,
                                                 marker_strategy = NULL,
                                                 minimum_markers = 4,
                                                 percentile_markers = NULL,
                                                 min_delta_cor_threshold = 0.05,
                                                 percentile_markers.min_corr = NULL,
                                                 verbose = FALSE){
  # 1. For each top cluster create the deconvolution for the subclusters
  # inside.
  total_markers <- NULL
  intersection_markers <- NULL
  for (top_cluster in top_clusters_list) {

    # 1.1. check the subcluster that are part of the top cluster analized.
    values_current_top_cluster <-
      single_cell_data_exp[[top_clusters_var]] == top_cluster
    values_subclusters <- single_cell_data_exp[[sub_clusters_var]]
    subclusters_list <- unique(values_subclusters[values_current_top_cluster])
    subclusters_list <-
      levels(factor(subclusters_list[!is.na(subclusters_list)]))


    reference_w_subclusters = NULL
    markers_subclusters = NULL
    # when we have a cluster with just 1 cell type, I CAN'T calculate a
    # reference model or a marker gene set.
    if(length(subclusters_list) > 1){
      # 1.2. Create the reference for the subclusters
      reference_w_subclusters <- decoflex_build_cell_reference(
        x = single_cell_data_exp,
        ct.sub = subclusters_list,
        ct.varname = sub_clusters_var,
        sample = sample,
        verbose = verbose)

      # 1.3. Calculate the marker genes for the separated deconvolution not
      # using the marker strategy.
      markers_subclusters.object <- calculate_markers(
        single_cell_data_exp = single_cell_data_exp,
        reference = reference_w_subclusters$basis,
        group_clusters_var = sub_clusters_var,
        use_min_cor_strategy = FALSE,
        param.logfc.threshold = param.logfc.threshold,
        param.p_val_adj = param.p_val_adj,
        marker_strategy = marker_strategy,
        minimum_markers = minimum_markers,
        percentile_markers = percentile_markers,
        min_delta_cor_threshold = min_delta_cor_threshold,
        percentile_markers.min_corr = percentile_markers.min_corr,
        verbose =  verbose)

      # 1.3.1. Assing the final markers and the top marker list
      markers_subclusters <- markers_subclusters.object$total_markers
      list.subcluster.marker <- markers_subclusters.object$list_markers

      # 1.4. Check the intersection between the markers to avoid shared genes.
      if(is.null(total_markers)){
        total_markers <- markers_subclusters
      }else{
        # I add more marker to the intersection.
        intersection_markers <-
          unique(c(intersection_markers,
                   intersect(total_markers, markers_subclusters)))
        total_markers <- unique(c(total_markers, markers_subclusters))
      }
    }else{
      if(verbose){
        message(
          'Group with 1 cell type. The process of ',
          'calculate_markers_level_intersection will not calculate a reference ',
          'and markers genes because it is not possible witht the cluster: ',
          subclusters_list)
      }
    }
  }

  # 2. I should check the shared genes
  print(paste0('Calculating shared level marker genes: ',
               length(intersection_markers)))

  # 3. Return the list of shared marker genes for the different groups in the
  # level.
  return(intersection_markers)
}

#' Rescale Clustering Proportions Based on Higher Level Data
#'
#' `calculate_prop_back_propagation` function computes rescaled proportions of
#' higher level clustering based on more detailed proportions. It assumes the
#' input `top_proportions` represent proportions for single-leaf samples in the
#' deconvolution.
#'
#' @param top_proportions A numeric or matrix of proportions representing higher
#'  level clustering.
#' @param level_proportions A numeric or matrix of proportions representing more
#'  detailed clustering.
#' @param verbose A logical value indicating whether to print detailed output
#'  during execution. Default is FALSE.
#'
#' @return A numeric or matrix of rescaled proportions.
#'
#' @details
#' The function applies a specific calculation depending on whether
#' `top_proportions` is numeric or a matrix. If it's numeric, it simply
#' multiplies `level_proportions` by `top_proportions`.
#' If it's a matrix, it applies a sweep operation to multiply each column of
#' `level_proportions` by the corresponding column in `top_proportions`. The
#' resulting proportions should sum up to the same values as `top_proportions`.
#'
#' @export
calculate_prop_back_propagation <- function(top_proportions,
                                            level_proportions,
                                            verbose = FALSE){
  # Just in case that we have just one sample and one cell type or two. In that
  # case the top_proportions will be just a number without columns or rows.
  rescaled_proportions = NULL

  if(is.numeric(top_proportions)){
    rescaled_proportions <- level_proportions * top_proportions
  }else{
    rescaled_proportions <- sweep(
      level_proportions,
      MARGIN=2,
      t(top_proportions[, colnames(level_proportions)]), `*`)
  }

  # Since the top is the reference proportions, the rescaled version must sum
  # up to the same values
  # colSums(rescaled_proportions) == top_proportions
  return(rescaled_proportions)
}

#' Join Back-Propagation Proportions from Multiple Deconvolution Results
#'
#' `join_back_propagation_proportions_leaves` function iterates over a list of
#' deconvolution results and merges the back-propagated proportions from each
#' result.
#'
#' @param deco_results_list A list containing deconvolution results. Each
#'  element in the list shouldbe a list with an element named
#'  "back_propagation_proportions_top_detailed" that contains the
#'  back-propagated proportions.
#' @param verbose A logical value indicating whether to print detailed output
#'  during execution. Default is FALSE.
#'
#' @return A matrix that contains the merged back-propagation proportions from
#' all the input deconvolution results.
#'
#' @details
#' The function iterates over the deconvolution results in `deco_results_list`.
#' From each result, it extracts the back-propagation proportions, then merges
#' them row-wise using the `rbind` function.
#' This process continues until all the back-propagation proportions from all
#' deconvolution results have been merged.
#'
#' @export
join_back_propagation_proportions_leaves <- function(deco_results_list,
                                                     verbose = FALSE){

  # iterate in the list of deconvolutions
  calculated_proportions <- NULL
  deco_result_list_names <- names(deco_results_list)
  for (cluster_counter in 1:length(deco_results_list)) {
    cluster <- deco_results_list[cluster_counter]
    cluster_name <- deco_result_list_names[cluster_counter]

    # I get the rescaled proportions from the current cluster. Since always is
    # a top proportion, I have to call the specific name.
    current_cluster_proportions <-
      cluster[[cluster_name]]$back_propagation_proportions_top_detailed

    # Integration of the rescaled proportions
    if(is.null(calculated_proportions)){
      calculated_proportions <- current_cluster_proportions
    }else{
      calculated_proportions <- rbind(calculated_proportions,
                                      current_cluster_proportions)
    }
  }

  calculated_proportions
}


#' Aggregate True Proportions for Top-Level Clustering
#'
#' `aggregate_true_top_clustering_prop` function calculates the true proportions
#' of top-level clustersbased on the true subclustering proportions.
#'
#' @param single_cell_data_exp A data frame or matrix containing the single
#'  cell data.
#' @param true_prop A matrix or data frame containing the true proportions of
#'  the subclusters.
#' @param top_clusters_var The variable name of the top clusters in the
#'  `single_cell_data_exp`.
#' @param top_clusters_list A list containing the identifiers of the top
#'  clusters.
#' @param sub_clusters_var The variable name of the subclusters in the
#'  `single_cell_data_exp`.
#'
#' @return A matrix that contains the aggregated true proportions of the top
#' clusters.
#'
#' @details
#' The function iterates over the top-level clusters in `top_clusters_list`.
#' For each top-level cluster, it identifies the corresponding subclusters and
#' calculates the sum of their true proportions.
#' The process continues until the aggregated true proportions for all top-level
#'  clusters have been calculated.
#'
#' @export
aggregate_true_top_clustering_prop <- function(single_cell_data_exp,
                                               true_prop,
                                               top_clusters_var,
                                               top_clusters_list,
                                               sub_clusters_var){


  # 1. Iterate on the set of top cluster and then check the total True
  # proportions for each one.
  aggregated_prop <- NULL
  for (top_cluster in top_clusters_list) {

    # 1.1. check the subcluster that are part of the top cluster analized.
    values_current_top_cluster <-
      single_cell_data_exp[[top_clusters_var]] == top_cluster
    values_subclusters <- single_cell_data_exp[[sub_clusters_var]]
    subclusters_list <- unique(values_subclusters[values_current_top_cluster])
    subclusters_list <-
      levels(factor(subclusters_list[!is.na(subclusters_list)]))

    # 1.2. Let's calculate the sum of the proportions for the set of cluster
    # that are into the top cluster.
    # It is possible that the actual group is composed by one cell type,
    # therefore, it will be a list not a matrix
    if(length(subclusters_list)==1){
      sum_subcluster_top <- true_prop[, subclusters_list]
    }else{
      sum_subcluster_top <- rowSums(true_prop[, subclusters_list])
    }

    # 1.3. Add the summarized proportion to the list.
    if(is.null(aggregated_prop)){
      aggregated_prop <- rbind(sum_subcluster_top)
    }else{
      aggregated_prop <- rbind(aggregated_prop, sum_subcluster_top)
    }
  }

  # 2. assign the rownames to the top proportion summary
  colnames(aggregated_prop) <- rownames(true_prop)
  rownames(aggregated_prop) <- top_clusters_list

  return(aggregated_prop)
}


#' Calculate Marker Genes
#'
#' The `calculate_markers` function identifies marker genes for specific cell
#' types based on given criteria and strategies. These strategies help ensure
#' that the marker genesidentified are more unique to each cell type, providing
#' better accuracy in downstream analysis.
#'
#' @param single_cell_data_exp A data frame or matrix containing the single cell
#'  data.
#' @param reference A data frame or matrix containing the reference data.
#' @param group_clusters_var The variable name of the group clusters in the
#'  `single_cell_data_exp`.
#' @param use_min_cor_strategy Logical value indicating whether to use the
#'  minimum correlation strategy. This strategy is used to identify marker genes
#'  that are least correlated with each other, thereby ensuring uniqueness.
#'  Default is FALSE.
#' @param delete_shared_level_markers Logical value indicating whether to
#'  delete shared level markers. This parameter enables the removal of marker
#'  genes shared across different levels to ensure uniqueness within each level.
#'  Default is TRUE.
#' @param shared_level_markers A vector containing shared level markers.
#'  NULL by default.
#' @param delete_shared_internal_markers Logical value indicating whether to
#'  delete shared internal markers. This parameter allows for the removal of
#'  marker genes that are common within a given level, to promote marker
#'  uniqueness. Default is TRUE.
#' @param p_value_attribute String indicating the attribute name of the p-value
#'  in the single cell data. Default is 'p_val_adj'.
#' @param param.logfc.threshold Threshold for the log-fold change, determining
#'  the minimum acceptable change for a gene to be considered a marker.
#'  Default is 2.0.
#' @param param.p_val_adj Adjusted p-value threshold. This threshold is used to
#'  control the false discovery rate. Default is 0.05.
#' @param test.use.value String indicating the statistical test used for
#' identifying markers. Default is 'wilcox'.
#' @param filter_markers A vector containing markers to be filtered out before
#' the analysis. NULL by default.
#' @param marker_strategy Strategy used to select marker genes. This could be a
#' specific algorithm or statistical method. NULL by default.
#' @param minimum_markers The minimum number of marker genes to be selected,
#'  default is set to 4. This can be adjusted based on the analysis
#'  requirements.
#' @param min_delta_cor_threshold Minimum difference between the correlation of
#'  t and t-1 steps. The idea is to skip genes that have small changes in
#'  the correlation. For default it is a 5%.
#' @param verbose Logical value indicating whether to display additional
#' information during the function's execution. Default is TRUE.
#'
#' @return A list of marker genes for each cell type.
#'
#' @import Seurat
#' @import plyr
#'
#' @details
#' The function calculates marker genes for each cell type in the provided data,
#' with options to apply different strategies and filters during the process.
#' It considers factors such as p-value, log-fold change, and uniqueness at
#' different levels to ensure the identified markers are most representative of
#' their respective cell types.
#'
#' @seealso \code{\link{Matrix}}
#'
#' @export
calculate_markers <- function(single_cell_data_exp,
                              reference,
                              group_clusters_var,
                              use_min_cor_strategy = FALSE,
                              delete_shared_level_markers = TRUE,
                              shared_level_markers = NULL,
                              delete_shared_internal_markers = TRUE,
                              p_value_attribute = 'p_val_adj',
                              param.logfc.threshold = 2.0,
                              param.p_val_adj = 0.05,
                              test.use.value = 'wilcox',
                              ordering_strategy = 'pvalue_foldchange',
                              filter_markers = NULL,
                              marker_strategy = 'keep_default_values',
                              minimum_markers = 4,
                              percentile_markers = NULL,
                              min_delta_cor_threshold = 0.05,
                              percentile_markers.min_corr = NULL,
                              verbose = TRUE){

  print(paste0('calculate_markers with verbose = ', verbose))

  # 0. To check the calculate_markers function input details.
  if(verbose){
    message('Parameters for calculate_markers function:')
    message('p_value_attribute: ', p_value_attribute)
    message('param.logfc.threshold: ', param.logfc.threshold)
    message('param.p_val_adj: ', param.p_val_adj)
    message('test.use.value: ', test.use.value)
    message('minimum_markers: ', minimum_markers)
    message('min_delta_cor_threshold: ', min_delta_cor_threshold)
  }

  # 1. Create seurat object
  seurat.reference <- Seurat::CreateSeuratObject(
    counts = single_cell_data_exp@assayData$exprs,
    project = "deconvolution_bulk_data",
    assay = "RNA")
  Seurat::Idents(seurat.reference) <- single_cell_data_exp[[group_clusters_var]]

  # We need to assing the assay and after normalize the data in order to use
  # MAST.
  if(test.use.value == 'MAST'){
    Seurat::DefaultAssay(seurat.reference) = "RNA"
    seurat.reference <- Seurat::NormalizeData(object = seurat.reference,
                                      normalization.method = 'LogNormalize')
  }

  # 2. For each group (top clusters or subcluster) I need to check the markers
  list_groups <- colnames(reference)
  intersection_markers <- NULL
  total_markers <-  NULL
  list_markers <- list()
  for (group in list_groups) {

    # 2.1. I need to calculate the rest of the groups to compare with.
    rest_groups <- list_groups[!(list_groups %in% c(group))]

    # To check the marker selection details.
    if(verbose){
      message('Finding markers (FindMarkers) for: group = ',
              group,', rest_groups = ',
              paste(shQuote(rest_groups), collapse=", "))
    }

    # 2.2 In case the FindMarkers function can't recover markers for a
    # specific group.
    markers.info <-  NULL
    markers.info.all <-  NULL
    tryCatch(
      {

        # 2.2.1 Now we can find the markers with the limited parameters
        # When the single cell data is big, let's say 150.000 cells, and then
        # the number of groups is two or three, the p-values can be NA if
        # max.cells.per.ident = Inf, therefore it is a good idea to limit this
        # to 1 third of the total of cells, or a high values as 50.000 in this
        # case. By default the value is infinite
        markers.info.all <- FindMarkers.proxy(
          seurat.reference = seurat.reference,
          group = group, rest_groups = rest_groups,
          test.use.value = test.use.value,
          param.logfc.threshold = param.logfc.threshold,
          marker_strategy = marker_strategy,
          verbose = verbose)

        markers.info <-  markers.info.all

        if(verbose){
          print(paste0('Number of NA: ', sum(is.na(markers.info))))
          print(paste0('Number of markers without filter: ',
                       nrow(markers.info)))
        }

      },
      error = function(e){
        # In case that the test in order to get markers doesn't reply any
        # marker based on the threshold or another error.
        warning(paste0('The FindMarker function is failing in the group: ',
                       group, '. The error is: ', e))

        if(!is.null(e)){
          markers.info <- data.frame()
        }else{
          stop('Error running the FindMarkers function.', e)
        }
      }
    )

    # 2.3. To filter the genes,
    if(p_value_attribute == 'p_val'){
      markers.info <- subset(markers.info,
                             markers.info$p_val < param.p_val_adj)
    }else{
      markers.info <- subset(markers.info,
                             markers.info$p_val_adj < param.p_val_adj)
    }

    # 2.2.1. Before doing something else, I should check if the object is
    # empty or null
    if(is.null(markers.info) | nrow(markers.info)==0){
      message("We couldn't find any marker for this combination of ",
              "groups/clusters")
      next
    }else{
      # The process continues and add these markerts to the analysis.
      print(paste0('is.null(markers.info): ', is.null(markers.info)))
      print(paste0('nrow(markers.info)==0: ', nrow(markers.info)==0))
      print(paste0('nrow(markers.info): ', nrow(markers.info)))
    }

    # Depend of the Seurat version we will find avg_logFC or avg_log2FC. I
    # wrote all cases to make the code more readable.
    foldchange_var <- 'avg_logFC'
    if(foldchange_var %in% colnames(markers.info)){
      foldchange_var <- 'avg_logFC'
    }else{
      foldchange_var <- 'avg_log2FC'
    }

    # I need to order first the p-values from least to greatest and then the
    # foldchange from biggest to smallest. I will use the variable
    # ordering_strategy [pvalue_foldchange, foldchange_pvalue]

    # Since arrange losses the rownames, I added it as columsn to fix the bug
    markers.info$gene_name <- rownames(markers.info)
    message(paste0("ordering_strategy", ordering_strategy))
    if(ordering_strategy == 'pvalue_foldchange'){

      # 1a. Possible solution 1: first order by p-values and then for
      # foldchange value.
      markers.info <- plyr::arrange(markers.info,
                              (markers.info[[p_value_attribute]]),
                              plyr::desc(markers.info[[foldchange_var]]))

    }else{
      # 1b. Possible solution 2: second order by foldchange and then p-values
      markers.info <- plyr::arrange(markers.info,
                                    plyr::desc(markers.info[[foldchange_var]]),
                                    (markers.info[[p_value_attribute]]))
    }

    #Anyway, I return the gene names to the df after the arrange.
    rownames(markers.info) <- markers.info$gene_name


    markers.info$group_name <- paste0('markers_', group)

    if(delete_shared_level_markers & !is.null(shared_level_markers)){

      if(verbose){
        print(paste0('1. Eliminating shared_level_markers: ',
                     length(shared_level_markers),
                     ' - current: ', nrow(markers.info)))
      }

      # Deleting the shared genes
      markers.info  <-
        markers.info[!(rownames(markers.info) %in% c(shared_level_markers)),]

      if(verbose){
        print(paste0('2. Eliminating shared_level_markers: ',
                     length(shared_level_markers),
                     ' - after: ',
                     nrow(markers.info)))
      }
    }

    # 2.3 Let's add the (key, value) pair to the list
    if(!is.null(markers.info$group_name)){
      list_markers[[ group ]] <- markers.info
    }

    # 2.4. Check the intersection between the markers to avoid shared genes.
    if(is.null(total_markers)){
      total_markers <- rownames(markers.info)
    }else{
      # I add more marker to the intersection.
      intersection_markers <- c(intersection_markers,
                                intersect(total_markers,
                                          rownames(markers.info)))
      total_markers <- c(total_markers, rownames(markers.info))
    }

    # 2.5. Printing of the total markers for the current group.
    if(verbose){
      message(paste0('Number of markers for the group: ', nrow(markers.info)))
    }
  }

  # 2.1. After calculating the markers, I will check if markers are zero
  if(length(total_markers)==0){
    stop('There are no marker genes identified for the current process. ',
         'Please verify the value assigned to the parameter ',
         '"param.logfc.threshold". If necessary, consider using a smaller ',
         'value. The current value is: ', param.logfc.threshold)
  }

  # 3. Remove shared genes between clusters if there are shared markers and
  # the parameter is TRUE.
  if(!is.null(intersection_markers) & delete_shared_internal_markers){

    if(verbose){
      print(paste0('Deleting shared internal marker genes: ',
                   length(intersection_markers)))
    }

    # 3.1. Deleting shared internal markers
    total_markers <-
      total_markers[!(total_markers %in% c(intersection_markers))]

    # 3.2. also I have to remove then from each list of markers
    for(counter in 1:length(list_markers)) {
      current_list_markers <- list_markers[[counter]]
      list_markers[[counter]] <-
        current_list_markers[!(rownames(current_list_markers) %in%
                                 c(intersection_markers)),]
    }
  }

  # 4. If the parameter is active I  will filter the genes that are in the
  # array. This sometimes is due to bid differences between the original bulk
  # data and the single cell, therefore deleting a gene that is highly
  # expressed in the bulk data we eliminate the bias.
  if(!is.null(filter_markers)){
    print(paste0('Deleting filtered marker genes passed in the ',
                 'parameter(filter_markers): ',
                 paste(shQuote(filter_markers), collapse=", ")
                 , ' (', length(filter_markers), ')'))
    total_markers <- total_markers[!(total_markers %in% c(filter_markers))]

    # also I have to remove then from each list of markers
    for(counter in 1:length(list_markers)) {
      current_list_markers <- list_markers[[counter]]
      list_markers[[counter]] <-
        current_list_markers[!(rownames(current_list_markers) %in%
                                 c(filter_markers)),]
    }
  }

  # 5. Use percentile filtering normally .75 to extract for each celltype the top percentile of markers.
  list_markers.percentile <- NULL
  if(!is.null(percentile_markers)){

    if(percentile_markers > 1){
      error(paste0('percentile_markers must be (0,1]. Current value: ', percentile_markers))
    }else{
      list_markers.percentile <- list_markers
      total_markers.percentile = NULL
      print(paste0("Percentile filtering: ", percentile_markers))

      #For each celltype I choose the markers below the percentile threshold.
      for(celltype in names(list_markers.percentile)){
        quartile.calculation <- quantile(as.data.frame(list_markers.percentile[[celltype]]$avg_log2FC), c(percentile_markers), na.rm = TRUE)
        print(quartile.calculation[1])
        filtered.data <- list_markers.percentile[[celltype]][list_markers.percentile[[celltype]]$avg_log2FC >= quartile.calculation[1],]
        print(paste0(celltype, '=', nrow(filtered.data)))
        list_markers.percentile[[celltype]] <- filtered.data
        total_markers.percentile <- append(total_markers.percentile, unlist(filtered.data$gene_name))
      }

      total_markers.percentile <- unlist(total_markers.percentile)
      total_markers <- total_markers.percentile
    }
  }

  # 6. Check with min_cor strategy.
  if(use_min_cor_strategy){
    results.min.corr <- calculate_min_correlation_incremental_markers(
      list_markers = list_markers,
      reference = reference,
      minimum_markers = minimum_markers,
      min_delta_cor_threshold = min_delta_cor_threshold,
      percentile_markers.min_corr = percentile_markers.min_corr,
      verbose = verbose)

    total_markers <- results.min.corr$markers
  }else{
    results.min.corr <- list()
    results.min.corr$accumulative_corr_markers <- NULL
  }

  # 7. I get just unique markers
  total_markers <- unique(total_markers)

  if(verbose){
    print(paste0('Total markers: ', length(total_markers)))
  }

  # Return both objects: final markers and markers by cell-type
  return(list(total_markers = total_markers,
              list_markers = list_markers,
              list_markers.percentile = list_markers.percentile,
              accumulative_corr_markers = results.min.corr$accumulative_corr_markers))
}

#' Calculate Marker Genes Proxy
#'
#' The `FindMarkers.proxy` function is a flexible interface for identifying
#' marker genes using various methods.It's designed to be readily extensible,
#' facilitating the integration of novel marker calculation techniques in the
#' future. The function calculates the markers for a specific group relative to
#' all other cells (rest_groups) in a given Seurat object (seurat.reference).
#'
#' @param seurat.reference A Seurat object used as the reference for marker gene
#'  calculation. This object contains the single-cell RNA-seq data.
#' @param group A string specifying the cell group (based on clustering) for
#'  which to calculate marker genes.
#' @param rest_groups A list or vector specifying the remaining cell groups to
#'  be considered as background in the marker genes calculation.
#' @param test.use.value String specifying the statistical test to use for
#'  marker identification. The methods currently supported include:
#' 'wilcox', 'bimod', 't', 'negbinom', 'poisson', 'LR', and 'MAST'. Please note
#'  that 'roc' and 'DESeq2' methods are currently not supported due to
#'  discrepancies in output values and support issues, respectively. Default is
#'  'wilcox'.
#' @param marker_strategy A string specifying the strategy to use in marker
#'  gene calculation.
#'  'keep_default_values' (default) retains the default values. Other strategies
#'  may be added in the future.
#' @param max.size.percentage Maximum acceptable size (as a proportion) for a
#'  group in the marker gene calculation. Default is 0.70.
#' @param param.logfc.threshold Threshold for the log-fold change, which
#'  determines the minimum change required for a gene to be considered as a
#'  marker. Default is 2.0.
#' @param verbose Logical value indicating whether to display detailed progress
#'  messages during function execution. Default is FALSE.
#'
#' @return A list of marker genes for the specified cell group.
#'
#' @import Seurat
#'
#' @details
#' This function serves as a proxy to different methods for marker gene
#' calculation. It currently supports several statistical tests and provides a
#' convenient point of extension for the inclusion of additional marker
#' calculation methods. Note that the support for different tests can be
#' influenced by the characteristics of the single-cell experiments, and it's
#' important to choose a method that's appropriate for your data.
#'
#' @export
FindMarkers.proxy <- function(seurat.reference, group, rest_groups,
                              test.use.value = 'wilcox',
                              marker_strategy = 'keep_default_values',
                              max.size.percentage = 0.70,
                              param.logfc.threshold = 2.0, verbose = TRUE){

  # Not supported test methods.
  if(test.use.value == 'roc' | test.use.value == 'DESeq2'){
    stop(paste0('The test: ', test.use.value,
                ' is not supported (roc and DESeq2)'))
  }

  max.cells.per.ident <- Inf

  # 1. zero_approach = I cut the number of cells for the sampling
  # 2. convert_inf_to_zero  = I will convert all inf values to zero.
  if(marker_strategy == 'threshold_cut'){
    message('Using threshold_cut strategy to set a lower value of ',
            'max.cells.per.ident different to Inf or NA (Default).')
    total_summary <- summary(seurat.reference@active.ident)
    total_sum <- sum(total_summary)
    percentages <- summary(seurat.reference@active.ident)/
      sum(summary(seurat.reference@active.ident))

    group.size <- total_summary[group]
    group.percentage <- group.size/total_sum

    # If the main group is more that 70% or the parameter, I calculate the
    # size of the sampling, if not, it will be infinite
    if(group.percentage >= max.size.percentage){
      max.cells.per.ident <- total_sum * max.size.percentage
    }else{
      max.cells.per.ident <- Inf
    }
  }

  if(verbose){
    message('Seurat::FindMarkers -> ',
            'ident.1 = group: ', group,
            ' ident.2 = rest_groups: ', paste(shQuote(rest_groups), collapse=", "))
  }

  # Calling Seurat function
  markers.info.all <-  NULL
  tryCatch(
    {
      markers.info.all <- Seurat::FindMarkers(seurat.reference,
                                              ident.1 = group, ident.2 = rest_groups,
                                              test.use = test.use.value,
                                              max.cells.per.ident = max.cells.per.ident,
                                              only.pos=TRUE,
                                              logfc.threshold = param.logfc.threshold,
                                              verbose = verbose)

      message('Group: ', group, ' vs rest_groups = ', length(markers.info.all))

    },
    error = function(e){
      warning(paste0('The FindMarker function is failing in the group: ',
                     group, '. The error is: ', e))

      markers.info.all <- data.frame()
    }
  )

  if(marker_strategy == 'convert_inf_to_zero'){
    message('Executing convert_inf_to_zero strategy')
    markers.info.all$p_val[is.na(markers.info.all$p_val)] = 0
    markers.info.all$p_val_adj[is.na(markers.info.all$p_val_adj)] = 0
  }else if(marker_strategy == 'keep_default_values'){
    # nothing to do, therefore if the p-value and p-adj-value are NA, the
    # library will skip the markers.
  }

  return(markers.info.all)
}

#' Calculate Minimum Correlation with Incremental Markers
#'
#' This function processes a list of marker genes for each group, sorted in
#' descending order by log2 fold change. The objective is to determine the
#' number of genes that must be selected to achieve the minimum correlation
#' between samples in the reference dataset. The function is designed to
#' incrementally add markers, starting with those that display the greatest
#' change (as indicated by the log2 fold change). For each new set of markers,
#' the function calculates the correlation between samples in the reference
#' dataset.
#'
#' @param list_markers A list of marker genes for each group, ordered by log2
#'  fold change in descending order.
#' @param reference The reference dataset to which the correlations will be
#'  calculated.
#' @param minimum_markers The minimum number of marker genes to be selected,
#'  default is set to 4. This can be adjusted based on the analysis
#'  requirements.
#' @param min_delta_cor_threshold Minimum difference between the correlation of
#'  t and t-1 steps. The idea is to skip genes that have small changes in
#'  the correlation. For default it is a 5%.
#' @param verbose If TRUE, prints detailed information during the function
#'  execution. Default is TRUE.
#' @param verbose_detailed If TRUE, prints even more detailed information
#'  during the function execution. Default is TRUE.
#'
#' @details
#' The function begins by including the most significant markers (those with
#' the highest log2 fold change) in the analysis and calculates the correlation
#' between samples in the reference dataset. This process is repeated, each
#' time adding more markers, until a satisfactory level of correlation is
#' reached. This way, the function attempts to balance the inclusion of
#' informative markers against the potential noise introduced by less
#' significant markers.
#'
#' @return A data frame with the minimum number of markers required to achieve
#' a satisfactory level of correlation.
#'
#'
#' @export
calculate_min_correlation_incremental_markers <- function(
    list_markers,
    reference,
    minimum_markers = 4,
    min_delta_cor_threshold = 0.05,
    percentile_markers.min_corr = NULL,
    verbose = TRUE,
    verbose_detailed =  TRUE){

  message('min_delta_cor_threshold: ', min_delta_cor_threshold)
  message('minimum_markers: ', minimum_markers)

  if(verbose){
    message('Begining of calculate_min_correlation_incremental_markers:')
    print(list_markers)
  }

  # 1. Limits by each sample
  list_limits_by_sample = list()

  # iteration on the list markers to calculate the maximum number of markers
  max.markers <- 0
  for(marker_set in list_markers){
    number_markers <- nrow(marker_set)
    if(number_markers > max.markers){
      max.markers <- number_markers
    }
  }

  # 2. Best value in terms of correlation with the reference. 1 is the worse
  # case scenario when everything is 100% correlated
  min_value <- 1
  optimum_n_marker = -1
  accumulative_corr_markers <- NULL

  # 3. I'm adding a marker each time and comparing the correlation between
  # samples
  current_markers <- NULL
  for(marker_count in 1:max.markers) {

    #for each group a add a marker and check
    for(couter_market_set in 1:length(list_markers)) {

      marker_set <- list_markers[[couter_market_set]]

      # let's initialize the markers with the minimum amount
      if(marker_count <= minimum_markers){
        # To make sure we don't try to get values in a empty position.
        if(marker_count <= nrow(marker_set)){
          current_markers <- c(current_markers,
                               rownames(marker_set)[marker_count])

          joint_cor <- 0
          # Correlation has to be done with more than 1 genes.
          if(length(current_markers)>=2){
            # Calculating the joint correlation even if I added the gene anyways
            joint_cor <- calculate_joint_cor(reference = reference,
                                             markers = current_markers)
          }

          # Logging even when is based on the minimum_markers
          if(verbose_detailed){
            print(paste0(joint_cor, ', marker_count:',
                         marker_count, ', ',
                         names(list_markers)[couter_market_set],
                         ', number_markers_added: ', length(current_markers),
                         '. Gene added because: minimum_markers parameter.'))
          }

          #Add line to the accumulative_corr_markers
          accumulative_corr_markers.temp <- data.frame('joint_cor'=joint_cor,
                                                       'marker_count'=marker_count,
                                                       'gene'=rownames(marker_set)[marker_count],
                                                       'celltype'=names(list_markers)[couter_market_set],
                                                       'reason_added'='minimum_markers_parameter',
                                                       'current_marker'=length(current_markers))

        }
      }else{
        # To make sure we don't try to get values in a empty position.
        if(marker_count <= nrow(marker_set)){
          current_markers <- c(current_markers,
                               rownames(marker_set)[marker_count])
          joint_cor <- calculate_joint_cor(reference = reference,
                                           markers = current_markers)

          # I can get a negative correlation, but it will be the minimum one.
          joint_cor <- abs(joint_cor)

          # I calculate the delta correlation t-1
          delta_correlation <- abs(min_value - joint_cor)

          reason_added <- 'no_added'

          # Getting the minimum value and the length of the marker vector when
          # that happen
          # Also I need to check that this value is greater than zero and
          # if the current delta correlation between the t and t-1 is greater
          # or equal than the min_delta_cor_threshold I can add the value.
          if(min_value > joint_cor & joint_cor >= 0 &
             delta_correlation >= min_delta_cor_threshold){
            min_value <- joint_cor
            optimum_n_marker <- length(current_markers)
            message('New optimum added with delta_correlation=',
                    delta_correlation, ', min_delta_cor_threshold: ',
                    min_delta_cor_threshold)

            reason_added <- paste0('added_optimum_delta_corr=',delta_correlation)
          }

          if(verbose_detailed){
            print(paste0(joint_cor, ', marker_count:',
                         marker_count, ', ',
                         names(list_markers)[couter_market_set],
                         ', number_markers_added: ', length(current_markers)))
          }

          #Add line to the accumulative_corr_markers
          accumulative_corr_markers.temp <- data.frame('joint_cor'=joint_cor,
                                                       'marker_count'=marker_count,
                                                       'gene'=rownames(marker_set)[marker_count],
                                                       'celltype'=names(list_markers)[couter_market_set],
                                                       'reason_added'=reason_added,
                                                       'current_marker'=length(current_markers))

        }
      }

      #Adding the new row
      accumulative_corr_markers <- rbind(accumulative_corr_markers,
                                         accumulative_corr_markers.temp)

    }
  }

  #Lets take the percentile 25% by default or the parameter value.
  if(!is.null(percentile_markers.min_corr)){

    if(percentile_markers.min_corr > 1){
      error(paste0('percentile_markers.min_corr must be (0,1]. Current value: ', percentile_markers.min_corr))
    }else{
      print(paste0("Percentile filtering in min correlation function: ", percentile_markers.min_corr))
      quartile.calculation <- quantile((accumulative_corr_markers$joint_cor), c(percentile_markers.min_corr), na.rm = TRUE)
      accumulative_corr_markers.filtered <- accumulative_corr_markers[accumulative_corr_markers$joint_cor<=quartile.calculation[1],]
      optimum_n_marker <- accumulative_corr_markers.filtered[nrow(accumulative_corr_markers.filtered),]$current_marker
      min_value <- accumulative_corr_markers.filtered[nrow(accumulative_corr_markers.filtered),]$joint_cor
      print(paste0('Opmimum n marker: ', optimum_n_marker))
    }
  }

  # In the case that the available markers is <<< than the minimum marker number
  # In that case I used all the markers.
  if(optimum_n_marker==-1 ){

    # In case that the number of markers is >> than the minimum markers.
    if(length(current_markers) > minimum_markers){
      optimum_n_marker <- minimum_markers
      message('There are not optimum marker since the minimum number of markers ',
              'is: ', minimum_markers, ' which is << than the size of the  ',
              'calculated markers: ', length(current_markers))
    }else{
      optimum_n_marker <- length(current_markers)
      message('There are not optimum marker since the minimum number of markers ',
              'is: ', minimum_markers, ' which is > than the size of the  ',
              'calculated markers: ', length(current_markers))
    }
  }

  # Check the final results of the analysis
  if(verbose){
    print(paste0('Minimum value of average correlation: ',
                 min_value, ' - optimum_n_marker: ', optimum_n_marker))
  }

  # 3. I have to return the optimum markers
  return(list(markers=current_markers[1:optimum_n_marker],
              current_markers = current_markers,
              accumulative_corr_markers=accumulative_corr_markers))
}


#' Calculate Joint Correlation
#'
#' This function computes the correlation among distinct groups within the
#' reference data based on a set of marker genes. The function then calculates
#' the average joint correlation, which reflects the average correlation
#' coefficient across all possible pairs of different groups.
#'
#' @param reference A matrix or dataframe of reference data, with each
#'  column representing a sample and each row corresponding to a gene.
#' @param markers A vector of marker genes, which will be used to subset the
#'  reference dataset for the correlation computation.
#' @param method a character string indicating which correlation coefficient
#'  (or covariance) is to be computed. One of "pearson" (default), "kendall",
#'  or "spearman": can be abbreviated.
#'
#'  @import stats
#'
#' @details
#' The function works by initially calculating the correlation matrix of the
#' reference data subsetted by the specified marker genes.
#' It then computes the joint correlation by summing all correlation
#' coefficients, excluding the self-correlations (i.e., correlation of a group
#' with itself).
#' This sum is divided by the number of unique pair combinations of different
#' groups, yielding the average joint correlation.
#' The result is a single numeric value representing the average degree of
#' correlation between the different groups in the reference dataset based on
#' the specified marker genes.
#'
#' @return A single numeric value representing the average joint correlation.
#'
#' @import stats
#'
#'
#' @export
calculate_joint_cor <- function(reference, markers, method ="pearson"){
  cor_matrix <- stats::cor(reference[markers,],
                           use = "everything",
                           method = method)
  n_col <- ncol(reference)
  joint_cor <- (sum(cor_matrix)-sum(diag(cor_matrix)))/(n_col*n_col-n_col)
  return(joint_cor)
}


#' Calculate and Print Performance Metrics
#'
#' This function computes various performance metrics to assess the accuracy of
#' calculated cell type proportions relative to the true proportions. The
#' performance metrics are computed both by clusters and by samples.
#' The results can be optionally printed.
#'
#' @param true_prop A matrix or dataframe of the true proportions. Rows
#'  correspond to cell types and columns to samples.
#' @param calc_prop A matrix or dataframe of the calculated proportions. It
#'  should have the same structure as true_prop.
#' @param title Optional title for the print output.
#' @param verbose If TRUE, the performance metrics will be printed.
#'  Default is FALSE.
#'
#' @details
#' The function computes performance metrics, such as mean absolute error,
#' root mean squared error, and correlation coefficient, using the
#' 'calculate_performance_metrics' function. These metrics are calculated for
#' each cell type (i.e., by clusters) and for each sample.
#' If verbose is set to TRUE, these metrics will be printed with the
#' 'print_performance_metrics' function.
#'
#' @return A list with three elements: 'pm_clusters' contains the performance
#' metrics by clusters, 'pm_samples' contains the performance metrics by
#' samples, true_prop' contains the true proportions, and 'calc_prop' contains
#' the calculated proportions.
#'
#'
#' @export
performance_metrics_general <- function(true_prop, calc_prop,
                                        title = NULL, verbose = FALSE){

  # 1. Calculate metrics.
  pm_clusters <- calculate_performance_metrics(true_prop = true_prop,
                                               calc_prop = calc_prop)
  pm_samples <- calculate_performance_metrics(true_prop = t(true_prop),
                                              calc_prop = t(calc_prop))

  # 2. Show the metrics.
  if(verbose){
    print_performance_metrics(list_metrics = pm_clusters,
                              title = paste0('Metrics by clusters - ', title))
    print('')
    print_performance_metrics(list_metrics = pm_samples,
                              paste0(title = 'Metrics by samples - ', title))
  }

  # 3. Return metrics.
  return(list(pm_clusters = pm_clusters, pm_samples = pm_samples,
              true_prop = true_prop, calc_prop = calc_prop))
}


#' Compute Performance Metrics
#'
#' This function calculates several performance metrics to assess the quality
#' of deconvolution results. These metrics quantify the accuracy of the
#' estimated cell type proportions compared to the true proportions.
#'
#' @param true_prop A matrix or dataframe of the true proportions. Rows
#'  correspond to cell types and columns to samples.
#' @param calc_prop A matrix or dataframe of the estimated proportions. It
#'  should have the same structure as true_prop.
#' @param digits The number of decimal places to which the computed metrics
#'  should be rounded. Default is 20.
#'
#' @details
#' This function computes performance metrics including
#' Mean Absolute Error (MAE), Root Mean Squared Error (RMSE),
#' and Pearson Correlation Coefficient (r) to assess the accuracy of
#' deconvolution results. The metrics are calculated separately for each cell
#' type and for each sample. The precision of the returned metrics can be
#' adjusted with the 'digits' parameter.
#'
#' @return A list of calculated performance metrics. Each element of the list
#' represents a different metric.
#'
#' @import stats
#'
#' @export
calculate_performance_metrics <- function(true_prop, calc_prop, digits = 20){

  # I order the columns and rows equally.
  true_prop <- true_prop[rownames(calc_prop), colnames(calc_prop)]
  true_prop <- as.data.frame(true_prop)

  # In case that we have just one sample I have to fix the dimensions.
  # This happens when it is just a sample and then the dataframe is converted
  # in an array, then the columns are rows and the other way around.
  if(nrow(calc_prop) != nrow(true_prop)){
    true_prop <- t(true_prop)
  }

  rownames(true_prop) <- rownames(calc_prop)
  colnames(true_prop) <- colnames(calc_prop)

  # let's calculate some metrics: RMSD (total, by sample), mAD,
  # Pearson (total, by sample)

  # 1. RMSD: Root Mean Square Deviation (RMSD) is the most commonly used
  # quantitative measure of the similarity between two superimposed atomic
  # coordinates.
  RMSD <- round(sqrt(mean(as.matrix((true_prop - calc_prop)^2), na.rm = T)),
                digits = digits)
  RMSD_by_sample <- round(sqrt(rowMeans((true_prop - calc_prop)^2)),
                          digits = digits)

  # 2. mAD: median absolute deviation (MAD) is a robust measure of the
  # variability of a univariate sample of quantitative data. It can also refer
  # to the population parameter that is estimated by the MAD calculated from a
  # sample.
  mAD <- round(mean(as.matrix(abs(true_prop - calc_prop)), na.rm = T),
               digits = digits)
  mAD_by_sample <- round(rowMeans(abs(true_prop - calc_prop)),
                         digits = digits)

  # 3. Pearson: The Pearson correlation measures the strength of the linear
  # relationship between two variables. It has a value between -1 to 1, with a
  # value of -1 meaning a total negative linear correlation, 0 being no
  # correlation, and + 1 meaning a total positive correlation.
  Pearson <- round(stats::cor(c(as.matrix(true_prop)), c(as.matrix(calc_prop))),
                   digits = digits)
  Pearson_by_sample <- sapply(1:nrow(true_prop), function(x) {
    round(stats::cor(c(as.matrix(true_prop[x, ])), c(as.matrix(calc_prop[x, ]))),
          digits = digits)
  })
  names(Pearson_by_sample) <- rownames(true_prop)

  # 4. RSS: The residual sum of squares (RSS) is a statistical technique used
  # to measure the amount of variance in a data set that is not explained by a
  # regression model itself. Instead, it estimates the variance in the
  # residuals, or error term.
  RSS_residuals <- round(sum((true_prop - calc_prop)^2), digits = digits)

  return(list(RMSD = RMSD,
              RMSD_by_sample = RMSD_by_sample,
              mAD = mAD,
              mAD_by_sample = mAD_by_sample,
              Pearson = Pearson,
              Pearson_by_sample = Pearson_by_sample,
              RSS_residuals = RSS_residuals))
}


#' Display Performance Metrics
#'
#' This function neatly prints the calculated performance metrics, facilitating
#' the interpretation of deconvolution results.
#'
#' @param list_metrics A list of calculated performance metrics, as produced
#' by the `calculate_performance_metrics()` function.
#' @param title An optional title for the printed metrics. Default is NULL.
#'
#' @details
#' The function takes as input a list of performance metrics, which should
#' include metrics such as Mean Absolute Error (MAE),
#' Root Mean Squared Error (RMSE), and Pearson Correlation Coefficient (r).
#' These metrics are then printed in a structured format for easy
#' interpretation. The optional 'title' parameter can be used to give the
#' printed metrics a descriptive title.
#'
#' @return NULL. The function is used for its side effect of printing the
#' metrics.
#'
#' @export
print_performance_metrics <- function(list_metrics, title = NULL){
  if(!is.null(title)){
    print(title)
  }

  # 1. Iteration on the methics and printing the results.
  for(counter in 1:length(list_metrics)) {
    print(paste0(names(list_metrics)[counter], ': ', list_metrics[counter]))
  }
}


#' Recalculate Cell Type Proportions
#'
#' This function recalculates the proportions of cell types based on the
#' top-level clustering proportions and the more granular sub-cluster
#' proportions.
#'
#' @param top_cluster_deco_list A list containing the deconvolution results for
#'  the top-level clusters.
#' @param top_clusters_list A vector containing the names or identifiers of the
#'  top-level clusters.
#' @param subcluster_deco_list A list of deconvolution results for the
#'  sub-clusters.
#' @param verbose A logical indicating whether or not to print detailed messages
#'  while the function is running. Default is FALSE.
#'
#' @details
#' This function takes as input the top-level cluster proportions and a list of
#' sub-cluster deconvolution results. It scales the sub-cluster proportions
#' by the corresponding top-level cluster proportions and then combines all the
#' proportions. This results in an updated set of cell type proportions
#' that reflects the hierarchical structure of the clustering.
#'
#' @return A matrix where each row represents a cell type and each column
#' represents a sample. Each entry in the matrix represents the recalculated
#' proportion of a cell type in a sample.
#'
#'
#' @export
recalculate_proportions <- function(top_cluster_deco_list, top_clusters_list,
                                    subcluster_deco_list, verbose = FALSE){

  top_clusters_prop <- top_cluster_deco_list$h[top_clusters_list, ]
  # iterate in the top clustering
  calculated_proportions <- NULL
  for (top_cluster_counter in 1:length(top_clusters_list)) {
    top_cluster_name = top_clusters_list[top_cluster_counter]
    subcluster_deco <- subcluster_deco_list[[top_cluster_name]]
    relative_prop <- subcluster_deco$h

    rescaled_proportions <- sweep(relative_prop,
                                  MARGIN=2,
                                  t(top_clusters_prop[top_cluster_name,]), `*`)

    # Integration of the rescaled proportions
    if(is.null(calculated_proportions)){
      calculated_proportions <- rescaled_proportions
    }else{
      calculated_proportions <- rbind(calculated_proportions,
                                      rescaled_proportions)
    }
  }

  calculated_proportions
}
