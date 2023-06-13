############### DecoFlex Package: Automatic Tree Discovery Module ##############
##
## This file, automatic_tree_discovery.R, is a crucial component of the DecoFlex
## R package. It encapsulates a collection of specialized functions designed for
## the automatic identification and optimization of hierarchical levels and
## clusters.
##
## The functions implemented within this module empower users to streamline
## their deconvolution analysis by automating tree discovery and related
## processes. Key functions included in this file are:
##
## 1. hierachy_selection_plus_deconvolution: Performs clustering method and
##    deconvolution.
## 2. calculate_best_cluster_global: Determines the optimal number of clusters
##    to minimize correlation.
## 3. filter_main_celltypes: Filters primary cell types using marker genes.
## 4. calculate_best_clustering_generic: Analyzes a level of clustering,
##    minimizing correlation.
## 5. calculate_k_clustering: Calculates the optimal k clustering cut for
##    hierarchical clustering.
## 6. plot_hierarchical_clustering: Visualizes the hierarchical clustering and
##    correlations.
## 7. plot_heatmap_correlation_reference: Generates a heatmap showcasing the
##    correlation among different cell-types.
## 8. add_clustering_single_cell_dataset: Enhances a single cell dataset by
##    appending metacluster information.
##
## Each function contributes to the seamless execution of automated tree
## discovery, ensuring that the most relevant, insightful hierarchical
## structures are generated for subsequent deconvolution. By utilizing these
## functions, users can harness the power of DecoFlex's deconvolution
## capabilities with minimal manual intervention.
##


#' Hierarchy Selection Plus Deconvolution
#'
#' Runs the clustering method and the deconvolution process on given single-cell
#' expression data.
#'
#' @name hierachy_selection_plus_deconvolution
#'
#' @param single_cell_data_exp A data frame or matrix containing single-cell
#'  expression data.
#' @param var_cell_type A vector or factor indicating the cell type for each
#'  cell.
#' @param subset_celltypes A vector of cell types to subset for the analysis.
#'  Default is NULL, implying all cell types are used.
#' @param var_sample A vector or factor indicating the sample source for each
#'  cell.
#' @param number_clusters_one_celltype An integer specifying the number of
#'  clusters for each cell type. Default is 1.
#' @param distance_method A string specifying the method to compute the
#'  distance. Default is 'euclidean'.
#' @param hclust_method A string specifying the method to perform hierarchical
#'  clustering. Default is 'average'.
#' @param min_size_leaf An integer specifying the minimum size of the leaf in
#'  the hierarchical clustering. Default is 3.
#' @param max_k_tried_hier_clustering An integer specifying the maximum number
#'  of clusters tried in hierarchical clustering. Default is 3.
#' @param random_seed An integer specifying the random seed for reproducibility.
#'  Default is NULL.
#' @param use_min_cor_strategy A logical value indicating whether to use the
#'  minimum correlation strategy. Default is TRUE.
#' @param delete_shared_level_markers A logical value indicating whether to
#'  delete markers shared at level of the hierarchy. Default is FALSE.
#' @param delete_shared_internal_markers A logical value indicating whether to
#'  delete markers shared internally. Default is FALSE.
#' @param filter_markers A vector specifying which markers to filter. Default
#' is NULL.
#' @param param.logfc.threshold A numeric value indicating the log-fold-change
#'  threshold for marker selection. Default is 2.0.
#' @param param.p_val_adj A numeric value indicating the adjusted p-value
#'  threshold for marker selection. Default is 0.05.
#' @param filter_main_markers A logical value indicating whether to filter
#'  main markers. Default is TRUE.
#' @param verbose A logical value indicating whether to print additional output
#'  during the computation. Default is FALSE.
#'
#' @return A list containing results of the clustering and deconvolution
#'  process.
#'
#' @export
hierachy_selection_plus_deconvolution <- function(
    single_cell_data_exp, var_cell_type,
    subset_celltypes = NULL, var_sample,
    number_clusters_one_celltype = 1,
    distance_method = 'euclidean',
    hclust_method = "average",
    min_size_leaf = 3,
    max_k_tried_hier_clustering = 3,
    random_seed = NULL,
    use_min_cor_strategy = TRUE,
    delete_shared_level_markers = FALSE,
    delete_shared_internal_markers = FALSE,
    filter_markers = NULL,
    param.logfc.threshold = 2.0,
    param.p_val_adj = 0.05,
    filter_main_markers = TRUE,
    verbose = FALSE){

  # 1. Calculation of optimum hierarchical clustering.
  if(verbose){
    message('1. Starting of process to find the optimum hierarchical ',
            'clustering...')
  }

  result.optimum.hierarchy <- calculate_best_cluster_global(
    single_cell_data_exp = single_cell_data_exp,
    var_cell_type = var_cell_type,
    subset_celltypes = subset_celltypes,
    number_clusters_one_celltype = number_clusters_one_celltype,
    min_size_leaf = min_size_leaf,
    var_sample = var_sample,
    use_min_cor_strategy = use_min_cor_strategy,
    delete_shared_level_markers = delete_shared_level_markers,
    delete_shared_internal_markers = delete_shared_internal_markers,
    random_seed = random_seed,
    filter_main_markers =  filter_main_markers,
    param.logfc.threshold = param.logfc.threshold,
    param.p_val_adj = param.p_val_adj,
    verbose = verbose)

  # 2. Deconvolution process with the input of the process 1.
  if(verbose){
    message('2. Starting of process to deconvoluting the data...')
  }

  deco.simulation.results <- run_deconvolution_simulation_generic_recursive(
    single_cell_data_exp = single_cell_data_exp,
    hierarchy = result.optimum.hierarchy$tree_result,
    sub_clusters_var = var_cell_type,
    sample = var_sample,
    use_min_cor_strategy = use_min_cor_strategy,
    delete_shared_level_markers = delete_shared_level_markers,
    delete_shared_internal_markers = delete_shared_internal_markers,
    filter_markers = filter_markers,
    param.logfc.threshold = param.logfc.threshold,
    param.p_val_adj = param.p_val_adj,
    verbose = verbose)

  # 3. Return both: the optimum hierarchy and the deconvolution
  return(list(result.optimum.hierarchy = result.optimum.hierarchy,
              deco.simulation.results = deco.simulation.results))
}

#' Calculate the Best Global Cluster
#'
#' Creates the tree of cluster groups and levels. Decides the number of clusters
#' that will minimize the correlation between them, given a list of cell types.
#'
#' @name calculate_best_cluster_global
#'
#' @param single_cell_data_exp A data frame or matrix containing single-cell
#'  expression data.
#' @param var_cell_type A vector or factor indicating the cell type for each
#'  cell.
#' @param subset_celltypes A vector of cell types to subset for the analysis.
#'  Default is NULL, implying all cell types are used.
#' @param var_sample A vector or factor indicating the sample source for each
#'  cell.
#' @param use_min_cor_strategy A logical value indicating whether to use the
#'  minimum correlation strategy. Default is TRUE.
#' @param delete_shared_level_markers A logical value indicating whether to
#'  delete markers shared at level of the hierarchy. Default is FALSE.
#' @param delete_shared_internal_markers A logical value indicating whether to
#'  delete markers shared internally. Default is FALSE.
#' @param number_clusters_one_celltype An integer specifying the number of
#'  clusters for each cell type. Default is 1.
#' @param distance_method A string specifying the method to compute the
#'  distance. Default is 'euclidean'.
#' @param hclust_method A string specifying the method to perform hierarchical
#'  clustering. Default is 'average'.
#' @param min_size_leaf An integer specifying the minimum size of the leaf in
#'  the hierarchical clustering. Default is 3.
#' @param max_k_tried_hier_clustering An integer specifying the maximum number
#'  of clusters tried in hierarchical clustering. Default is 3.
#' @param random_seed An integer specifying the random seed for reproducibility.
#'  Default is NULL.
#' @param param.logfc.threshold A numeric value indicating the log-fold-change
#'  threshold for marker selection. Default is 2.0.
#' @param param.p_val_adj A numeric value indicating the adjusted p-value
#'  threshold for marker selection. Default is 0.05.
#' @param filter_main_markers A logical value indicating whether to filter main
#'  markers. Default is TRUE.
#' @param verbose A logical value indicating whether to print additional output
#'  during the computation. Default is FALSE.
#'
#' @return A list containing results of the clustering process.
#'
#' @export
calculate_best_cluster_global <- function(
    single_cell_data_exp,
    var_cell_type,
    subset_celltypes = NULL,
    var_sample,
    use_min_cor_strategy = TRUE,
    delete_shared_level_markers = FALSE,
    delete_shared_internal_markers = FALSE,
    number_clusters_one_celltype = 1,
    distance_method = 'euclidean',
    hclust_method = "average",
    min_size_leaf = 3,
    max_k_tried_hier_clustering = 3,
    random_seed = NULL,
    param.logfc.threshold = 2.0,
    param.p_val_adj = 0.05,
    filter_main_markers =  TRUE,
    verbose = FALSE){

  message("Starting tree-guided deconvolution...")

  # In order to reproduce results I can use a seed for random values like
  # the hierarchical clustering
  if(is.null(random_seed)){
    rm(.Random.seed, envir=globalenv())
  }else{
    # For testing: 240
    set.seed(random_seed)
  }

  # Just in case we want to run the generic analysis over a subset of cell
  # types, but still with a recursive method.
  subset_celltypes_filter = NULL
  if(!is.null(subset_celltypes)){
    subset_celltypes_filter <- subset_celltypes
    # I filter the single cell to only have the subset of cell types that I
    # want to analyze.
    single_cell_data_exp <-
      single_cell_data_exp[,single_cell_data_exp[[var_cell_type]] %in%
                             subset_celltypes]
  }else{
    subset_celltypes_filter <- unique(single_cell_data_exp[[var_cell_type]])
  }

  # I create the reference of all celltypes with the filtering
  reference_celltypes_top <- decoflex_build_cell_reference(
    x = single_cell_data_exp,
    ct.sub = subset_celltypes_filter,
    ct.varname = var_cell_type,
    sample = var_sample, verbose = verbose)

  # I get all the markers for each cell and then I'm using those for the analysis
  reference_celltypes_top.filtered = NULL
  if(filter_main_markers){
    message('Filtering the main dataset with just the marker ',
            'genes of the whole set of celltypes.')
    reference_celltypes_top.filtered.object <- filter_main_celltypes(
      single_cell_data_exp = single_cell_data_exp,
      reference_celltypes_top = reference_celltypes_top,
      group_clusters_var = var_cell_type,
      use_min_cor_strategy = use_min_cor_strategy,
      delete_shared_level_markers = delete_shared_level_markers,
      delete_shared_internal_markers = delete_shared_internal_markers,
      param.logfc.threshold = param.logfc.threshold,
      param.p_val_adj = param.p_val_adj,
      verbose = verbose)

    reference_celltypes_top.filtered <-
      reference_celltypes_top.filtered.object$reference_celltypes_top.filtered


    # I filter the single cell to only have the subset of markers genes that I
    # want to analyze.
    single_cell_data_exp <- single_cell_data_exp[
      reference_celltypes_top.filtered.object$markers_subclusters.object$
        total_markers, ]


  }else{
    # In case that the parameter is FALSE I will use the dataset without
    # filtering it
    reference_celltypes_top.filtered <- reference_celltypes_top
  }

  # I send all cell types as first step with a generic function. The tree-level
  # will be zero and the leaf_number is 1, meaning the root
  results_top <- calculate_best_clustering_generic(
    single_cell_data_exp = single_cell_data_exp,
    reference = reference_celltypes_top.filtered$basis,
    var_cell_type = var_cell_type,
    cell_types = subset_celltypes_filter,
    number_clusters_one_celltype = number_clusters_one_celltype,
    distance_method = distance_method,
    hclust_method = hclust_method,
    min_size_leaf = min_size_leaf,
    max_k_tried_hier_clustering = max_k_tried_hier_clustering,
    var_sample = var_sample,
    verbose = verbose)

  return(list(reference_celltypes_top = reference_celltypes_top,
              tree_result = list('1'= results_top)))
}

#' Filter Main Cell Types with Markers
#'
#' Filters the main cell types using marker genes and performs clustering
#' optimization on the filtered data.
#'
#' @name filter_main_celltypes
#'
#' @param single_cell_data_exp A data frame or matrix containing single-cell
#'  expression data.
#' @param reference_celltypes_top A list of reference cell types to be used for
#'  marker identification.
#' @param group_clusters_var A string or variable indicating the clusters for
#'  each cell.
#' @param use_min_cor_strategy A logical value indicating whether to use the
#'  minimum correlation strategy. Default is TRUE.
#' @param delete_shared_level_markers A logical value indicating whether to
#'  delete markers shared at level of the hierarchy. Default is FALSE.
#' @param delete_shared_internal_markers A logical value indicating whether to
#'  delete markers shared internally. Default is FALSE.
#' @param param.logfc.threshold A numeric value indicating the log-fold-change
#'  threshold for marker selection. Default is 2.0.
#' @param param.p_val_adj A numeric value indicating the adjusted p-value
#'  threshold for marker selection. Default is 0.05.
#' @param verbose A logical value indicating whether to print additional output
#'  during the computation. Default is FALSE.
#'
#' @return A data frame or matrix containing the filtered single-cell data.
#'
#' @export
filter_main_celltypes <- function(single_cell_data_exp,
                                  reference_celltypes_top,
                                  group_clusters_var,
                                  use_min_cor_strategy = TRUE,
                                  delete_shared_level_markers = FALSE,
                                  delete_shared_internal_markers = FALSE,
                                  param.logfc.threshold = 2.0,
                                  param.p_val_adj = 0.05,
                                  verbose = FALSE){

  # Calculate the marker genes for the separated deconvolution not using the
  # marker strategy.
  markers_subclusters.object <- calculate_markers(
    single_cell_data_exp = single_cell_data_exp,
    reference = reference_celltypes_top$basis,
    group_clusters_var = group_clusters_var,
    use_min_cor_strategy = use_min_cor_strategy,
    delete_shared_level_markers = delete_shared_level_markers,
    delete_shared_internal_markers = delete_shared_internal_markers,
    param.logfc.threshold = param.logfc.threshold,
    param.p_val_adj = param.p_val_adj,
    verbose =  verbose)

  # Assign the final markers and the top marker list
  markers_subclusters <- markers_subclusters.object$total_markers
  list.subcluster.marker <- markers_subclusters.object$list_markers

  # Filtering the markers genes
  reference_celltypes_top.filtered <- reference_celltypes_top
  reference_celltypes_top.filtered$basis <-
    reference_celltypes_top.filtered$basis[markers_subclusters, ]

  if(verbose){

    # Plot of the correlation heatmap with just the marker genes.
    plot_heatmap_correlation_reference(
      reference_w.k.celltypes = reference_celltypes_top.filtered,
      title = paste0('Heatmap for the main dataset: all cell-types'))
  }

  return(list(
    reference_celltypes_top.filtered = reference_celltypes_top.filtered,
    markers_subclusters.object = markers_subclusters.object))
}

#' Calculate Optimal Clustering at Generic Level
#'
#' This function recursively analyses a level of clustering by attempting to
#' separate it into 2, 3, or 4 clusters. The aim is to minimize the correlation
#' between the groups. For each subgroup, the function is called again until
#' each cluster contains 2 or 3 cell types. This enables granular and optimized
#' clustering across different levels of hierarchy.
#'
#' @name calculate_best_clustering_generic
#'
#' @param single_cell_data_exp A data frame or matrix containing single-cell
#'  expression data.
#' @param reference A reference data set or object to guide the clustering
#'  process.
#' @param var_cell_type A string or variable indicating the cell type for each
#'  cell.
#' @param cell_types A vector of cell types present in the data.
#' @param var_sample A string or variable indicating the sample for each cell.
#' @param tree_level An integer indicating the current level of the tree
#'  hierarchy in the clustering process. Default is 0.
#' @param leaf_number An integer indicating the number of leaf nodes in the tree
#'  hierarchy. Default is 1.
#' @param number_clusters_one_celltype An integer indicating the number of
#'  clusters per cell type. Default is 1.
#' @param distance_method A string indicating the distance measure to be used
#'  for clustering. Default is 'euclidean'.
#' @param hclust_method A string indicating the method to be used for
#'  hierarchical clustering. Default is 'average'.
#' @param min_size_leaf An integer indicating the minimum size of the leaf
#'  nodes. Default is 3.
#' @param max_k_tried_hier_clustering An integer indicating the maximum number
#'  of clusters tried in the hierarchical clustering. Default is 3.
#' @param verbose A logical value indicating whether to print additional output
#'  during the computation. Default is FALSE.
#'
#' @return A list containing the optimal clustering results.
#'
#' @export
calculate_best_clustering_generic <- function(
    single_cell_data_exp, reference, var_cell_type, cell_types, var_sample,
    tree_level = 0, leaf_number = 1, number_clusters_one_celltype = 1,
    distance_method = 'euclidean', hclust_method = "average", min_size_leaf = 3,
    max_k_tried_hier_clustering = 3, verbose = FALSE){

  optimum_clustering <- NULL
  minimum_value_joint_cor <- 1.0

  if(verbose){
    message("Level: ", tree_level, ' leaf: ', leaf_number)
  }

  # Iteration between 2 and 3, since we want low rank problems. The best case
  # scenario is to have k=2.
  k_result_list = list()
  for (counter in 2:max_k_tried_hier_clustering) {

    # Checking if the number of cells is compatible with the k analysis.
    # Therefore, it is not going to try with k=4 when there are 5 celltypes,
    # since it will raise a error (k must be between 2 and celltypes-1)
    if(counter <= (length(cell_types)-1)){
      k_result <- calculate_k_clustering(
        single_cell_data_exp = single_cell_data_exp,
        reference = reference,
        var_cell_type = var_cell_type,
        cell_types = cell_types,
        k_value = counter,
        var_sample = var_sample,
        extra_title = paste0(' level_', tree_level, '_leaf_', leaf_number),
        distance_method = distance_method, hclust_method = hclust_method,
        verbose = verbose)


      # We don't want to have more than 1 or zero cluster with just one cell
      # type. There is a parameter by default in 1.
      k_fit_summary <- table(k_result$fit)
      # Clusters with more than 1 celltypes
      k_fit_summary_summary <- table(k_fit_summary > 1)

      # If the option TRUE exist, which means that there are clusters with 1
      # cell type, if it doesn't exist there are none cluster with 1 cell type.
      add_registry <-  FALSE
      if('FALSE' %in% names(k_fit_summary_summary))
      {
        # Now let's see how many 1-celltype cluster we have. I recommend to
        # have zero or 1 group with one cell type.
        if(k_fit_summary_summary[rownames(k_fit_summary_summary) %in%
                                 c(FALSE)] <= number_clusters_one_celltype){
          add_registry <- TRUE
        }else{
          add_registry <-  FALSE
        }
      }else{
        # All the clusters have more than one cell-type.
        add_registry <- TRUE
      }

      # Now we can add or not based on the rules applied already
      if(add_registry){

        # Add the model to the list of possible models.
        k_result_list[[counter]] <- k_result

        # Changing the minimum value of the average of correlation between
        # groups.
        if(minimum_value_joint_cor > k_result$value_joint_cor){
          minimum_value_joint_cor <- k_result$value_joint_cor
          optimum_clustering <- k_result
        }
      }else{
        print(paste0('The model is not compatible with the policy:',
                     ' level_', tree_level, '_leaf_', leaf_number))
      }
    }
  }

  # If after the analysis I can't find a good clustering, that will means that
  # the policy used is not compatible with the reference given.
  if(is.null(optimum_clustering)){
    stop('Based on the parameters, it is impossible to find cell-type ',
         'clustering at a more detailed level. Check with a different',
         ' value of number_clusters_one_celltype. Right now is: ',
         number_clusters_one_celltype)
  }

  # I assing the cell-type names to the list to complete the informative
  # structure of the method
  celltype_list <- names(optimum_clustering$fit)

  if(verbose){
    message('Cell-types for the level ', tree_level,
            ' (', length(celltype_list), '): ',
            paste0(celltype_list, collapse=', '))
    print(paste0('Optimum k: ', optimum_clustering$k_value,
                 ', with minimum correlation value: ',
                 optimum_clustering$value_joint_cor))

    # Heatmap with the optimum reference with the selected k.
    plot_hierarchical_clustering(
      hierar_clust = optimum_clustering$hierar_clust,
      k_value = optimum_clustering$k_value,
      fit_celltypes = optimum_clustering$fit,
      title = paste0('Hierarchical clustering with optimum level_',
                     tree_level, '_leaf_', leaf_number))
    # Plot of the optimum reference with the k clustering selected.
    plot_heatmap_correlation_reference(
      reference_w.k.celltypes = optimum_clustering$reference_w.k.celltypes,
      title = paste0('Heatmap with optimum level_',
                     tree_level, '_leaf_', leaf_number))
  }

  # I add the level models to a list. For each submodel I will save the level
  # and the results for the level with the list of subclusters.
  next_level_clustering <- list()

  # adding the iterative/recursive part
  # I analyze each of the leaves of the clustering
  freq.clusters <- table(optimum_clustering$fit)
  for (leaf_counter in rownames(freq.clusters)) {

    # I get the list of cell types that compose the analized leaf
    list_celltypes_leaf <-
      names(optimum_clustering$fit[optimum_clustering$fit==leaf_counter])

    # If is greater than 3 I will analyze deeply, if not, the process just
    # stops
    if(freq.clusters[leaf_counter] >= min_size_leaf){
      print(paste0('Creation of additional level with celltypes: ',
                   length(list_celltypes_leaf)))

      # I create the reference of all celltypes with the filtering
      reference_celltypes_leaf <- decoflex_build_cell_reference(
        x = single_cell_data_exp,
        ct.sub = list_celltypes_leaf,
        ct.varname = var_cell_type,
        sample = var_sample)

      # Let's call the generic function recursively to calculate the optimum
      # k hierarchical clustering minimizing correlation between groups of
      # cell-types.
      results_leaf <- calculate_best_clustering_generic(
        single_cell_data_exp = single_cell_data_exp,
        reference = reference_celltypes_leaf$basis,
        var_cell_type = var_cell_type, #TODO: change a new var with the level.
        cell_types = list_celltypes_leaf,
        tree_level = (tree_level + 1),
        leaf_number = leaf_counter,
        number_clusters_one_celltype = number_clusters_one_celltype,
        min_size_leaf = min_size_leaf,
        var_sample = var_sample)
      # Adding the leaf result to the list.
      next_level_clustering[[leaf_counter]] <- results_leaf

      if(verbose){
        message("Process ended for the next level for leaf part: ",
                leaf_counter, ' of the model level: ', tree_level)
      }
    }else{
      # I will show the cell-types just for the final leaf
      if(verbose){
        message("Final part of level ", tree_level, " leaf: ", leaf_counter)
        message('Final cell-types for level ', tree_level,
                " leaf: ", leaf_counter,
                ' (', length(list_celltypes_leaf), '): ',
                paste0(list_celltypes_leaf, collapse=', '))
      }

      # Even though it is the last part for the left, I will add to keep the
      # structure of the tree
      final_leaf_summary <- list(tree_level = (tree_level+1),
                                 leaf_number = leaf_counter,
                                 number_clusters_one_celltype = 0,
                                 celltype_list = list_celltypes_leaf,
                                 next_level_clustering = NULL)
      next_level_clustering[[leaf_counter]] <- final_leaf_summary
    }
  }

  return(list(tree_level = tree_level, leaf_number = leaf_number,
              number_clusters_one_celltype = number_clusters_one_celltype,
              celltype_list = celltype_list,
              next_level_clustering = next_level_clustering,
              optimum = optimum_clustering,
              k_result_list = k_result_list))
}


#' Optimal k Clustering for Hierarchical Clustering
#'
#' This function calculates the optimal number (k) of clusters for a
#' hierarchical clustering given a specific set of cell types. It can be
#' customized with different distance measures and agglomeration methods.
#'
#' @name calculate_k_clustering
#'
#' @param single_cell_data_exp A data frame or matrix containing single-cell
#'  expression data.
#' @param reference A reference data set or object to guide the clustering
#'  process.
#' @param var_cell_type A string or variable indicating the cell type for each
#'  cell.
#' @param cell_types A vector of cell types present in the data.
#' @param k_value An integer indicating the number of clusters to be created.
#' @param var_sample A string or variable indicating the sample for each cell.
#' @param distance_method A string indicating the distance measure to be used
#'  for clustering. This must be one of "euclidean", "maximum", "manhattan",
#'  "canberra", "binary" or "minkowski". Any unambiguous substring can be given.
#'  Default is 'euclidean'.
#' @param hclust_method A string indicating the method to be used for
#'  hierarchical clustering. This should be (an unambiguous abbreviation of)
#'  one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA),
#'  "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC). Default
#'  is 'average'.
#' @param extra_title A string used for the additional title in the plotting.
#' @param verbose A logical value indicating whether to print additional output
#'  during the computation. Default is TRUE.
#'
#' @import stats
#'
#' @return A list containing the optimal k clustering results.
#'
#' @export
calculate_k_clustering <- function(single_cell_data_exp, reference,
                                   var_cell_type, cell_types, k_value,
                                   var_sample,
                                   distance_method = 'euclidean',
                                   hclust_method = "average",
                                   extra_title = NULL, verbose = TRUE){
  # Finding distance matrix
  distance_mat <- stats::dist(t(stats::cor(reference)), method = distance_method)

  # Fitting Hierarchical clustering Model to training dataset
  # Setting seed,  "complete", "average" ,"ward.D"
  hierar_clust <- stats::hclust(distance_mat, method = hclust_method)

  # Cutting tree by no. of clusters and choosing no. of clusters
  fit_celltypes <- stats::cutree(hierar_clust, k = k_value)

  if(verbose){
    title = paste0('Dendrogram - ', extra_title, ' - distance_method: ',
                   distance_method, ', hclust_method: ',
                   hclust_method, ', k_value: ', k_value)

    # Let's call the plot for the hierarchical clustering
    plot_hierarchical_clustering(hierar_clust = hierar_clust,
                                 k_value = k_value,
                                 fit_celltypes = fit_celltypes,
                                 title = title)
  }

  single_cell_data_exp.metacluster <- add_clustering_single_cell_dataset(
    sc.eset = single_cell_data_exp,
    var_cluster_name = var_cell_type,
    fit_hierarchical = fit_celltypes)

  # I create the reference of all celltypes
  reference_w.k.celltypes <- decoflex_build_cell_reference(
    x = single_cell_data_exp.metacluster,
    ct.sub = NULL,
    ct.varname = 'metacluster',
    sample = var_sample,
    verbose = verbose)

  value_joint_cor <- calculate_joint_cor(
    reference = reference_w.k.celltypes$basis,
    markers = rownames(reference_w.k.celltypes$basis))

  # Return the fit model
  return(list(fit = fit_celltypes, hierar_clust = hierar_clust,
              reference_w.k.celltypes = reference_w.k.celltypes,
              value_joint_cor = value_joint_cor, k_value = k_value))
}

#' Plot hierarchical clustering
#'
#' This function generates a plot of the hierarchical clustering and includes
#' the option for a custom title. The plot visualizes the dendrogram resulting
#' from the clustering, with the height of the cut denoted by a green line.
#' Clusters resulting from the cut are highlighted with green borders.
#' The function also outputs a table summarizing the resulting cluster
#' assignments.
#'
#' @name plot_hierarchical_clustering
#' @param hierar_clust A hierarchical clustering object, which is the result of
#'  hclust function or similar.
#' @param k_value An integer value indicating the number of clusters to
#'  highlight in the plot.
#' @param fit_celltypes A vector containing the assignment of each data point
#'  to a cluster.
#' @param title Optional; a character string representing the title of the plot.
#'  If NULL, no title is added.
#' @return This function generates a plot as a side-effect. The returned value
#'  is a table summarizing the cluster assignments.
#'
#' @import graphics
#' @import stats
#'
#' @export
plot_hierarchical_clustering <- function(hierar_clust, k_value,
                                         fit_celltypes, title = NULL){
  # This will affect also legend title font size
  graphics::par(cex.main=0.5)

  # Plotting dendrogram
  plot(hierar_clust, main = title)

  # Cutting tree by height
  graphics::abline(h = 110, col = "green")

  # green line with the cluster selected
  stats::rect.hclust(hierar_clust, k = k_value, border = "green")

  table(fit_celltypes)
}


#' Plot Heatmap of Correlation Among Cell-Types or Groups
#'
#' This function creates a heatmap to visually display the correlation between
#' different cell-types or groups of cell-types. The color gradient of the
#' heatmap represents the correlation coefficient values, providing a clear
#' visual summary of the correlations in the data.
#'
#' @name plot_heatmap_correlation_reference
#' @param reference_w.k.celltypes A data frame containing the cell-types or
#' groups of cell-types and their corresponding basis values.
#' @param title Optional; a character string representing the title of the
#'  heatmap. If not provided, the title defaults to 'Heatmap'.
#' @return This function generates a heatmap as a side-effect. It does not
#'  return any value.
#'
#'  @import grDevices
#'  @import ggplot2
#'  @import gplots
#'  @import graphics
#'  @import stats
#'
#' @export
plot_heatmap_correlation_reference <- function(reference_w.k.celltypes,
                                               title = 'Heatmap'){

  color_gradient <- grDevices::colorRampPalette(c("white","steelblue"))

  # Calculation of correlation between groups.
  cor.reference_matrix_w.k.celltypes <- stats::cor(reference_w.k.celltypes$basis)

  graphics::par(cex.main=0.5) ## this will affect also legend title font size

  # Plot of correlations between the groups or cell types
  gplots::heatmap.2(cor.reference_matrix_w.k.celltypes,
                    cellnote = round(cor.reference_matrix_w.k.celltypes, 2),
                    trace="none",
                    col=color_gradient(10),
                    breaks=seq(0,1,0.1),
                    margins=c(12,12),
                    cexRow = 1, cexCol = 1,
                    dendrogram = "both",
                    key=FALSE,
                    main=title)
}


#' Add Clustering Information to Single Cell Dataset
#'
#' This function enhances a single cell dataset by adding metacluster
#' information to it. It achieves this by iterating through each unique value
#' in the hierarchical fit, identifying the celltypes associated with that
#' value, and tagging the corresponding cells in the single cell dataset with a
#' new metacluster identifier.
#'
#' @name add_clustering_single_cell_dataset
#' @param sc.eset A single cell ExpressionSet object. This dataset is modified
#'  in-place to include metacluster identifiers.
#' @param var_cluster_name A string specifying the column name in sc.eset where
#'  the cluster identifiers are stored. Default is 'cluster_normalized'.
#' @param fit_hierarchical A named vector where the names correspond to
#'  celltypes and the values represent the hierarchical clustering fit for each
#'  celltype.
#' @return An ExpressionSet object (sc.eset) with additional metacluster
#'  information.
#'
#' @export
add_clustering_single_cell_dataset <- function(
    sc.eset,
    var_cluster_name='cluster_normalized',
    fit_hierarchical){

  list_values <- unique(fit_hierarchical)

  for (counter in list_values) {

    list_celltypes <- names(fit_hierarchical[fit_hierarchical==counter])
    sc.eset$metacluster[sc.eset[[var_cluster_name]] %in% list_celltypes] <-
      paste0("Subcluster_", counter)
  }

  return(sc.eset)
}
