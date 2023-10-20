library(DecoFlex)

run.decoflex <- function(single_cell_data_exp,  sub_clusters_var, hierarchy, bulk_data, max_iterations = 15000){
  
  #Parameters for the simulation
  single_cell_data_exp <- single_cell_data_exp
  hierarchy = hierarchy
  bulk_data = bulk_data
  sub_clusters_var <- sub_clusters_var
  sample <- 'Disorder' #just one value: Control
  use_min_cor_strategy = TRUE
  delete_shared_level_markers = FALSE
  delete_shared_internal_markers = FALSE
  ordering_strategy = 'foldchange_pvalue'
  
  start_time <- Sys.time()
  
  #let's run de simulation
  deco.actual.data <- 
    run_deconvolution_tree_guided_recursive(bulk_data = bulk_data,
                                            single_cell_data_exp = single_cell_data_exp, 
                                            hierarchy = hierarchy,
                                            sub_clusters_var = sub_clusters_var,
                                            sample = sample, 
                                            use_min_cor_strategy = use_min_cor_strategy,
                                            ordering_strategy = ordering_strategy,
                                            delete_shared_level_markers = delete_shared_level_markers, 
                                            delete_shared_internal_markers = delete_shared_internal_markers,
                                            param.logfc.threshold = 1.0,
                                            param.p_val_adj = 0.05,
                                            minimum_markers = 4,
                                            min_delta_cor_threshold = 0.0,
                                            max_iterations = max_iterations,
                                            verbose = TRUE)
  
  end_time <- Sys.time()
  end_time - start_time
  
  return(deco.actual.data)
}


#Function to delete a given celltype name in all the hierarchical structure
delete_celltype_hierarchy <- function(hierarchy, celltype_name_delete){
  current_hierarchy <- hierarchy
  
  #replace celltype name in the current level of the hierarchy
  current_hierarchy$celltype_list <- current_hierarchy$celltype_list[!current_hierarchy$celltype_list == celltype_name_delete]
  
  if(!is.null(current_hierarchy$next_level_clustering)){
    
    #for the next hierarchy or level
    current_hierarchy$next_level_clustering <- delete_celltype_hierarchy(hierarchy = current_hierarchy$next_level_clustering,
                                                                         celltype_name_delete = celltype_name_delete)
  }
  
  return(current_hierarchy)
}


library(DecoFlex)

run_deco_CA_orig.celltype <- function(data.correlation.all.CA.orig.celltype, marker_value, max_iterations=25000){
  
  markers.clusters.object <- list()
  markers.clusters.object$accumulative_corr_markers <- data.correlation.all.CA.orig.celltype
  markers.clusters.object$list_markers <- NULL
  markers.clusters.object$list_markers.percentile <- NULL
  
  #Adding just the markers that I'm interesting on.
  marker_filtered <- data.correlation.all.CA.orig.celltype[data.correlation.all.CA.orig.celltype$current_marker == marker_value,]
  id_selected <- rownames(marker_filtered)[nrow(marker_filtered)]
  markers.clusters.object$total_markers <- unique(data.correlation.all.CA.orig.celltype$gene[1:id_selected])
  
  
  #Parameters for the simulation
  single_cell_data_exp <- eset.sc.sparse_CA
  hierarchy = hierarchical_clustering_CA.orig.celltype$`1`
  bulk_data = data.frame(as.matrix(eset.sparse_pseudobulk_CA@assayData$exprs))
  sub_clusters_var <- 'orig.celltype'
  sample <- 'Disorder' #just one value: Control
  use_min_cor_strategy = TRUE
  delete_shared_level_markers = FALSE
  delete_shared_internal_markers = FALSE
  ordering_strategy = 'foldchange_pvalue'
  
  start_time <- Sys.time()
  
  #let's run de simulation
  deco.actual.data.CA.orig.celltype_value <- 
    run_deconvolution_tree_guided_recursive(bulk_data = bulk_data,
                                            single_cell_data_exp = single_cell_data_exp, 
                                            hierarchy = hierarchy,
                                            sub_clusters_var = sub_clusters_var,
                                            sample = sample, 
                                            use_min_cor_strategy = use_min_cor_strategy, 
                                            ordering_strategy = ordering_strategy,
                                            delete_shared_level_markers = delete_shared_level_markers, 
                                            delete_shared_internal_markers = delete_shared_internal_markers,
                                            param.logfc.threshold = 1.0,
                                            param.p_val_adj = 0.05,
                                            minimum_markers = 4,
                                            min_delta_cor_threshold = 0.0,
                                            delta_threshold = 1e-15,
                                            max_iterations = max_iterations,
                                            markers.clusters.parameter = markers.clusters.object,
                                            verbose = TRUE)
  
  
  end_time <- Sys.time()
  end_time - start_time
  
  return(deco.actual.data.CA.orig.celltype_value)
}


