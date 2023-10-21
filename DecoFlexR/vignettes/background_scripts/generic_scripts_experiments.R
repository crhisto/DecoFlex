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


#Calculate final correlation
calculate_correlation_CA.orig.celltype <- function(list_cell_types_CA.orig.celltype, deco.actual.data.CA.orig.celltype, true.CA){
  
  proportions.temp <- t(deco.actual.data.CA.orig.celltype$back_propagation_proportions_top_detailed)
  colnames(proportions.temp) <- list_cell_types_CA.orig.celltype
  cor.calc.vs.true.CA.orig.celltype <- cor(x=proportions.temp[, list_cell_types_CA.orig.celltype], 
                                           y=true.CA$orig.celltype[,list_cell_types_CA.orig.celltype], 
                                           method='pearson')
  
  cor.calc.vs.true.CA.orig.celltype
  
  return(mean(diag(cor.calc.vs.true.CA.orig.celltype)))
}

list_cell_types_CA.orig.celltype <- c('Astro L1-6 FGFR3 SLC14A1', 'Exc L2 LAMP5 LTK', 'Exc L2-3 LINC00507 FREM3', 'Exc L3-4 RORB CARM1P1', 'Exc L3-5 RORB ESR1', 'Exc L4-5 RORB FOLH1B', 'Exc L4-6 FEZF2 IL26', 'Exc L4-6 RORB SEMA3E', 'Exc L5-6 FEZF2 ABO', 'Exc L5-6 FEZF2 EFTUD1P1', 'Exc L5-6 THEMIS C1QL3', 'Inh L1 SST NMBR', 'Inh L1-3 SST CALB1', 'Inh L1-4 LAMP5 LCP2', 'Inh L2-4 PVALB WFDC2', 'Inh L2-6 LAMP5 CA1', 'Oligo L1-6 OPALIN', 'OPC L1-6 PDGFRA')

metrics.CA.orig.celltype <- NULL

for(model_counter in 1:length(list_deco_CA_orig.celltype)){
  model_name <- names(list_deco_CA_orig.celltype)[model_counter]
  current_model <- list_deco_CA_orig.celltype[[model_name]]
  print(model_name)
  correlation_model <- calculate_correlation_CA.orig.celltype(list_cell_types_CA.orig.celltype, current_model, true.CA)
  print(paste0('', correlation_model))
  
  last_iteration <- current_model$result_deco_top_cluster$running_info[nrow(current_model$result_deco_top_cluster$running_info),]
  info_model <- data.frame("marker_count" = length(current_model$markers.top.clusters.object$total_markers),
                           "r_correlation" = correlation_model, 
                           "iteration" = last_iteration$iteration, 
                           "delta_divergence_value" = last_iteration$delta_divergence_value)
  
  if(is.null(metrics.CA.orig.celltype)){
    metrics.CA.orig.celltype <- info_model
  }else{
    metrics.CA.orig.celltype <- rbind(metrics.CA.orig.celltype, info_model)
  }
}


create_semi_reference_objects <- function(extra_unknown_celltypes = 1, cell_type_names, bulk.data_mixtures.brain, 
                                          w_fixed.value, marker_genes, extra_marker_genes_semireference.value = c(), version = 'three'){
  
  
  #1. First, I will delete the marker genes from the additional ones to avoid problems in the x, w construction
  extra_marker_genes_semireference <- extra_marker_genes_semireference.value[!(extra_marker_genes_semireference.value %in% marker_genes)]
  
  
  #2. I will get the final w signature matrix, with the fixed markers such: markers + extra markers
  # The problem is that the reference could have all markers and all the extra markers...but maybe not in the extra marker case.
  shared_markers_vs_w_reference <- intersect(rownames(w_fixed.value), extra_marker_genes_semireference)
  extra_marker_genes_semireference.known <- shared_markers_vs_w_reference
  #In case that all extra markers are in the w, then the unknown will be zero
  extra_marker_genes_semireference.unknown <- extra_marker_genes_semireference[!(extra_marker_genes_semireference %in% extra_marker_genes_semireference.known)]
  
  w_fixed <- w_fixed.value[c(marker_genes, 
                             extra_marker_genes_semireference.known,
                             extra_marker_genes_semireference.unknown
  ), ]
  
  
  #3. Now I build the bulk data with the marker genes: first the normal ones and then the additionals.
  bulk.data_mixtures.brain.filtered <- bulk.data_mixtures.brain[c(marker_genes, 
                                                                  extra_marker_genes_semireference.known,
                                                                  extra_marker_genes_semireference.unknown), ]
  
  #4. Calculation of some numbers important for the construction.
  number_cell_types <- length(cell_type_names)
  samples_bulk_data <- ncol(bulk.data_mixtures.brain) #100
  bulk.column_names <- colnames(bulk.data_mixtures.brain)
  
  # celltypes=5 x 100 columns (samples)
  # I have n genes as rows and m+1 columns, being the last one the unknown one. 
  mask_h <- matrix(FALSE, nrow = (number_cell_types + extra_unknown_celltypes), ncol = samples_bulk_data)
  colnames(mask_h) <- bulk.column_names
  rownames(mask_h) <- c(cell_type_names, paste0('unknown_', 1:extra_unknown_celltypes))
  
  
  # In this case, I can do two things: 1. Fix just the celltypes that I have as reference and the unknown leave it as black. 2. Assign the unknown cell-type as zero, since it suppose to be a WT, and then in the other conditions it will appear a proportion. Also the Tumor DNMT1 KO, will be similar to WT.
  if(version=='one'){
    #Version 1.
    #mask_h[1:3,1:3] <- TRUE
  }else if(version=='two'){
    #Version 1.
    #mask_h[1:(3+extra_unknown_celltypes),1:3] <- TRUE
    #Fixing all the proportions
    #mask_h[1:3,4:12] <- TRUE
  }else if(version=='three'){
    #nothing is fixed in h, therefore all values will be optimized
  }
  
  
  #Now let's put the proportions that I needed, filling with zeros the rest
  partial_h_fixed <- matrix(0, nrow = (number_cell_types + extra_unknown_celltypes), ncol = samples_bulk_data)
  colnames(partial_h_fixed) <- bulk.column_names
  rownames(partial_h_fixed) <- c(cell_type_names, paste0('unknown_', 1:extra_unknown_celltypes))
  
  #In case we add known proportions.
  if(version!='three'){
    #let's take the proportion
    # proportions_wt_single_cell <- as.data.frame(deco.simulation.thymus.tree_guided.recursive.reduced.results.top$pseudo.bulk.data$truep) 
    # colnames(proportions_wt_single_cell) <- c('sample_6_8W')
    # proportions_wt_single_cell <- t(proportions_wt_single_cell)
    # proportions_wt_single_cell['sample_6_8W', paste0('cluster_', c(1:3)) ]
    # 
    # #Since I have three replicates I have to fix the proportions three times
    # partial_h_fixed[1:3,1] <- as.matrix(proportions_wt_single_cell['sample_6_8W', paste0('cluster_', c(1:3)) ])
    # partial_h_fixed[1:3,2] <- as.matrix(proportions_wt_single_cell['sample_6_8W', paste0('cluster_', c(1:3)) ])
    # partial_h_fixed[1:3,3] <- as.matrix(proportions_wt_single_cell['sample_6_8W', paste0('cluster_', c(1:3)) ])
  }
  
  
  
  
  #Now I will create the W
  mask_w <- matrix(TRUE, 
                   nrow = (length(marker_genes) + length(extra_marker_genes_semireference)), 
                   ncol = (number_cell_types + extra_unknown_celltypes))
  colnames(mask_w) <- c(cell_type_names, paste0('unknown_', 1:extra_unknown_celltypes))
  rownames(mask_w) <- c(marker_genes, extra_marker_genes_semireference.known, extra_marker_genes_semireference.unknown)
  
  #Now I will mark as FALSE (CALCULATE) the marker genes corresponding to the unknown celltypes values. I took 
  #w_fixed since it has all the markers.
  mask_w[1:nrow(w_fixed), (number_cell_types + 1):(number_cell_types + extra_unknown_celltypes)] <- FALSE
  
  #If I have extra genes "unknown", I have to add false for the unknown celltytes
  if(length(extra_marker_genes_semireference.unknown) > 0){
    mask_w[(length(marker_genes)+length(extra_marker_genes_semireference.known)+ 1):nrow(mask_w), 1:(number_cell_types + extra_unknown_celltypes)] <- FALSE
  }
  
  # I have k cell-types rows by m samples.
  partial_w_fixed <- matrix(0, 
                            nrow = (length(marker_genes) + length(extra_marker_genes_semireference)), 
                            ncol = (number_cell_types + extra_unknown_celltypes))
  colnames(partial_w_fixed) <- c(cell_type_names, paste0('unknown_', 1:extra_unknown_celltypes))
  rownames(partial_w_fixed) <- c(marker_genes, extra_marker_genes_semireference.known, extra_marker_genes_semireference.unknown)
  
  partial_w_fixed[1:nrow(w_fixed), 1:number_cell_types] <- as.matrix(w_fixed)  
  
  return(list(mask_h = mask_h, partial_h_fixed = partial_h_fixed,
              mask_w = mask_w, partial_w_fixed = partial_w_fixed, 
              bulk_data = bulk.data_mixtures.brain.filtered))
}


run_deconvolution_decoflex_semi_reference <- function(bulk.data_mixtures.brain, 
                                                      partial_w_fixed, partial_h_fixed, 
                                                      mask_w, mask_h,
                                                      number_cell_types, extra_unknown_celltypes,
                                                      proportion_constraint_h = TRUE, max_iterations = 1000, delta_threshold = 1e-15){
  
  total_celltypes <- (number_cell_types + extra_unknown_celltypes)
  print(paste0('k=', total_celltypes))
  
  start_time <- Sys.time()
  
  # I run the deconvolution without any other parameter
  deco.results <- run_complete_deconvolution(
    x_matrix = data.frame(bulk.data_mixtures.brain),
    y_matrix = NULL,
    z_matrix = NULL,
    k = total_celltypes,
    alpha = 0, 
    beta = 0,
    max_iterations = max_iterations, 
    delta_threshold = delta_threshold, 
    proportion_constraint_h = proportion_constraint_h,
    partial_w_fixed=data.frame(partial_w_fixed),
    partial_h_fixed=data.frame(partial_h_fixed),
    w_mask_fixed=data.frame(mask_w),
    h_mask_fixed=data.frame(mask_h),
    batches_partial_fixed=1)
  
  end_time <- Sys.time()
  end_time - start_time
  
  return(deco.results)
}

plot_correlation_iterations <- function(data.correlation.all, 
                                        zoom.xlim = c(0, 55),
                                        red_line = 28,
                                        black_line = 155, 
                                        color_lines = c('red', 'black')){
  correlation.vs.markers.selected.plot <- ggplot(data.correlation.all, aes(x = marker_count, y = joint_cor, group = celltype, color = celltype, fill = celltype)) +
    #geom_line() +
    geom_area()+
    theme_bw() +
    theme(panel.border = invis, axis.line = element_line(), axis.title = element_text(size = 8), 
          axis.text = element_text(size = 8), legend.position = "right") + 
    ylim(-1,1) + 
    labs(x = 'Batch number of marker genes selected incrementally', y = 'Join correlation between celltypes')
  
  if(!is.null(red_line)){
    correlation.vs.markers.selected.plot <- correlation.vs.markers.selected.plot + geom_vline(xintercept=red_line,linetype = "dashed", colour = color_lines[1])
  }
  
  if(!is.null(black_line)){
    correlation.vs.markers.selected.plot <- correlation.vs.markers.selected.plot + geom_vline(xintercept=black_line,linetype = "dashed", colour = color_lines[2])
  }
  
  if(!is.null(zoom.xlim)){
    correlation.vs.markers.selected.plot <- correlation.vs.markers.selected.plot + ggforce::facet_zoom(xlim = zoom.xlim)
  }
  
  return(correlation.vs.markers.selected.plot)
}


create_dataset_markers_multiple_levels <- function(result_level_object, level_list){
  
  #all the genes por level
  top.level.marker.genes <- result_level_object$markers.top.clusters.object$total_markers
  
  top.markers.detailed_all = NULL
  top.markers.detailed_with_repetition = NULL
  for (counter_level in level_list) {
    #Getting the indivigual groups for level
    top.level.marker.genes.subcluster <- rownames(result_level_object$markers.top.clusters.object$list_markers[[counter_level]])
    #Creation of one dataframe with all genes repited
    top.markers.detailed <- data.frame(cbind(gene_name = intersect(top.level.marker.genes.subcluster, top.level.marker.genes), group = paste0('', counter_level)))
    
    
    #Version without repetition and deleting that from subcluster 2 to 1.
    top.markers.detailed.filtered <- top.markers.detailed[!(top.markers.detailed$gene_name %in% top.markers.detailed_all$gene_name), ]
    
    #without_repetition
    if(is.null(top.markers.detailed_all)){
      top.markers.detailed_all <-data.frame(top.markers.detailed.filtered)
      
    }else{
      top.markers.detailed_all <- rbind(top.markers.detailed_all, top.markers.detailed.filtered)
    }
    
    #Version with repetition to check the sharing genes between celltypes
    if(is.null(top.markers.detailed_with_repetition)){
      top.markers.detailed_with_repetition <-data.frame(top.markers.detailed)
      
    }else{
      top.markers.detailed_with_repetition <- rbind(top.markers.detailed_with_repetition, top.markers.detailed)
    }
  }
  
  
  #Checking repetition
  top.markers.detailed_all$repetition <- !duplicated(top.markers.detailed_all$gene_name)
  
  #Those genes doesn't exist in bulk data.
  delete.genes <- c()
  top.markers.detailed_all <- top.markers.detailed_all[!(top.markers.detailed_all$gene_name %in% delete.genes),]
  
  #Order of groups
  top.markers.detailed_all$group <- factor(x = top.markers.detailed_all$group, levels = unique(top.markers.detailed_all$group))
  rownames(top.markers.detailed_all) <- top.markers.detailed_all$gene_name
  
  return(list(top.markers.detailed_all = top.markers.detailed_all, 
              top.markers.detailed_with_repetition = top.markers.detailed_with_repetition))
}



run_deco_DM_merge1 <- function(data.correlation.all.DM.merge1, marker_value, max_iterations=25000){
  
  markers.clusters.object <- list()
  markers.clusters.object$accumulative_corr_markers <- data.correlation.all.DM.merge1
  markers.clusters.object$list_markers <- NULL
  markers.clusters.object$list_markers.percentile <- NULL
  
  #Adding just the markers that I'm interesting on.
  marker_filtered <- data.correlation.all.DM.merge1[data.correlation.all.DM.merge1$current_marker == marker_value,]
  id_selected <- rownames(marker_filtered)[nrow(marker_filtered)]
  markers.clusters.object$total_markers <- unique(data.correlation.all.DM.merge1$gene[1:id_selected])
  
  
  #Parameters for the simulation
  single_cell_data_exp <- eset.sc.sparse_dm
  hierarchy = hierarchical_clustering_dm$`1`
  bulk_data = data.frame(as.matrix(eset.sparse_pseudobulk_dm@assayData$exprs))
  sub_clusters_var <- 'orig.celltype'
  sample <- 'Disorder' #just one value: Control
  use_min_cor_strategy = TRUE
  delete_shared_level_markers = FALSE
  delete_shared_internal_markers = FALSE
  ordering_strategy = 'foldchange_pvalue'
  
  start_time <- Sys.time()
  
  #let's run de simulation
  deco.actual.data.DM.merge1_value <- 
    run_deconvolution_tree_guided_recursive(bulk_data = bulk_data,
                                            single_cell_data_exp = single_cell_data_exp, 
                                            hierarchy = hierarchy,
                                            sub_clusters_var = sub_clusters_var,
                                            sample = sample, 
                                            use_min_cor_strategy = use_min_cor_strategy, 
                                            ordering_strategy = ordering_strategy,
                                            delete_shared_level_markers = delete_shared_level_markers, 
                                            delete_shared_internal_markers = delete_shared_internal_markers,
                                            param.logfc.threshold = 0.1,
                                            param.p_val_adj = 0.05,
                                            minimum_markers = 4,
                                            min_delta_cor_threshold = 0.0,
                                            delta_threshold = 1e-15,
                                            max_iterations = max_iterations,
                                            markers.clusters.parameter = markers.clusters.object,
                                            verbose = TRUE)
  
  
  end_time <- Sys.time()
  end_time - start_time
  
  return(deco.actual.data.DM.merge1_value)
}


#Calculate final correlation
calculate_correlation_DM.merge1 <- function(list_cell_types_DM.merge1, deco.actual.data.dm, true.dm){
  cor.calc.vs.true.DM.merge1 <- cor(x=t(deco.actual.data.dm$back_propagation_proportions_top_detailed)[, list_cell_types_DM.merge1], 
                                    y=true.dm$merge2[,list_cell_types_DM.merge1], 
                                    method='pearson')
  
  cor.calc.vs.true.DM.merge1
  
  return(mean(diag(cor.calc.vs.true.DM.merge1)))
}

list_cell_types_DM.merge1 <- c("Astrocytes", "Endothelia", "Microglia", "Oligodendrocytes", "Neurons")

metrics.DM.merge1 <- NULL

for(model_counter in 1:length(list_deco_DM_merge1)){
  model_name <- names(list_deco_DM_merge1)[model_counter]
  current_model <- list_deco_DM_merge1[[model_name]]
  print(model_name)
  correlation_model <- calculate_correlation_DM.merge1(list_cell_types_DM.merge1, current_model, true.dm)
  print(paste0('', correlation_model))
  
  last_iteration <- current_model$result_deco_top_cluster$running_info[nrow(current_model$result_deco_top_cluster$running_info),]
  info_model <- data.frame("marker_count" = length(current_model$markers.top.clusters.object$total_markers),
                           "r_correlation" = correlation_model, 
                           "iteration" = last_iteration$iteration, 
                           "delta_divergence_value" = last_iteration$delta_divergence_value)
  
  if(is.null(metrics.DM.merge1)){
    metrics.DM.merge1 <- info_model
  }else{
    metrics.DM.merge1 <- rbind(metrics.DM.merge1, info_model)
  }
}


run_deco_VL_merge2 <- function(data.correlation.all.VL.merge2, marker_value, max_iterations=25000){
  
  markers.clusters.object <- list()
  markers.clusters.object$accumulative_corr_markers <- data.correlation.all.VL.merge2
  markers.clusters.object$list_markers <- NULL
  markers.clusters.object$list_markers.percentile <- NULL
  
  #Adding just the markers that I'm interesting on.
  marker_filtered <- data.correlation.all.VL.merge2[data.correlation.all.VL.merge2$current_marker == marker_value,]
  id_selected <- rownames(marker_filtered)[nrow(marker_filtered)]
  markers.clusters.object$total_markers <- unique(data.correlation.all.VL.merge2$gene[1:id_selected])
  
  
  #Parameters for the simulation
  single_cell_data_exp <- eset.sc.sparse_VL
  hierarchy = hierarchical_clustering_VL.merge2$`1`
  bulk_data = data.frame(as.matrix(eset.sparse_pseudobulk_VL@assayData$exprs))
  sub_clusters_var <- 'merge2'
  sample <- 'Disorder' #just one value: Control
  use_min_cor_strategy = TRUE
  delete_shared_level_markers = FALSE
  delete_shared_internal_markers = FALSE
  ordering_strategy = 'foldchange_pvalue'
  
  start_time <- Sys.time()
  
  #let's run de simulation
  deco.actual.data.VL.merge2_value <- 
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
  
  return(deco.actual.data.VL.merge2_value)
}


#Calculate final correlation
calculate_correlation_VL <- function(list_cell_types_VL.merge2, deco.actual.data.VL.merge2, true.VL){
  #fixing output
  rownames(deco.actual.data.VL.merge2$back_propagation_proportions_top_detailed)[which(rownames(deco.actual.data.VL.merge2$back_propagation_proportions_top_detailed) == "Neurons_Exc")] <-"Excitatory"
  rownames(deco.actual.data.VL.merge2$back_propagation_proportions_top_detailed)[which(rownames(deco.actual.data.VL.merge2$back_propagation_proportions_top_detailed) == "Neurons_Inh")] <-"Inhibitory"
  
  
  cor.calc.vs.true.VL.merge2 <- cor(x=t(deco.actual.data.VL.merge2$back_propagation_proportions_top_detailed)[, list_cell_types_VL.merge2], 
                                    y=true.VL$merge2[,list_cell_types_VL.merge2], 
                                    method='pearson')
  
  cor.calc.vs.true.VL.merge2
  
  return(mean(diag(cor.calc.vs.true.VL.merge2)))
}

list_cell_types_VL.merge2 <- c("Astrocytes", "Excitatory", "Oligodendrocytes", "OPCs", "Microglia", "Inhibitory", "Endothelial")

metrics.VL.merge2 <- NULL

for(model_counter in 1:length(list_deco_VL_merge2)){
  model_name <- names(list_deco_VL_merge2)[model_counter]
  current_model <- list_deco_VL_merge2[[model_name]]
  print(model_name)
  correlation_model <- calculate_correlation_VL(list_cell_types_VL.merge2, current_model, true.VL)
  print(paste0('', correlation_model))
  
  last_iteration <- current_model$result_deco_top_cluster$running_info[nrow(current_model$result_deco_top_cluster$running_info),]
  info_model <- data.frame("marker_count" = length(current_model$markers.top.clusters.object$total_markers),
                           "r_correlation" = correlation_model, 
                           "iteration" = last_iteration$iteration, 
                           "delta_divergence_value" = last_iteration$delta_divergence_value)
  
  if(is.null(metrics.VL.merge2)){
    metrics.VL.merge2 <- info_model
  }else{
    metrics.VL.merge2 <- rbind(metrics.VL.merge2, info_model)
  }
}
