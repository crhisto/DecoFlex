
# Note: this execution using background jobs test the functionlity where
# DecoFlex caclulates in the to level of the hierarchy the marker genes for 
# all Celltypes using OMiC process and then in each level, the same set of 
# markers is used fintering the specific celtypes that the level has.

list_cell_types_VL.merge1 <- c("OPCs", "Endothelial", "Oligodendrocytes", "Astrocytes", "Microglia", "Neurons")

# Initially for the top celltypes: 
hierarchical_clustering_VL.merge1 <- list(
  '1'=list(tree_level=0, leaf_number=1, celltype_list = list_cell_types_VL.merge1, 
           next_level_clustering=NULL
  )
)


library(DecoFlex)

#Parameters for the simulation
single_cell_data_exp <- eset.sc.sparse_VL
hierarchy = hierarchical_clustering_VL.merge1$`1`
bulk_data = data.frame(as.matrix(eset.sparse_pseudobulk_VL@assayData$exprs))
sub_clusters_var <- 'merge1'
sample <- 'Disorder' #just one value: Control
use_min_cor_strategy = TRUE
delete_shared_level_markers = FALSE
delete_shared_internal_markers = FALSE
ordering_strategy = 'foldchange_pvalue'

start_time <- Sys.time()

#let's run de simulation
deco.actual.data.VL.merge1.global.markers <- 
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
                                          max_iterations = 1,
                                          use_global_makers = TRUE,
                                          verbose = TRUE)
end_time <- Sys.time()
end_time - start_time