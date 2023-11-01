library(DecoFlex)

#Parameters for the simulation
single_cell_data_exp <- eset.sc.sparse_BR.normal
hierarchy = hierarchical_clustering_BR.merge1$`1`
bulk_data = data.frame(as.matrix(eset.sparse_pseudobulk_BR@assayData$exprs))
sub_clusters_var <- 'merge1'
sample <- 'Disorder' #just one value: Control
use_min_cor_strategy = TRUE
delete_shared_level_markers = FALSE
delete_shared_internal_markers = FALSE
ordering_strategy = 'foldchange_pvalue'
markers.clusters.parameter = markers.clusters.object

start_time <- Sys.time()

#let's run de simulation
deco.actual.data.BR.merge1.opt <- 
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
                                          max_iterations = 15000,
                                          delta_threshold = 1e-15,
                                          markers.clusters.parameter = markers.clusters.parameter,
                                          verbose = TRUE)
#  mins
end_time <- Sys.time()
end_time - start_time

root.dir <- '/mnt_volumen/GIT_REPOSITORIES/DecoFlex/DecoFlexR/vignettes/figures_manuscript/'
save(deco.actual.data.BR.merge1.opt, file = paste0(root.dir, "deco.actual.data.BR.merge1.opt.rDAta"))