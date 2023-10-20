list_cell_types_VL.merge1 <- c("OPCs", "Endothelial", "Oligodendrocytes", "Astrocytes", "Microglia", "Neurons")

# Initially for the top celltypes: 
hierarchical_clustering_VL.merge1 <- list(
  '1'=list(tree_level=0, leaf_number=1, celltype_list = list_cell_types_VL.merge1, 
           next_level_clustering=NULL
  )
)

list_cell_types_VL.merge2 <- c("OPCs", "Endothelial", "Neurons_Inh", "Neurons_Exc", "Oligodendrocytes", "Astrocytes", "Microglia")

# Initially for the top celltypes: 
hierarchical_clustering_VL.merge2 <- list(
  '1'=list(tree_level=0, leaf_number=1, celltype_list = list_cell_types_VL.merge2, 
           next_level_clustering=NULL
  )
)



source("/mnt_volumen/GIT_REPOSITORIES/DecoFlex/DecoFlexR/vignettes/background_scripts/generic_scripts_experiments.R")
root.dir <- "/mnt_volumen/GIT_REPOSITORIES/DecoFlex/DecoFlexR/vignettes/results/"
iterations_number <- 2

drop.exp0.complete <- drop.sub.exp0.complete <- list()

celltypes.exp0 <- c("OPCs", "Endothelial", "Neurons", "Oligodendrocytes", "Astrocytes", "Microglia")

#For each cell-type dropped in the merge1 group, I will run the deconvolution
for (j in celltypes.exp0) {
  print(paste0('Running Decoflex deleting: ', j))
  hierarchy.celltype.deleted <- delete_celltype_hierarchy(hierarchy = hierarchical_clustering_VL.merge1$`1`,
                                                          celltype_name_delete = j)
  print(Sys.time())
  drop.exp0.complete[[j]] <- run.decoflex(single_cell_data_exp = eset.sc.sparse_VL,
                                     sub_clusters_var = 'merge1',
                                     hierarchy = hierarchy.celltype.deleted,
                                     bulk_data = data.frame(as.matrix(eset.sparse_pseudobulk_VL@assayData$exprs)),
                                     max_iterations = iterations_number)
  
  #Saving the object with the list of results
  save(drop.exp0.complete, file = paste0(root.dir, "Dropped_sub.exp0_Estimates.rda"))
}


#Deleting the celltype in the hierarchy: Neurons_Exc
hierarchy.celltype.deleted_Neurons_Exc  <- delete_celltype_hierarchy(hierarchy = hierarchical_clustering_VL.merge2$`1`,
                                                                     celltype_name_delete = "Neurons_Exc")
drop.sub.exp0.complete$Neurons_Exc <- run.decoflex(single_cell_data_exp=eset.sc.sparse_VL,
                                              sub_clusters_var = 'merge2',
                                              hierarchy = hierarchy.celltype.deleted_Neurons_Exc,
                                              bulk_data = data.frame(as.matrix(eset.sparse_pseudobulk_VL@assayData$exprs)),
                                              max_iterations = iterations_number)

#Saving the object with the list of results
save(drop.exp0.complete, file = paste0(root.dir, "Dropped_sub.exp0_Estimates.rda"))

#Deleting the celltype in the hierarchy: Neurons_Inh
hierarchy.celltype.deleted_Neurons_Inh  <- delete_celltype_hierarchy(hierarchy = hierarchical_clustering_VL.merge2$`1`,
                                                                     celltype_name_delete = "Neurons_Inh")
drop.sub.exp0.complete$Neurons_Inh <- run.decoflex(single_cell_data_exp=eset.sc.sparse_VL,
                                              sub_clusters_var = 'merge2',
                                              hierarchy = hierarchy.celltype.deleted_Neurons_Inh,
                                              bulk_data = data.frame(as.matrix(eset.sparse_pseudobulk_VL@assayData$exprs)),
                                              max_iterations = iterations_number)

#Saving the object with the list of results
save(drop.exp0.complete, file = paste0(root.dir, "Dropped_sub.exp0_Estimates.rda"))
