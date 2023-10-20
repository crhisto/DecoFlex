
#parameteres:
root.dir <- "/mnt_volumen/GIT_REPOSITORIES/DecoFlex/DecoFlexR/vignettes/results/exp_0/"
iterations_number <- 15000
core_number <- 10

#importing functinos to parallelize
library(parallel)
library(doParallel)

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

#Including generics functions.
source("/mnt_volumen/GIT_REPOSITORIES/DecoFlex/DecoFlexR/vignettes/background_scripts/generic_scripts_experiments.R")

drop.exp0.complete <- drop.sub.exp0.complete <- list()

celltypes.exp0.merge_1_2 <- c("OPCs", "Endothelial", "Neurons", "Oligodendrocytes", "Astrocytes", "Microglia", 'Neurons_Exc', 'Neurons_Inh')

celltypes.exp0 <- c("OPCs", "Endothelial", "Neurons", "Oligodendrocytes", "Astrocytes", "Microglia")


print('Executing parallelized function...')

if(is.null(core_number)){
  # Calculate the number of cores
  no_cores <- detectCores() - 1
}else{
  no_cores <- core_number
}

print(paste0('Number of cores to use: ', no_cores))

registerDoParallel(no_cores)

drop.exp0.complete <- NULL

# For each cell-type dropped in the merge1 group, I will run the deconvolution
# Definition of the parallel function
drop.exp0.complete <- foreach(j = 1:length(celltypes.exp0.merge_1_2),
                   .combine = c)  %dopar%
  {
    
    celltype_name <- celltypes.exp0.merge_1_2[j]

    drop.exp0.complete.temp <- NULL
    
    if(celltype_name %in% celltypes.exp0){
      
      print(paste0('Running Decoflex deleting: ', celltype_name))
      
      hierarchy.celltype.deleted <- delete_celltype_hierarchy(hierarchy = hierarchical_clustering_VL.merge1$`1`,
                                                              celltype_name_delete = celltype_name)
      print(Sys.time())
      drop.exp0.complete.temp <- run.decoflex(single_cell_data_exp = eset.sc.sparse_VL,
                                              sub_clusters_var = 'merge1',
                                              hierarchy = hierarchy.celltype.deleted,
                                              bulk_data = data.frame(as.matrix(eset.sparse_pseudobulk_VL@assayData$exprs)),
                                              max_iterations = iterations_number)
      
    }else if(celltype_name == 'Neurons_Exc'){

      #Deleting the celltype in the hierarchy: Neurons_Exc
      hierarchy.celltype.deleted_Neurons_Exc  <- delete_celltype_hierarchy(hierarchy = hierarchical_clustering_VL.merge2$`1`,
                                                                           celltype_name_delete = "Neurons_Exc")
      drop.exp0.complete.temp <- run.decoflex(single_cell_data_exp=eset.sc.sparse_VL,
                                                         sub_clusters_var = 'merge2',
                                                         hierarchy = hierarchy.celltype.deleted_Neurons_Exc,
                                                         bulk_data = data.frame(as.matrix(eset.sparse_pseudobulk_VL@assayData$exprs)),
                                                         max_iterations = iterations_number)
      
    }else if(celltype_name == 'Neurons_Inh'){
      
      #Deleting the celltype in the hierarchy: Neurons_Inh
      hierarchy.celltype.deleted_Neurons_Inh  <- delete_celltype_hierarchy(hierarchy = hierarchical_clustering_VL.merge2$`1`,
                                                                           celltype_name_delete = "Neurons_Inh")
      drop.exp0.complete.temp<- run.decoflex(single_cell_data_exp=eset.sc.sparse_VL,
                                                         sub_clusters_var = 'merge2',
                                                         hierarchy = hierarchy.celltype.deleted_Neurons_Inh,
                                                         bulk_data = data.frame(as.matrix(eset.sparse_pseudobulk_VL@assayData$exprs)),
                                                         max_iterations = iterations_number)
    }
    
    #Saving the object with the list of results
    save(drop.exp0.complete.temp, file = paste0(root.dir, celltype_name, "_Dropped_sub.exp0_Estimates.rda"))
    
    assign(celltype_name, drop.exp0.complete.temp)
    
    #Returning the object
    get0(celltype_name, ifnotfound = paste0('Object not found: ', celltype_name)) 
  }
stopImplicitCluster()

#Saving everything
save(drop.exp0.complete, file = paste0(root.dir, "all_Dropped_sub.exp0_Estimates.rda"))

