# Deconvolution using different references to deconvolute: dm
#parameteres:
root.dir <- "/mnt_volumen/GIT_REPOSITORIES/DecoFlex/DecoFlexR/vignettes/results/exp_7/"
iterations_number <- 15000
core_number <- 6

#importing functions to parallelize
library(parallel)
library(doParallel)
library(DecoFlex)

#Loading the needed files: eset.sc.sparse_VL, deco.actual.data.VL.merge1, deco.actual.data.VL.merge2, eset.sparse_pseudobulk_VL ...
load(paste0(root.dir, "resources_exp7.rData"))


list_cell_types_dm <- c("Astrocytes", "Endothelial", "Microglia", "Oligodendrocytes", "Neurons")


# Initially for the top celltypes: 
hierarchical_clustering_dm <- list(
  '1'=list(tree_level=0, leaf_number=1, celltype_list = list_cell_types_dm, 
           next_level_clustering=NULL
  )
)


#2. Including generics functions.
source("/mnt_volumen/GIT_REPOSITORIES/DecoFlex/DecoFlexR/vignettes/background_scripts/generic_scripts_experiments.R")

drop.exp7.complete <- drop.sub.exp7.complete <- list()

#Order based on celltype size
celltypes.exp7.merge1 <- c("Astrocytes", "Endothelial", "Microglia", "Oligodendrocytes", "Neurons")

print('Executing parallelized function...')

if(is.null(core_number)){
  # dmlculate the number of cores
  no_cores <- detectCores() - 1
}else{
  no_cores <- core_number
}

print(paste0('Number of cores to use: ', no_cores))

registerDoParallel(no_cores)

drop.exp7.complete <- NULL

# For each cell-type dropped in the merge1 group, I will run the deconvolution
# Definition of the parallel function
drop.exp7.complete <- foreach(j = 1:length(sc.ref.object),
                   .combine = c)  %dopar%
  {
    
    dataset_name <- names(sc.ref.object)[j]
    dataset_value <- sc.ref.object[[dataset_name]]

    drop.exp7.complete.temp <- NULL
    
    print(paste0('Running Decoflex merge1 dataset: ', dataset_name))
    
    #1. Checking the celltypes that are missing in the dataset
    missing_celltypes <- setdiff(celltypes.exp7.merge1, unique(dataset_value$merge1))
    
    #2. Modifying the hierarchy based on the celltypes that are present in the dataset
    hierarchy.celltype.modified <- NULL
    if(length(missing_celltypes)>0){
      print(paste0('Missing celltypes: ', paste(shQuote(missing_celltypes), collapse=", ")))
      hierarchy.celltype.modified <- delete_celltype_hierarchy(hierarchy = hierarchical_clustering_dm$`1`,
                                                               celltype_name_delete = missing_celltypes)
    }else{
      print(paste0('None missing celltypes.'))
      hierarchy.celltype.modified <- hierarchical_clustering_dm$`1`
    }
    
    #3. Creation of single cell data object (reference)
    sparse_matrix_dataset <- as(as.matrix(dataset_value@assays$RNA@counts), "dgCMatrix")
    fdata_split_dataset_dataset <- rownames(dataset_value@assays$RNA@counts)
    pdata_split_dataset_dataset <- cbind(orig.celltype = as.character(dataset_value$orig.celltype),
                                    Disorder = as.character("Control"),
                                    merge1 = as.character(dataset_value$merge1))
    eset.sc.sparse_dataset <- getESET(sparse_matrix_dataset, fdata = fdata_split_dataset_dataset, pdata = pdata_split_dataset_dataset)  
      
    
    #4. Running the deconvolution: I use eset.sparse_pseudobulk_dm.symbol because the original one has IDs instead of symbols.
    drop.exp7.complete.temp <- run.decoflex(single_cell_data_exp = eset.sc.sparse_dataset,
                                            sub_clusters_var = 'merge1',
                                            hierarchy = hierarchy.celltype.modified,
                                            bulk_data = data.frame(as.matrix(eset.sparse_pseudobulk_dm.symbol@assayData$exprs)),
                                            param.logfc.threshold = 0.1,
                                            max_iterations = iterations_number)
    
    #Saving the object with the list of results
    save(drop.exp7.complete.temp, file = paste0(root.dir, dataset_name, "_reference_dm_input.exp7_Estimates.rda"))
    
    list_values_model <- list()
    list_values_model[[dataset_name]] <- drop.exp7.complete.temp
    
    
    assign(dataset_name, list_values_model)
    
    #Returning the object
    get0(dataset_name, ifnotfound = paste0('Object not found: ', dataset_name)) 
  }
stopImplicitCluster()

#Saving everything
save(drop.exp7.complete, file = paste0(root.dir, "all_reference_dm_input.exp7_Estimates.rda"))

