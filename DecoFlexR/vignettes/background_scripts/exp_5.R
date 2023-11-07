# Deconvolution using different references to deconvolute: VL

#parameteres:
root.dir <- "/mnt_volumen/GIT_REPOSITORIES/DecoFlex/DecoFlexR/vignettes/results/exp_5/"
iterations_number <- 15000
core_number <- 6

#importing functions to parallelize
library(parallel)
library(doParallel)
library(DecoFlex)

#Loading the needed files: eset.sc.sparse_VL, deco.actual.data.VL.merge1, deco.actual.data.VL.merge2, eset.sparse_pseudobulk_VL
load(paste0(root.dir, "resources_exp5.rData"))


list_cell_types_VL.merge1 <- c("OPCs", "Endothelial", "Oligodendrocytes", "Astrocytes", "Microglia", "Neurons")

#1. Adding initially for the top celltypes: 
hierarchical_clustering_VL.merge1 <- list(
  '1'=list(tree_level=0, leaf_number=1, celltype_list = list_cell_types_VL.merge1, 
           next_level_clustering=NULL
  )
)


#2. Including generics functions.
source("/mnt_volumen/GIT_REPOSITORIES/DecoFlex/DecoFlexR/vignettes/background_scripts/generic_scripts_experiments.R")

drop.exp5.complete <- drop.sub.exp5.complete <- list()

#Order based on celltype size
celltypes.exp5.merge1 <- c("OPCs", "Microglia", "Endothelial", "Astrocytes", "Oligodendrocytes", "Neurons")

print('Executing parallelized function...')

if(is.null(core_number)){
  # Calculate the number of cores
  no_cores <- detectCores() - 1
}else{
  no_cores <- core_number
}

print(paste0('Number of cores to use: ', no_cores))

registerDoParallel(no_cores)

drop.exp5.complete <- NULL

# For each cell-type dropped in the merge1 group, I will run the deconvolution
# Definition of the parallel function
drop.exp5.complete <- foreach(j = 1:length(sc.ref.object),
                   .combine = c)  %dopar%
  {
    
    dataset_name <- names(sc.ref.object)[j]
    dataset_value <- sc.ref.object[[dataset_name]]

    drop.exp5.complete.temp <- NULL
    
    print(paste0('Running Decoflex merge1 dataset: ', dataset_name))
    
    #1. Checking the celltypes that are missing in the dataset
    missing_celltypes <- setdiff(celltypes.exp5.merge1, unique(dataset_value$merge1))
    
    #2. Modifying the hierarchy based on the celltypes that are present in the dataset
    hierarchy.celltype.modified <- NULL
    if(length(missing_celltypes)>0){
      print(paste0('Missing celltypes: ', paste(shQuote(missing_celltypes), collapse=", ")))
      hierarchy.celltype.modified <- delete_celltype_hierarchy(hierarchy = hierarchical_clustering_VL.merge1$`1`,
                                                               celltype_name_delete = missing_celltypes)
    }else{
      print(paste0('None missing celltypes.'))
      hierarchy.celltype.modified <- hierarchical_clustering_VL.merge1$`1`
    }
    
    #3. Creation of single cell data object (reference)
    sparse_matrix_dataset <- as(as.matrix(dataset_value@assays$RNA@counts), "dgCMatrix")
    fdata_split_dataset_dataset <- rownames(dataset_value@assays$RNA@counts)
    pdata_split_dataset_dataset <- cbind(orig.celltype = as.character(dataset_value$orig.celltype),
                                    Disorder = as.character("Control"),
                                    merge1 = as.character(dataset_value$merge1))
    eset.sc.sparse_dataset <- getESET(sparse_matrix_dataset, fdata = fdata_split_dataset_dataset, pdata = pdata_split_dataset_dataset)  
      
    
    #4. Running the deconvolution
    drop.exp5.complete.temp <- run.decoflex(single_cell_data_exp = eset.sc.sparse_dataset,
                                            sub_clusters_var = 'merge1',
                                            hierarchy = hierarchy.celltype.modified,
                                            bulk_data = data.frame(as.matrix(eset.sparse_pseudobulk_VL@assayData$exprs)),
                                            param.logfc.threshold = 0.1,
                                            max_iterations = iterations_number)
    
    #Saving the object with the list of results
    save(drop.exp5.complete.temp, file = paste0(root.dir, dataset_name, "_reference_VL_input.exp5_Estimates.rda"))
    
    list_values_model <- list()
    list_values_model[[dataset_name]] <- drop.exp5.complete.temp
    
    
    assign(dataset_name, list_values_model)
    
    #Returning the object
    get0(dataset_name, ifnotfound = paste0('Object not found: ', dataset_name)) 
  }
stopImplicitCluster()

#Saving everything
save(drop.exp5.complete, file = paste0(root.dir, "all_reference_VL_input.exp5_Estimates.rda"))

