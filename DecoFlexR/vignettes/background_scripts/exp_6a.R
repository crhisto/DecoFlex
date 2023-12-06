# Deconvolution using different bulk references (F5, MM, IP) to deconvolute: CA

#parameteres:
root.dir <- "/mnt_volumen/GIT_REPOSITORIES/DecoFlex/DecoFlexR/vignettes/results/exp_6a/"
iterations_number <- 200
core_number <- 3

#importing functions to parallelize
library(parallel)
library(doParallel)
library(DecoFlex)

#Loading the needed files: eset.sc.sparse_CA, deco.actual.data.CA.merge1, deco.actual.data.CA.merge2, eset.sparse_pseudobulk_CA
load(paste0(root.dir, "resources_exp6.rData"))
load(paste0(root.dir, "resources_exp6_extra.rData"))
load(paste0(root.dir, "resources_exp0.rData"))


list_cell_types_CA.merge1 <- c("Neurons", "Oligodendrocytes", "OPCs", "Astrocytes")

#1. Adding initially for the top celltypes: 
hierarchical_clustering_CA.merge1 <- list(
  '1'=list(tree_level=0, leaf_number=1, celltype_list = list_cell_types_CA.merge1, 
           next_level_clustering=NULL
  )
)


#2. Including generics functions.
source("/mnt_volumen/GIT_REPOSITORIES/DecoFlex/DecoFlexR/vignettes/background_scripts/generic_scripts_experiments.R")

drop.exp6a.complete <- drop.sub.exp6a.complete <- list()

#Order based on celltype size
celltypes.exp6a.merge1 <- c("Neurons", "Oligodendrocytes", "OPCs", "Astrocytes")

print('Executing parallelized function...')

if(is.null(core_number)){
  # Calculate the number of cores
  no_cores <- detectCores() - 1
}else{
  no_cores <- core_number
}

print(paste0('Number of cores to use: ', no_cores))

registerDoParallel(no_cores)

drop.exp6a.complete <- NULL

list_bulk_signatures <- c('F5', 'MM', 'IP')

# For each cell-type dropped in the merge1 group, I will run the deconvolution
# Definition of the parallel function
drop.exp6a.complete <- foreach(j = 1:length(list_bulk_signatures),
                   .combine = c)  %dopar%
  {
    
    dataset_name <- list_bulk_signatures[j]
    dataset_value <- sigsSNME[[dataset_name]]
    
    #Fixing names dataset
    
    #replacing endothelia for endothelial
    colnames(dataset_value) <- replace(colnames(dataset_value), colnames(dataset_value)=='Endothelia', 'Endothelial')

    drop.exp6a.complete.temp <- NULL
    
    print(paste0('Running Decoflex merge1 dataset: ', dataset_name))
    
    #1. Checking the celltypes that are missing in the dataset
    missing_celltypes <- setdiff(celltypes.exp6a.merge1, colnames(dataset_value))

    
    #2. Modifying the hierarchy based on the celltypes that are present in the dataset
    hierarchy.celltype.modified <- NULL
    if(length(missing_celltypes)>0){
      print(paste0('Missing celltypes: ', paste(shQuote(missing_celltypes), collapse=", ")))
      hierarchy.celltype.modified <- delete_celltype_hierarchy(hierarchy = hierarchical_clustering_CA.merge1$`1`,
                                                               celltype_name_delete = missing_celltypes)
    }else{
      print(paste0('None missing celltypes.'))
      hierarchy.celltype.modified <- hierarchical_clustering_CA.merge1$`1`
    }
    
    #Celltypes that are going to be deconvoluted
    celltype_present <- hierarchy.celltype.modified$celltype_list
    print(paste0('Celltypes present: ', paste(shQuote(celltype_present), collapse=", ")))
    
    #Standart reference
    rownames(dataset_value) <- toupper(rownames(dataset_value))
    
    # 3. Gene Intersection between bulk and reference
    gene_intersection <- intersect(rownames(eset.sparse_pseudobulk_CA@assayData$exprs), rownames(dataset_value))

    # 4. Deconvolution with the partially fixed h and w
    drop.exp6a.complete.temp <- run_deconvolution_decoflex_semi_reference(
      bulk.data_mixtures.brain = data.frame(as.matrix(eset.sparse_pseudobulk_CA@assayData$exprs[gene_intersection,])),
      fixed_w = data.frame(as.matrix(dataset_value[gene_intersection, celltype_present])),
      number_cell_types = length(celltype_present), 
      proportion_constraint_h = TRUE,
      max_iterations = iterations_number,
      delta_threshold = 1e-15)
    

    #Saving the object with the list of results
    save(drop.exp6a.complete.temp, file = paste0(root.dir, dataset_name, "_reference_CA_input.exp6a_Estimates.rda"))
    
    list_values_model <- list()
    list_values_model[[dataset_name]] <- drop.exp6a.complete.temp
    
    
    assign(dataset_name, list_values_model)
    
    #Returning the object
    get0(dataset_name, ifnotfound = paste0('Object not found: ', dataset_name)) 
  }
stopImplicitCluster()

#Saving everything
save(drop.exp6a.complete, file = paste0(root.dir, "all_reference_CA_input.exp6a_Estimates.rda"))

