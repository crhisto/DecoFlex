
#parameteres:
root.dir <- "/mnt_volumen/GIT_REPOSITORIES/DecoFlex/DecoFlexR/vignettes/results/exp_1/"
iterations_number <- 15000
core_number <- 8

#importing functions to parallelize
library(parallel)
library(doParallel)
library(DecoFlex)

#1.oading the needed files: eset.sc.sparse_VL, deco.actual.data.VL.merge1, deco.actual.data.VL.merge2, eset.sparse_pseudobulk_VL
load(paste0(root.dir, "resources_exp1.rData"))

#2. Including generics functions.
source("/mnt_volumen/GIT_REPOSITORIES/DecoFlex/DecoFlexR/vignettes/background_scripts/generic_scripts_experiments.R")

celltypes.exp1.merge_1_2 <- c("OPCs", "Endothelial", "Neurons", "Oligodendrocytes", "Astrocytes", "Microglia", 'Neurons_Exc', 'Neurons_Inh')
celltypes.exp1.merge1 <- c("OPCs", "Endothelial", "Neurons", "Oligodendrocytes", "Astrocytes", "Microglia")

print('Executing parallelized function...')

if(is.null(core_number)){
  # Calculate the number of cores
  no_cores <- detectCores() - 1
}else{
  no_cores <- core_number
}

print(paste0('Number of cores to use: ', no_cores))

registerDoParallel(no_cores)

drop.exp1.complete <- NULL

# For each cell-type dropped in the merge1 group, I will run the deconvolution
# Definition of the parallel function
drop.exp1.complete <- foreach(j = 1:length(celltypes.exp1.merge_1_2),
                              .combine = c)  %dopar%
  {
    
    celltype_name <- celltypes.exp1.merge_1_2[j]
    
    drop.exp1.complete.temp <- NULL
    
    #1. Getting the actual drop exp0 model
    current.drop <- exp_0_results$drop.exp0.complete[[celltype_name]]
    
    if(celltype_name %in% celltypes.exp1.merge1){
      
      print(paste0('Running Decoflex deleting merge1: ', celltype_name))
      
      # 2. Creation of the semi-reference objects
      semi_reference_objects <- create_semi_reference_objects(
        extra_unknown_celltypes = 1,
        cell_type_names = rownames(current.drop$result_deco_top_cluster$h),
        bulk.data_mixtures.brain = current.drop$result_deco_top_cluster$x,
        w_fixed = current.drop$result_deco_top_cluster$w,
        marker_genes = rownames(current.drop$result_deco_top_cluster$w),
        version = 'one')
      
      # 3. Deconvolution with the partially fixed h and w
      drop.exp1.complete.temp <- run_deconvolution_decoflex_semi_reference(bulk.data_mixtures.brain = current.drop$result_deco_top_cluster$x, 
                                                                           partial_w_fixed = semi_reference_objects$partial_w_fixed, 
                                                                           partial_h_fixed = semi_reference_objects$partial_h_fixed, 
                                                                           mask_w = semi_reference_objects$mask_w,
                                                                           mask_h = semi_reference_objects$mask_h,
                                                                           number_cell_types = nrow(current.drop$result_deco_top_cluster$h), 
                                                                           extra_unknown_celltypes = 1 ,
                                                                           proportion_constraint_h = TRUE,
                                                                           max_iterations = iterations_number,
                                                                           delta_threshold = 1e-15)
      
      
    }else if(celltype_name == 'Neurons_Exc' | celltype_name == 'Neurons_Inh'){
      
      print(paste0('Running Decoflex deleting merge2: ', celltype_name))
    
      
      # 2. Creation of the semi-reference objects
      semi_reference_objects <- create_semi_reference_objects(
        extra_unknown_celltypes = 1,
        cell_type_names = rownames(current.drop$result_deco_top_cluster$h),
        bulk.data_mixtures.brain = current.drop$result_deco_top_cluster$x,
        w_fixed = current.drop$result_deco_top_cluster$w,
        marker_genes = rownames(current.drop$result_deco_top_cluster$w),
        version = 'one')
      
      # 3. Deconvolution with the partially fixed h and w
      drop.exp1.complete.temp <- run_deconvolution_decoflex_semi_reference(bulk.data_mixtures.brain = current.drop$result_deco_top_cluster$x, 
                                                                           partial_w_fixed = semi_reference_objects$partial_w_fixed, 
                                                                           partial_h_fixed = semi_reference_objects$partial_h_fixed, 
                                                                           mask_w = semi_reference_objects$mask_w,
                                                                           mask_h = semi_reference_objects$mask_h,
                                                                           number_cell_types = nrow(current.drop$result_deco_top_cluster$h), 
                                                                           extra_unknown_celltypes = 1 ,
                                                                           proportion_constraint_h = TRUE,
                                                                           max_iterations = iterations_number,
                                                                           delta_threshold = 1e-15)
    }
    
    #Saving the object with the list of results
    save(drop.exp1.complete.temp, file = paste0(root.dir, celltype_name, "_Dropped_sub.exp1_Estimates.rda"))
    
    list_values_model <- list()
    list_values_model[[celltype_name]] <- drop.exp1.complete.temp
    
    
    assign(celltype_name, list_values_model)
    
    #Returning the object
    get0(celltype_name, ifnotfound = paste0('Object not found: ', celltype_name)) 
  }
stopImplicitCluster()

#Saving everything
save(drop.exp1.complete, file = paste0(root.dir, "all_Dropped_sub.exp1_Estimates.rda"))

