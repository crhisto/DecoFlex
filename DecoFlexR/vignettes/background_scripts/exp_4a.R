# The idea is to increment the celltypes deleted, therefore the references are going to be delete from the celltype with larger size to the higher one. (contraty to exp_4)

#parameteres:
root.dir <- "/mnt_volumen/GIT_REPOSITORIES/DecoFlex/DecoFlexR/vignettes/results/exp_4a/"
iterations_number <- 15000
core_number <- 6

#importing functions to parallelize
library(parallel)
library(doParallel)
library(DecoFlex)

#Loading the needed files: eset.sc.sparse_VL, deco.actual.data.VL.merge1, deco.actual.data.VL.merge2, eset.sparse_pseudobulk_VL
load(paste0(root.dir, "resources_exp4.rData"))

# 1. I need to create the bulk data with the additional information
#creation of pseudo bulk data with the same information for merge1 and merge2
pseudo.eset.sc.sparse_VL.all.merge1 <- generateBulk_allcells(eset.sc.sparse_VL, ct.varname = "merge1", sample = "Disorder", ct.sub = NULL)
pseudo.eset.sc.sparse_VL.all.merge2 <- generateBulk_allcells(eset.sc.sparse_VL, ct.varname = "merge2", sample = "Disorder", ct.sub = NULL)


#2. Including generics functions.
source("/mnt_volumen/GIT_REPOSITORIES/DecoFlex/DecoFlexR/vignettes/background_scripts/generic_scripts_experiments.R")

drop.exp4a.complete <- drop.sub.exp4a.complete <- list()

#Order based on celltype size
celltypes.exp4a.merge1 <- c("Neurons", "Oligodendrocytes", "Astrocytes", "Endothelial", "Microglia", "OPCs")

celltypes.exp4a.merge1.obj <- list()
celltypes.exp4a.merge1.obj$one <- c("Neurons")
celltypes.exp4a.merge1.obj$two <- c("Neurons", "Oligodendrocytes")
celltypes.exp4a.merge1.obj$three <- c("Neurons", "Oligodendrocytes", "Astrocytes")
celltypes.exp4a.merge1.obj$four <- c("Neurons", "Oligodendrocytes", "Astrocytes", "Endothelial")
celltypes.exp4a.merge1.obj$five <- c("Neurons", "Oligodendrocytes", "Astrocytes", "Endothelial", "Microglia")
celltypes.exp4a.merge1.obj$six <- c("Neurons", "Oligodendrocytes", "Astrocytes", "Endothelial", "Microglia", "OPCs")

print('Executing parallelized function...')

if(is.null(core_number)){
  # Calculate the number of cores
  no_cores <- detectCores() - 1
}else{
  no_cores <- core_number
}

print(paste0('Number of cores to use: ', no_cores))

registerDoParallel(no_cores)

drop.exp4a.complete <- NULL

# For each cell-type dropped in the merge1 group, I will run the deconvolution
# Definition of the parallel function
drop.exp4a.complete <- foreach(j = 1:length(celltypes.exp4a.merge1.obj),
                   .combine = c)  %dopar%
  {
    
    celltype_name <- names(celltypes.exp4a.merge1.obj)[j]
    celltype_list_missing <- celltypes.exp4a.merge1.obj[[celltype_name]]

    
    
    drop.exp4a.complete.temp <- NULL
    
    print(paste0('Running Decoflex deleting merge1: ', paste(shQuote(celltype_list_missing), collapse=", ")))
    
    # 1. Filtering the deleted celltype
    rest_celltypes <- celltypes.exp4a.merge1[!(celltypes.exp4a.merge1 %in% celltype_list_missing)]
    
    # 2. Creation of bulk data with extra proportions
    #I get the real proportions for all celltypes even the one that I don't know. The unknown are at the end.
    pseudo.eset.sc.sparse_VL.all.prop <- data.frame(Control=pseudo.eset.sc.sparse_VL.all.merge1$truep[c(rest_celltypes, celltype_list_missing)])
    rownames(pseudo.eset.sc.sparse_VL.all.prop)[which(rownames(pseudo.eset.sc.sparse_VL.all.prop) %in% celltype_list_missing)] <- paste0('unknown_', 1:length(celltype_list_missing))
    
    #2.1 Fixing proportions
    eset.sparse_pseudobulk_VL.plus.fixed <- data.frame(round((pseudo.eset.sc.sparse_VL.all.merge1$pseudo_eset@assayData$exprs)/100),
                                                       eset.sparse_pseudobulk_VL@assayData$exprs)
    
    # 3. Creation of the semi-reference objects
    semi_reference_objects.exp4a.temp <- create_semi_reference_objects(
      extra_unknown_celltypes = length(celltype_list_missing), #this will increased from 1 to 6 (all celltypes)
      cell_type_names = rest_celltypes,
      bulk.data_mixtures.brain = eset.sparse_pseudobulk_VL.plus.fixed,
      w_fixed.value = subset(data.frame(deco.actual.data.VL.merge1$result_deco_top_cluster$w), select=rest_celltypes),
      marker_genes = rownames(deco.actual.data.VL.merge1$result_deco_top_cluster$w),
      extra_marker_genes_semireference.value = NULL,
      fixed_h_values = pseudo.eset.sc.sparse_VL.all.prop, 
      version = 'two')
    
    
    #3.1 Reordering the data.
    gene_names_order <- rownames(deco.actual.data.VL.merge1$result_deco_top_cluster$w)
    semi_reference_objects.exp4a.temp$partial_w_fixed <- data.frame(semi_reference_objects.exp4a.temp$partial_w_fixed)
    semi_reference_objects.exp4a.temp$partial_w_fixed[gene_names_order, rest_celltypes] <- data.frame(deco.actual.data.VL.merge1$result_deco_top_cluster$w[gene_names_order , rest_celltypes]) 
    
    # 4. Deconvolution with the partially fixed h and w
    drop.exp4a.complete.temp <- run_deconvolution_decoflex_semi_reference(
      bulk.data_mixtures.brain = semi_reference_objects.exp4a.temp$bulk_data,
      partial_w_fixed = semi_reference_objects.exp4a.temp$partial_w_fixed, 
      partial_h_fixed = semi_reference_objects.exp4a.temp$partial_h_fixed, 
      mask_w = semi_reference_objects.exp4a.temp$mask_w,
      mask_h = semi_reference_objects.exp4a.temp$mask_h,
      scale_w_unfixed_col = TRUE,
      number_cell_types = length(rest_celltypes), 
      extra_unknown_celltypes = length(celltype_list_missing) ,
      proportion_constraint_h = TRUE,
      max_iterations = iterations_number,
      delta_threshold = 1e-15)
    

    #Saving the object with the list of results
    save(drop.exp4a.complete.temp, file = paste0(root.dir, celltype_name, "_Dropped_sub.exp4a_Estimates.rda"))
    
    list_values_model <- list()
    list_values_model[[celltype_name]] <- drop.exp4a.complete.temp
    
    
    assign(celltype_name, list_values_model)
    
    #Returning the object
    get0(celltype_name, ifnotfound = paste0('Object not found: ', celltype_name)) 
  }
stopImplicitCluster()

#Saving everything
save(drop.exp4a.complete, file = paste0(root.dir, "all_Dropped_sub.exp4a_Estimates.rda"))

