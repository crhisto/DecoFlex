# The main idea is to make empty the genes that below to the missing celltype but 
# in the known celltypes, mimicking the missing values that we could have
# in a reference datasets, even for the known celltypes. For this I will 
# check for each missed celltype the marker genes and then I will keep those 
# unknown for both the known and unknown celltypes.

#parameteres:
root.dir <- "/mnt_volumen/GIT_REPOSITORIES/DecoFlex/DecoFlexR/vignettes/results/exp_3c/"
iterations_number <- 15000
core_number <- 8

#importing functions to parallelize
library(parallel)
library(doParallel)
library(DecoFlex)
library(reticulate)

#Loading the needed files: eset.sc.sparse_VL, deco.actual.data.VL.merge1, deco.actual.data.VL.merge2, eset.sparse_pseudobulk_VL
load(paste0(root.dir, "resources_exp3.rData"))
load(paste0(root.dir, "extra_resources_exp3c.rda"))

# 1. I need to create the bulk data with the additional information
#creation of pseudo bulk data with the same information for merge1 and merge2
pseudo.eset.sc.sparse_VL.all.merge1 <- generateBulk_allcells(eset.sc.sparse_VL, ct.varname = "merge1", sample = "Disorder", ct.sub = NULL)
pseudo.eset.sc.sparse_VL.all.merge2 <- generateBulk_allcells(eset.sc.sparse_VL, ct.varname = "merge2", sample = "Disorder", ct.sub = NULL)


#2. Including generics functions.
source("/mnt_volumen/GIT_REPOSITORIES/DecoFlex/DecoFlexR/vignettes/background_scripts/generic_scripts_experiments.R")

drop.exp3c.complete <- drop.sub.exp3c.complete <- list()

celltypes.exp3c.merge_1_2 <- c("OPCs", "Endothelial", "Neurons", "Oligodendrocytes", "Astrocytes", "Microglia", 'Neurons_Exc', 'Neurons_Inh')

celltypes.exp3c.merge1 <- c("OPCs", "Endothelial", "Neurons", "Oligodendrocytes", "Astrocytes", "Microglia")
celltypes.exp3c.merge2 <- c("OPCs", "Endothelial", "Oligodendrocytes", "Astrocytes", "Microglia", 'Neurons_Exc', 'Neurons_Inh')

print('Executing parallelized function...')

if(is.null(core_number)){
  # Calculate the number of cores
  no_cores <- detectCores() - 1
}else{
  no_cores <- core_number
}

print(paste0('Number of cores to use: ', no_cores))

registerDoParallel(no_cores)

drop.exp3c.complete <- NULL

# For each cell-type dropped in the merge1 group, I will run the deconvolution
# Definition of the parallel function
drop.exp3c.complete <- foreach(j = 1:length(celltypes.exp3c.merge_1_2),
                   .combine = c)  %dopar%
  {
    
    celltype_name <- celltypes.exp3c.merge_1_2[j]

    drop.exp3c.complete.temp <- NULL
    
    if(celltype_name %in% celltypes.exp3c.merge1){
      
      print(paste0('Running Decoflex deleting merge1: ', celltype_name))
      
      # 1. Filtering the deleted celltype
      rest_celltypes <- celltypes.exp3c.merge1[!celltypes.exp3c.merge1 == celltype_name]
      
      
      # 2. Creation of bulk data with extra proportions
      #I get the real proportions for all celltypes even the one that I don't know.
      pseudo.eset.sc.sparse_VL.all.prop <- data.frame(Control=pseudo.eset.sc.sparse_VL.all.merge1$truep[c(rest_celltypes, celltype_name)])
      rownames(pseudo.eset.sc.sparse_VL.all.prop)[which(rownames(pseudo.eset.sc.sparse_VL.all.prop) == celltype_name)] <- 'unknown_1'
      
      #2.1 Fixing proportions
      eset.sparse_pseudobulk_VL.plus.fixed <- data.frame(round((pseudo.eset.sc.sparse_VL.all.merge1$pseudo_eset@assayData$exprs)/100),
                                                         eset.sparse_pseudobulk_VL@assayData$exprs)
      
      # 3. Getting the celltype markers
      markers_celltype_dropped <- markers.list.VL.merge1[[celltype_name]]
      all_markers <- rownames(deco.actual.data.VL.merge1$result_deco_top_cluster$w)
      rest_markers <- setdiff(all_markers, markers_celltype_dropped)
      
      # 4. Creation of the semi-reference objects. I will send just the rest_markers
      # to be able to construct the semi-reference for markers from the dropped
      # celltype.
      semi_reference_objects.exp3c.temp <- create_semi_reference_objects(
        extra_unknown_celltypes = 1,
        cell_type_names = rest_celltypes,
        bulk.data_mixtures.brain = eset.sparse_pseudobulk_VL.plus.fixed,
        w_fixed.value = data.frame(deco.actual.data.VL.merge1$result_deco_top_cluster$w[rest_markers, rest_celltypes]), 
        marker_genes = rest_markers,
        extra_marker_genes_semireference.value = markers_celltype_dropped,
        fixed_h_values = pseudo.eset.sc.sparse_VL.all.prop, 
        version = 'two')
      
      
      #4.1 Reordering the data.
      gene_names_order <- c(rest_markers, markers_celltype_dropped)
      semi_reference_objects.exp3c.temp$partial_w_fixed <- data.frame(semi_reference_objects.exp3c.temp$partial_w_fixed)
      semi_reference_objects.exp3c.temp$partial_w_fixed[rest_markers, rest_celltypes] <- data.frame(deco.actual.data.VL.merge1$result_deco_top_cluster$w[rest_markers , rest_celltypes]) 

      # 5. Deconvolution with the partially fixed h and w
      drop.exp3c.complete.temp <- run_deconvolution_decoflex_semi_reference(
        bulk.data_mixtures.brain = semi_reference_objects.exp3c.temp$bulk_data,
        partial_w_fixed = semi_reference_objects.exp3c.temp$partial_w_fixed, 
        partial_h_fixed = semi_reference_objects.exp3c.temp$partial_h_fixed, 
        mask_w = semi_reference_objects.exp3c.temp$mask_w,
        mask_h = semi_reference_objects.exp3c.temp$mask_h,
        scale_w_unfixed_col = TRUE,
        number_cell_types = length(rest_celltypes), 
        extra_unknown_celltypes = 1 ,
        proportion_constraint_h = TRUE,
        max_iterations = iterations_number,
        delta_threshold = 1e-15)
      

    }else if(celltype_name == 'Neurons_Exc' | celltype_name == 'Neurons_Inh'){

      print(paste0('Running Decoflex deleting merge2: ', celltype_name))
      
      # 1. Filtering the deleted celltype
      rest_celltypes <- celltypes.exp3c.merge2[!celltypes.exp3c.merge2 == celltype_name]
      
      
      # 2. Creation of bulk data with extra proportions
      #I get the real proportions for all celltypes even the one that I don't know.
      pseudo.eset.sc.sparse_VL.all.prop <- data.frame(Control=pseudo.eset.sc.sparse_VL.all.merge2$truep[c(rest_celltypes, celltype_name)])
      rownames(pseudo.eset.sc.sparse_VL.all.prop)[which(rownames(pseudo.eset.sc.sparse_VL.all.prop) == celltype_name)] <- 'unknown_1'
      
      #2.1 Fixing proportions
      eset.sparse_pseudobulk_VL.plus.fixed <- data.frame(round((pseudo.eset.sc.sparse_VL.all.merge2$pseudo_eset@assayData$exprs)/100),
                                                         eset.sparse_pseudobulk_VL@assayData$exprs)
      
      # 3. Getting the celltype markers
      markers_celltype_dropped <- markers.list.VL.merge2[[celltype_name]]
      all_markers <- rownames(deco.actual.data.VL.merge2$result_deco_top_cluster$w)
      rest_markers <- setdiff(all_markers, markers_celltype_dropped)
      
      
      # 4. Creation of the semi-reference objects. I will send just the rest_markers
      # to be able to construct the semi-reference for markers from the dropped
      # celltype.
      semi_reference_objects.exp3c.temp <- create_semi_reference_objects(
        extra_unknown_celltypes = 1,
        cell_type_names = rest_celltypes,
        bulk.data_mixtures.brain = eset.sparse_pseudobulk_VL.plus.fixed,
        w_fixed.value = data.frame(deco.actual.data.VL.merge2$result_deco_top_cluster$w[rest_markers, rest_celltypes]), 
        marker_genes = rest_markers,
        extra_marker_genes_semireference.value = markers_celltype_dropped,
        fixed_h_values = pseudo.eset.sc.sparse_VL.all.prop, 
        version = 'two')
      
      #4.1 Reordering the data.
      gene_names_order <- c(rest_markers, markers_celltype_dropped)
      semi_reference_objects.exp3c.temp$partial_w_fixed <- data.frame(semi_reference_objects.exp3c.temp$partial_w_fixed)
      semi_reference_objects.exp3c.temp$partial_w_fixed[rest_markers, rest_celltypes] <- data.frame(deco.actual.data.VL.merge2$result_deco_top_cluster$w[rest_markers , rest_celltypes]) 
      

      # 5. Deconvolution with the partially fixed h and w
      drop.exp3c.complete.temp <- run_deconvolution_decoflex_semi_reference(
        bulk.data_mixtures.brain = semi_reference_objects.exp3c.temp$bulk_data,
        partial_w_fixed = semi_reference_objects.exp3c.temp$partial_w_fixed, 
        partial_h_fixed = semi_reference_objects.exp3c.temp$partial_h_fixed, 
        mask_w = semi_reference_objects.exp3c.temp$mask_w,
        mask_h = semi_reference_objects.exp3c.temp$mask_h,
        scale_w_unfixed_col = TRUE,
        number_cell_types = length(rest_celltypes), 
        extra_unknown_celltypes = 1 ,
        proportion_constraint_h = TRUE,
        max_iterations = iterations_number,
        delta_threshold = 1e-15)
      
    }
    
    #Saving the object with the list of results
    save(drop.exp3c.complete.temp, file = paste0(root.dir, celltype_name, "_Dropped_sub.exp3c_Estimates.rda"))
    
    list_values_model <- list()
    list_values_model[[celltype_name]] <- drop.exp3c.complete.temp
    
    
    assign(celltype_name, list_values_model)
    
    #Returning the object
    get0(celltype_name, ifnotfound = paste0('Object not found: ', celltype_name)) 
  }
stopImplicitCluster()

#Saving everything
save(drop.exp3c.complete, file = paste0(root.dir, "all_Dropped_sub.exp3c_Estimates.rda"))

