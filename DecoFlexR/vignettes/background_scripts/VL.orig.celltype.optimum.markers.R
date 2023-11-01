
#parameteres:
root.dir <- "/mnt_volumen/GIT_REPOSITORIES/DecoFlex/DecoFlexR/vignettes/results/VL.orig.celltype/"
iterations_number <- 15000
core_number <-5
list_markers <- c('2000', '3000', '4000', '5000', '6000')

#Including generics functions.
source("/mnt_volumen/GIT_REPOSITORIES/DecoFlex/DecoFlexR/vignettes/background_scripts/generic_scripts_experiments.R")

#load needed objects to run the script.
load(paste0(root.dir, "resources_deco_VL_orig.celltype.rDAta"))

#importing functions to parallelize
library(parallel)
library(doParallel)


#1. First, let's calculate the OMiC data for the  VL orig.cell

list_cell_types_VL.orig.celltype <- c('OPC', 'Endothelial', 'IN-VIP', 'L5/6-CC', 'L4', 'IN-PV', 'L2/3', 'Oligodendrocytes', 'AST-PP', 'Microglia', 'Neu-mat', 'IN-SV2C', 'L5/6', 'IN-SST', 'Neu-NRGN-I', 'AST-FB', 'Neu-NRGN-II')

# Initially for the top celltypes: .71
hierarchical_clustering_VL.orig.celltype <- list(
  '1'=list(tree_level=0, leaf_number=1, celltype_list = list_cell_types_VL.orig.celltype, 
           next_level_clustering=NULL
  )
)

#Parameters for the simulation
single_cell_data_exp <- eset.sc.sparse_VL
hierarchy = hierarchical_clustering_VL.orig.celltype$`1`
bulk_data = data.frame(as.matrix(eset.sparse_pseudobulk_VL@assayData$exprs))
sub_clusters_var <- 'orig.celltype'
sample <- 'Disorder' #just one value: Control
use_min_cor_strategy = TRUE
delete_shared_level_markers = FALSE
delete_shared_internal_markers = FALSE
ordering_strategy = 'foldchange_pvalue'

run_markers <- TRUE
if(run_markers){
  
  start_time <- Sys.time()
  # #let's run de simulation
  marker.genes.OMiC.VL.orig.celltype <-
    run_marker_selection_OMiC(bulk_data = bulk_data,
                              single_cell_data_exp = single_cell_data_exp,
                              hierarchy = hierarchy,
                              sub_clusters_var = sub_clusters_var,
                              sample = sample,
                              use_min_cor_strategy = use_min_cor_strategy,
                              ordering_strategy = ordering_strategy,
                              delete_shared_level_markers = delete_shared_level_markers,
                              delete_shared_internal_markers = delete_shared_internal_markers,
                              param.logfc.threshold = 0,
                              param.p_val_adj = 0.1,
                              minimum_markers = 4,
                              min_delta_cor_threshold = 0.0,
                              verbose = TRUE)
    
}

data.correlation.all.VL.orig.celltype <- marker.genes.OMiC.VL.orig.celltype$accumulative_corr_markers




#2. I will use the marker optimum object to 

print('Executing parallelized function...')

if(is.null(core_number)){
  # Calculate the number of cores
  no_cores <- detectCores() - 1
}else{
  no_cores <- core_number
}

print(paste0('Number of cores to use: ', no_cores))

registerDoParallel(no_cores)

#Runnning DecoFlex with different values of markers to see behaviour of R , iterations and convergence. Note: 4890=272 batch
all_deco_VL_orig.celltype <- NULL

# For each cell-type dropped in the merge1 group, I will run the deconvolution
# Definition of the parallel function
all_deco_VL_orig.celltype <- foreach(j = 1:length(list_markers),
                              .combine = c)  %dopar%
  {
    
    marker_value <- list_markers[j]
    
    print(paste0('Model with fixed markers: ', marker_value))
    
    deco_VL_orig.celltype.temp <- run_deco_VL(data.correlation.all.VL = data.correlation.all.VL.orig.celltype, marker_value = marker_value, 
                                              max_iterations = iterations_number,
                                              hierarchical_clustering_VL = hierarchical_clustering_VL.orig.celltype, 
                                              sub_clusters_var = 'orig.celltype', param.logfc.threshold = 0.1)
    
    name_model <- paste0('markers_', marker_value)
    
    #Save the model independently
    save(deco_VL_orig.celltype.temp, file = paste0(root.dir, name_model, "_deco_VL_orig.celltype.temp.rDAta"))
    
    list_values_model <- list()
    list_values_model[[name_model]] <- deco_VL_orig.celltype.temp
    
    assign(name_model, list_values_model)
    
    #Returning the object
    get0(name_model, ifnotfound = paste0('Object not found: ', name_model)) 
  }
stopImplicitCluster()

#Saving everything
#save(all_deco_VL_orig.celltype, file = paste0(root.dir, "all_deco_VL_orig.celltype.rda"))
