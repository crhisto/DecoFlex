
#parameteres:
root.dir <- "/mnt_volumen/GIT_REPOSITORIES/DecoFlex/DecoFlexR/vignettes/results/CA.orig.celltype/"
iterations_number <- 10000
core_number <- 3

#Including generics functions.
source("/mnt_volumen/GIT_REPOSITORIES/DecoFlex/DecoFlexR/vignettes/background_scripts/generic_scripts_experiments.R")

#importing functions to parallelize
library(parallel)
library(doParallel)


#1. First, let's calculate the OMiC data for the  VL orig.cell

list_cell_types_CA.orig.celltype <- c('Astro L1-6 FGFR3 SLC14A1', 'Exc L2 LAMP5 LTK', 'Exc L2-3 LINC00507 FREM3', 'Exc L3-4 RORB CARM1P1', 'Exc L3-5 RORB ESR1', 'Exc L4-5 RORB FOLH1B', 'Exc L4-6 FEZF2 IL26', 'Exc L4-6 RORB SEMA3E', 'Exc L5-6 FEZF2 ABO', 'Exc L5-6 FEZF2 EFTUD1P1', 'Exc L5-6 THEMIS C1QL3', 'Inh L1 SST NMBR', 'Inh L1-3 SST CALB1', 'Inh L1-4 LAMP5 LCP2', 'Inh L2-4 PVALB WFDC2', 'Inh L2-6 LAMP5 CA1', 'Oligo L1-6 OPALIN', 'OPC L1-6 PDGFRA')

# Initially for the top celltypes: .71
hierarchical_clustering_CA.orig.celltype <- list(
  '1'=list(tree_level=0, leaf_number=1, celltype_list = list_cell_types_CA.orig.celltype, 
           next_level_clustering=NULL
  )
)

#Parameters for the simulation
single_cell_data_exp <- eset.sc.sparse_CA
hierarchy = hierarchical_clustering_CA.orig.celltype$`1`
bulk_data = data.frame(as.matrix(eset.sparse_pseudobulk_CA@assayData$exprs))
sub_clusters_var <- 'orig.celltype'
sample <- 'Disorder' #just one value: Control
use_min_cor_strategy = TRUE
delete_shared_level_markers = FALSE
delete_shared_internal_markers = FALSE
ordering_strategy = 'foldchange_pvalue'


start_time <- Sys.time()

# #let's run de simulation
# marker.genes.OMiC.CA.orig.celltype <- 
#   run_marker_selection_OMiC(bulk_data = bulk_data,
#                             single_cell_data_exp = single_cell_data_exp, 
#                             hierarchy = hierarchy,
#                             sub_clusters_var = sub_clusters_var,
#                             sample = sample, 
#                             use_min_cor_strategy = use_min_cor_strategy, 
#                             ordering_strategy = ordering_strategy, 
#                             delete_shared_level_markers = delete_shared_level_markers, 
#                             delete_shared_internal_markers = delete_shared_internal_markers,
#                             param.logfc.threshold = 0,
#                             param.p_val_adj = 1,
#                             minimum_markers = 4,
#                             min_delta_cor_threshold = 0.0,
#                             verbose = TRUE)  
# 
# data.correlation.all.CA.orig.celltype <- marker.genes.OMiC.CA.orig.celltype$accumulative_corr_markers




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
list_markers <- c('1000', '1200', '1500')
all_deco_CA_orig.celltype <- NULL

# For each cell-type dropped in the merge1 group, I will run the deconvolution
# Definition of the parallel function
all_deco_CA_orig.celltype <- foreach(j = 1:length(list_markers),
                              .combine = c)  %dopar%
  {
    
    marker_value <- list_markers[j]
    
    print(paste0('Model with fixed markers: ', marker_value))
    
    deco_CA_orig.celltype.temp <- run_deco_CA_orig.celltype(data.correlation.all.CA.orig.celltype, marker_value, max_iterations = iterations_number)
    
    name_model <- paste0('markers_', marker_value)
    
    #Save the model independently
    save(deco_CA_orig.celltype.temp, file = paste0(root.dir, name_model, "_deco_CA_orig.celltype.temp.rDAta"))
    
    list_values_model <- list()
    list_values_model[[name_model]] <- deco_CA_orig.celltype.temp
    
    assign(name_model, list_values_model)
    
    #Returning the object
    get0(name_model, ifnotfound = paste0('Object not found: ', name_model)) 
  }
stopImplicitCluster()

#Saving everything
#save(all_deco_CA_orig.celltype, file = paste0(root.dir, "all_deco_CA_orig.celltype.rda"))
