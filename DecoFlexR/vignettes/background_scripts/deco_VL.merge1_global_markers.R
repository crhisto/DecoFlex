
# Note: this execution using background jobs test the functionlity where
# DecoFlex caclulates in the to level of the hierarchy the marker genes for 
# all Celltypes using OMiC process and then in each level, the same set of 
# markers is used fintering the specific celtypes that the level has.

list_cell_types_VL.merge1 <- c("OPCs", "Endothelial", "Oligodendrocytes", "Astrocytes", "Microglia", "Neurons")

# Initially for the top celltypes: 
hierarchical_clustering_VL.merge1 <- list(
  '1'=list(tree_level=0, leaf_number=1, celltype_list = list_cell_types_VL.merge1, 
           next_level_clustering=NULL
  )
)

list_cell_types_VL.orig.celltype <- c('OPC', 'Endothelial', 'IN-VIP', 'L5/6-CC', 'L4', 'IN-PV', 'L2/3', 'Oligodendrocytes', 'AST-PP', 'Microglia', 'Neu-mat', 'IN-SV2C', 'L5/6', 'IN-SST', 'Neu-NRGN-I', 'AST-FB', 'Neu-NRGN-II')

hierarchical_clustering_VL.total.v2.1.top.names <- c('OPC', 'Endothelial', 'Oligodendrocytes', 'Microglia', 'Astrocytes', 'Neurons')
hierarchical_clustering_VL.total.v2.1 <- list(
  '1'=list(tree_level=0, leaf_number=1, celltype_list = list_cell_types_VL.orig.celltype, 
           next_level_clustering=list(
             '1'=list(tree_level=1, leaf_number=1, name='OPC', celltype_list=c('OPC'), next_level_clustering=NULL), 
             '2'=list(tree_level=1, leaf_number=2, name='Endothelial', celltype_list=c('Endothelial'), next_level_clustering=NULL),
             '3'=list(tree_level=1, leaf_number=3, name='Oligodendrocytes', celltype_list=c('Oligodendrocytes'), next_level_clustering=NULL),
             '4'=list(tree_level=1, leaf_number=4, name='Microglia', celltype_list=c('Microglia'), next_level_clustering=NULL),
             '5'=list(tree_level=1, leaf_number=5, name='Astrocytes', celltype_list=c('AST-FB', 'AST-PP'), next_level_clustering=NULL),
             '6'=list(tree_level=1, leaf_number=6, name='Neurons', celltype_list=c('L2/3', 'L4', 'L5/6', 'L5/6-CC', 'Neu-mat', 'IN-PV', 'IN-SST', 'IN-SV2C', 'IN-VIP'), next_level_clustering=list(
               '1'=list(tree_level=2, leaf_number=1, name='Neurons_Exc', celltype_list=c('L2/3', 'L4', 'L5/6', 'L5/6-CC', 'Neu-mat'), next_level_clustering=NULL), 
               '2'=list(tree_level=2, leaf_number=2, name='Neurons_Inh', celltype_list=c('IN-PV', 'IN-SST', 'IN-SV2C', 'IN-VIP'), next_level_clustering=NULL)
             )
             )
           )
  )
)


#All cell types: neurons groupped differently. Neu-mat is now gropped with the rest of Neu
hierarchical_clustering_VL.total.v4 <- list(
  '1'=list(tree_level=0, leaf_number=1, celltype_list = list_cell_types_VL.orig.celltype, 
           next_level_clustering=list(
             '1'=list(tree_level=1, leaf_number=1, name='OPC', celltype_list=c('OPC'), next_level_clustering=NULL), 
             '2'=list(tree_level=1, leaf_number=2, name='Endothelial', celltype_list=c('Endothelial'), next_level_clustering=NULL),
             '3'=list(tree_level=1, leaf_number=3, name='Oligodendrocytes', celltype_list=c('Oligodendrocytes'), next_level_clustering=NULL),
             '4'=list(tree_level=1, leaf_number=4, name='Microglia', celltype_list=c('Microglia'), next_level_clustering=NULL),
             '5'=list(tree_level=1, leaf_number=5, name='Astrocytes', celltype_list=c('AST-FB', 'AST-PP'), next_level_clustering=NULL),
             '6'=list(tree_level=1, leaf_number=6, name='Neurons', celltype_list=c('L2/3', 'L4', 'L5/6', 'L5/6-CC', 'IN-PV', 'IN-SST', 'IN-SV2C', 'IN-VIP', 'Neu-mat', 'Neu-NRGN-I', 'Neu-NRGN-II'), next_level_clustering=list(
               '1'=list(tree_level=2, leaf_number=1, name='Neurons_Exc_Inh', celltype_list=c('L2/3', 'L4', 'L5/6', 'L5/6-CC', 'IN-PV', 'IN-SST', 'IN-SV2C', 'IN-VIP'), next_level_clustering=list(
                 '1'=list(tree_level=3, leaf_number=1, name='Neurons_Exc', celltype_list=c('L2/3', 'L4', 'L5/6', 'L5/6-CC'), next_level_clustering=NULL), 
                 '2'=list(tree_level=3, leaf_number=2, name='Neurons_Inh', celltype_list=c('IN-PV', 'IN-SST', 'IN-SV2C', 'IN-VIP'), next_level_clustering=NULL)
               )),
               '2'=list(tree_level=2, leaf_number=2, name='Neurons_NRGN', celltype_list=c('Neu-mat', 'Neu-NRGN-I', 'Neu-NRGN-II'), next_level_clustering=NULL)
             )
             )
           )
  )
)


library(DecoFlex)

#Parameters for the simulation
single_cell_data_exp <- eset.sc.sparse_VL
hierarchy = hierarchical_clustering_VL.total.v4$`1`
bulk_data = data.frame(as.matrix(eset.sparse_pseudobulk_VL@assayData$exprs))
sub_clusters_var <- 'orig.celltype'
sample <- 'Disorder' #just one value: Control
use_min_cor_strategy = TRUE
delete_shared_level_markers = FALSE
delete_shared_internal_markers = FALSE
ordering_strategy = 'foldchange_pvalue'

start_time <- Sys.time()

#let's run de simulation
deco.actual.data.VL.merge1.global.markers <- 
  run_deconvolution_tree_guided_recursive(bulk_data = bulk_data,
                                          single_cell_data_exp = single_cell_data_exp, 
                                          hierarchy = hierarchy,
                                          sub_clusters_var = sub_clusters_var,
                                          sample = sample, 
                                          use_min_cor_strategy = use_min_cor_strategy, 
                                          ordering_strategy = ordering_strategy, 
                                          delete_shared_level_markers = delete_shared_level_markers, 
                                          delete_shared_internal_markers = delete_shared_internal_markers,
                                          param.logfc.threshold = 1.0,
                                          param.p_val_adj = 0.05,
                                          minimum_markers = 4,
                                          min_delta_cor_threshold = 0.0,
                                          max_iterations = 1,
                                          use_global_makers = TRUE,
                                          verbose = TRUE)
end_time <- Sys.time()
end_time - start_time