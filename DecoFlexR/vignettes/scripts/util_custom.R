#' Plot Heatmap of Correlation Between Deconvolution Results and True Values
#'
#' This function plots a heatmap showing the correlation between deconvolution
#' results and true values. It first calculates the correlation matrix between
#' the estimated methylation profiles and reference methylation profiles. Then
#' it identifies the references that had the highest correlation with each of
#' the estimated methylation profiles. Finally, it plots a heatmap of these
#' correlations.
#'
#' @param reference_meth Methylation reference for the cell-types in the
#'  analysis. It could be  the one present in: EDecExampleData::reference_meth
#' @param reference_meth_class Class on the methylation reference. Could be
#'  the one present in EDecExampleData::reference_meth_class.
#' @param deco.results.w A matrix where rows are the overrepresented markers
#' and columns are the estimated methylation profiles. The matrix contains
#' the deconvolution results.
#' @param markers_ovr A vector containing the overrepresented markers.
#' @param k The number of clusters to consider when identifying the references
#' with the highest correlation. Defaults to 4.
#'
#' @import ggplot2
#' @import grDevices
#' @import stats
#'
#' @return This function doesn't return a value. It generates a heatmap plot
#' in the graphics device.
#'
#' @export
plot_heatmap_correlation_deconvolution_vs_true_Edec <- function(
    reference_meth,
    reference_meth_class,
    deco.results.w,
    markers_ovr,
    k = 4){

  # Defining the w matrix
  rownames(deco.results.w) <- markers_ovr

  # Let's see the correlation of X with the x_hat calculated in the
  # deconvolution. Compute correlation between estimated methylation profiles,
  # and reference methylation profiles
  cors_deconv_refs_4ct.vs.python.imp <-
    stats::cor(reference_meth[markers_ovr,],
               deco.results.w[markers_ovr,])

  # Check what references had the highest correlation with each
  # of the estimated methylation profiles
  best_cors = rbind(apply(cors_deconv_refs_4ct.vs.python.imp,2,which.max),
                    apply(cors_deconv_refs_4ct.vs.python.imp,2,max))

  best_cor_labels = matrix("",nrow=nrow(cors_deconv_refs_4ct.vs.python.imp),
                           ncol=ncol(cors_deconv_refs_4ct.vs.python.imp))
  for (i in 1:k){
    best_cor_labels[best_cors[1,i],i] = as.character(round(best_cors[2,i],2))
  }

  # Create a vector of colors representing the class of each reference
  ref_class_colors <- as.factor(reference_meth_class)
  levels(ref_class_colors) <- RColorBrewer::brewer.pal(4,"Accent")
  ref_class_colors <- as.character(ref_class_colors)

  # Create a color gradient to be used in a heatmap of correlations
  color_gradient <- grDevices::colorRampPalette(c("white","steelblue"))

  # Plot correlation matrix
  gplots::heatmap.2(cors_deconv_refs_4ct.vs.python.imp,
                    trace="none",
                    col=color_gradient(10),
                    breaks=seq(0,1,0.1),
                    margins=c(4,12),
                    RowSideColors = ref_class_colors,
                    cellnote = best_cor_labels,
                    notecol="black")
}

#' @title Plot Estimated Profiles vs Truly Mixtures
#'
#' @description Compute and plot correlation between estimated methylation
#' profiles
#' and methylation profiles truly used to build the mixtures
#'
#' @param true_cell_type_meth A matrix or dataframe, contains the true
#'  methylation profiles of cell types.
#' @param w.matrix A matrix or dataframe, contains the estimated methylation
#'  profiles.
#' @param markers A vector of character, names of the markers for methylation.
#'
#' @details This function first computes the correlation between true and
#' estimated methylation profiles, and then generates a heatmap to visualize
#' the correlation. A color gradient is used in the heatmap to indicate the
#' strength of the correlation. A cell note is added to each cell of the
#' heatmap to show the round-off correlation value.
#'
#' @return A heatmap plot showing the correlation between estimated methylation
#' profiles and methylation profiles truly used to build the mixtures.
#'
plot_Edec.estimated.profiles.vs.truly.mixtures <- function(
    true_cell_type_meth,
    w.matrix,
    markers){
  # Compute and plot correlation between estimated methylation profiles,
  # and methylation profiles truly used to build the mixtures
  cors_deconv_true = cor(true_cell_type_meth[markers,],
                         w.matrix[markers,])

  # Create a color gradient to be used in a heatmap of correlations
  color_gradient <- colorRampPalette(c("white","steelblue"))
  #color_gradient <- colorRampPalette(brewer.pal(10, "RdYlBu"))

  gplots::heatmap.2(cors_deconv_true,
                    cellnote = round(cors_deconv_true, 2),
                    trace="none",
                    col=color_gradient(10),
                    breaks=seq(0,1,0.1),
                    margins=c(12,12),
                    cexRow =1, cexCol = 1)

}

