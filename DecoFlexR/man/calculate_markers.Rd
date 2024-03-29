% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_guided_deconvolution.R
\name{calculate_markers}
\alias{calculate_markers}
\title{Calculate Marker Genes}
\usage{
calculate_markers(
  single_cell_data_exp,
  reference,
  group_clusters_var,
  use_min_cor_strategy = FALSE,
  delete_shared_level_markers = TRUE,
  shared_level_markers = NULL,
  delete_shared_internal_markers = TRUE,
  p_value_attribute = "p_val_adj",
  param.logfc.threshold = 2,
  param.p_val_adj = 0.05,
  test.use.value = "wilcox",
  filter_markers = NULL,
  marker_strategy = "keep_default_values",
  minimum_markers = 4,
  min_delta_cor_threshold = 0.05,
  verbose = TRUE
)
}
\arguments{
\item{single_cell_data_exp}{A data frame or matrix containing the single cell
data.}

\item{reference}{A data frame or matrix containing the reference data.}

\item{group_clusters_var}{The variable name of the group clusters in the
\code{single_cell_data_exp}.}

\item{use_min_cor_strategy}{Logical value indicating whether to use the
minimum correlation strategy. This strategy is used to identify marker genes
that are least correlated with each other, thereby ensuring uniqueness.
Default is FALSE.}

\item{delete_shared_level_markers}{Logical value indicating whether to
delete shared level markers. This parameter enables the removal of marker
genes shared across different levels to ensure uniqueness within each level.
Default is TRUE.}

\item{shared_level_markers}{A vector containing shared level markers.
NULL by default.}

\item{delete_shared_internal_markers}{Logical value indicating whether to
delete shared internal markers. This parameter allows for the removal of
marker genes that are common within a given level, to promote marker
uniqueness. Default is TRUE.}

\item{p_value_attribute}{String indicating the attribute name of the p-value
in the single cell data. Default is 'p_val_adj'.}

\item{param.logfc.threshold}{Threshold for the log-fold change, determining
the minimum acceptable change for a gene to be considered a marker.
Default is 2.0.}

\item{param.p_val_adj}{Adjusted p-value threshold. This threshold is used to
control the false discovery rate. Default is 0.05.}

\item{test.use.value}{String indicating the statistical test used for
identifying markers. Default is 'wilcox'.}

\item{filter_markers}{A vector containing markers to be filtered out before
the analysis. NULL by default.}

\item{marker_strategy}{Strategy used to select marker genes. This could be a
specific algorithm or statistical method. NULL by default.}

\item{minimum_markers}{The minimum number of marker genes to be selected,
default is set to 4. This can be adjusted based on the analysis
requirements.}

\item{min_delta_cor_threshold}{Minimum difference between the correlation of
t and t-1 steps. The idea is to skip genes that have small changes in
the correlation. For default it is a 5\%.}

\item{verbose}{Logical value indicating whether to display additional
information during the function's execution. Default is TRUE.}
}
\value{
A list of marker genes for each cell type.
}
\description{
The \code{calculate_markers} function identifies marker genes for specific cell
types based on given criteria and strategies. These strategies help ensure
that the marker genesidentified are more unique to each cell type, providing
better accuracy in downstream analysis.
}
\details{
The function calculates marker genes for each cell type in the provided data,
with options to apply different strategies and filters during the process.
It considers factors such as p-value, log-fold change, and uniqueness at
different levels to ensure the identified markers are most representative of
their respective cell types.
}
\seealso{
\code{\link{Matrix}}
}
