% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/automatic_tree_discovery.R
\name{plot_hierarchical_clustering}
\alias{plot_hierarchical_clustering}
\title{Plot hierarchical clustering}
\usage{
plot_hierarchical_clustering(
  hierar_clust,
  k_value,
  fit_celltypes,
  title = NULL
)
}
\arguments{
\item{hierar_clust}{A hierarchical clustering object, which is the result of
hclust function or similar.}

\item{k_value}{An integer value indicating the number of clusters to
highlight in the plot.}

\item{fit_celltypes}{A vector containing the assignment of each data point
to a cluster.}

\item{title}{Optional; a character string representing the title of the plot.
If NULL, no title is added.}
}
\value{
This function generates a plot as a side-effect. The returned value
is a table summarizing the cluster assignments.
}
\description{
This function generates a plot of the hierarchical clustering and includes
the option for a custom title. The plot visualizes the dendrogram resulting
from the clustering, with the height of the cut denoted by a green line.
Clusters resulting from the cut are highlighted with green borders.
The function also outputs a table summarizing the resulting cluster
assignments.
}
