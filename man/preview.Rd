% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preview.R
\name{preview}
\alias{preview}
\title{Preview cell clusters in a reduced dimensional space}
\usage{
preview(obj, opt, reduction = "umap", clusters = obj@meta.data$seurat_clusters)
}
\arguments{
\item{obj}{A Seurat object.}

\item{opt}{A list containing a `mapping` element. This parameter is
expected to be the return value of one of the optimization functions
in this package, such as \code{optimize.fast}, but it only requires
the `mapping` element to be present.}

\item{reduction}{The name of the reduction method to use for extracting cell
coordinates, defaults to "umap".}

\item{clusters}{A vector of cell cluster IDs, defaults to the 'seurat_clusters'
column in the object's meta.data.}
}
\value{
A list containing two elements:
  `cells`: A data frame of the original cell coordinates in the specified
    reduced dimensional space, with an additional 'type' column indicating the
    cell type.
  `labs`: A data frame of centroid positions for each cell type in the
    specified reduced dimensional space, with columns 'x' and 'y' for the
    coordinates and a 'type' column indicating the cell type.
}
\description{
This function generates a preview of cell clusters in a specified reduced dimensional
space (e.g., UMAP or tSNE) by calculating the centroid positions of each cell type
and returning both the centroid positions and the original cell coordinates.
}
