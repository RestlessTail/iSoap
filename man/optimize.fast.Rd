% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimize.fast.R
\name{optimize.fast}
\alias{optimize.fast}
\title{Optimize clustering resolution for single-cell data}
\usage{
optimize.fast(obj, markers, tolerance = 3, itr = 5, resolution = c(0.001, 1))
}
\arguments{
\item{obj}{A Seurat object.}

\item{markers}{A data frame specifying marker genes and their corresponding
cell types. It should have at least two columns: `markers` (gene names) and `type`
(cell type labels).}

\item{tolerance}{A numeric value specifying the tolerance level for evaluating the
clustering results. Higher values indicate a more lenient evaluation.}

\item{itr}{The number of iterations to perform in the optimization process.}

\item{resolution}{A numeric vector of length 2, specifying the initial search range
for the clustering resolution.}
}
\value{
A list containing the optimized clustering resolution, the evaluation value,
  and the mapping of clusters to cell types.
}
\description{
This function optimizes the clustering resolution parameter for single-cell RNA-seq
data, aiming to find a resolution that best matches the number of cell types present
in the data. It is particularly suitable when each cell type is relatively abundant
and the differences between cell types are clear.
}
\details{
The function uses an iterative approach to narrow down the optimal resolution range
based on the clustering results and the expression of marker genes. It is designed
to run efficiently, but may force the clustering of cells into distinct groups even
when the boundaries between cell types are blurry, or misclassify rare cell types
into more abundant ones.
}
