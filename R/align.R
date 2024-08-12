#' Assign main marker types to clusters in single-cell RNA-seq data
#'
#' This function takes a Seurat object and a list of markers, along with cluster assignments,
#' and assigns a main marker type to each cluster based on the expression of the markers.
#'
#' @param obj A Seurat object.
#' @param markers A data frame specifying marker genes and their corresponding
#'   cell types. It should have at least two columns: `markers` (gene names) and `type`
#'   (cell type labels).
#' @param clusters A vector of cluster assignments for each cell, defaulting to the
#'   'seurat_clusters' column in the metadata of the Seurat object.
#'
#' @return A list with two elements:
#' `value`: A numeric vector of the maximum weighted expression ratios for each cluster.
#' `mapping`: A character vector of the dominant marker type assigned to each cluster.
#'
#' @export
align <- function (obj, markers, clusters = obj@meta.data$seurat_clusters) {
  expr <- GetAssayData(obj)[markers$markers, ]
  expr <- (expr - min(expr)) / (max(expr) - min(expr))
  weights <- lapply(split(markers$markers, markers$type), function (X) {
    return(sum(expr[X, ]))
  }) %>% unlist()
  u.clust <- sort(unique(clusters))
  mapping <- sapply(u.clust, function (X) {
    s <- sapply(markers$markers, function (Y) {
      return(sum(expr[Y, clusters == X] / sum(expr[Y, ])))
    }) %>% split(markers$type) %>% lapply(mean) %>% unlist()
    s <- s * weights
    return(list(names(s)[which.max(s)], max(s) / sum(s)))
  })
  return(list(value = unlist(mapping[2, ]), mapping = unlist(mapping[1, ])))
}
