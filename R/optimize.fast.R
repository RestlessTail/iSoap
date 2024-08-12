#' Optimize clustering resolution for single-cell data
#'
#' This function optimizes the clustering resolution parameter for single-cell RNA-seq
#' data, aiming to find a resolution that best matches the number of cell types present
#' in the data. It is particularly suitable when each cell type is relatively abundant
#' and the differences between cell types are clear.
#'
#' The function uses an iterative approach to narrow down the optimal resolution range
#' based on the clustering results and the expression of marker genes. It is designed
#' to run efficiently, but may force the clustering of cells into distinct groups even
#' when the boundaries between cell types are blurry, or misclassify rare cell types
#' into more abundant ones.
#'
#' @param obj A Seurat object.
#' @param markers A data frame specifying marker genes and their corresponding
#'   cell types. It should have at least two columns: `markers` (gene names) and `type`
#'   (cell type labels).
#' @param tolerance A numeric value specifying the tolerance level for evaluating the
#'   clustering results. Higher values indicate a more lenient evaluation.
#' @param itr The number of iterations to perform in the optimization process.
#' @param resolution A numeric vector of length 2, specifying the initial search range
#'   for the clustering resolution.
#'
#' @return A list containing the optimized clustering resolution, the evaluation value,
#'   and the mapping of clusters to cell types.
#'
#' @export
optimize.fast <- function (obj, markers, tolerance = 3, itr = 5, resolution = c(0.001, 1)) {
  expr <- GetAssayData(obj)[markers$markers, ]
  expr <- (expr - min(expr)) / (max(expr) - min(expr))
  weights <- lapply(split(markers$markers, markers$type), function (X) {
    return(sum(expr[X, ]))
  }) %>% unlist()

  x.l <- resolution[1]
  x.r <- resolution[2]
  x.m <- exp(mean(log(c(x.l, x.r))))
  y.l <- eval.fast(obj, markers, weights, tolerance, expr, x.l)
  y.r <- eval.fast(obj, markers, weights, tolerance, expr, x.r)
  y.m <- eval.fast(obj, markers, weights, tolerance, expr, x.m)
  if (y.l$value * y.r$value > 0) { stop("Too narrow range of resolution.") }
  for (i in 1:(itr - 1)) {
    if (y.l$value * y.m$value < 0) {
      x.r <- x.m
    } else {
      x.l <- x.m
    }
    x.m <- exp(mean(log(c(x.l, x.r))))
    y.m <- eval.fast(obj, markers, weights, tolerance, expr, x.m)
  }
  return(y.m)
}

eval.fast <- function (obj, markers, weights, tolerance, expr, resolution) {
  obj %<>% FindClusters(resolution = resolution)
  clusters <- obj@meta.data$seurat_clusters
  u.clust <- sort(unique(clusters))
  scores <- sapply(u.clust, function (X) {
    s <- sapply(markers$markers, function (Y) {
      return(sum(expr[Y, clusters == X] / sum(expr[Y, ])))
    }) %>% split(markers$type) %>% lapply(mean) %>% unlist()
    s <- s * weights
    return(list(names(s)[which.max(s)], max(s) / sum(s)))
  })
  d <- length(unique(clusters)) - length(unique(unlist(scores[1, ])))
  return(list(resolution = resolution, value = ifelse(d >= 0, 1, -1) * exp(abs(d) - tolerance) / exp(mean(log(unlist(scores[2, ])))), mapping = unlist(scores[1, ])))
}
