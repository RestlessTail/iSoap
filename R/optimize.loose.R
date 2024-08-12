#' Optimize Cluster Resolution for Enhanced Cell Type Identification
#'
#' `optimize.loose` automatically adjusts the resolution parameter within a Seurat object's
#' single-cell RNA-seq data to achieve more precise matching of cell types, especially for
#' rare and difficult-to-detect populations. This function employs an optimization strategy
#' that may extend beyond the initially specified resolution range to find the optimal
#' clustering configuration.
#'
#' @param obj A \code{Seurat} object.
#' @param markers A data frame specifying marker genes and their corresponding
#'   cell types. It should have at least two columns: `markers` (gene names) and `type`
#'   (cell type labels).
#' @param itr The number of iterations for the optimization algorithm. Higher values may
#'   improve the accuracy but increase runtime.
#' @param resolution A numeric vector of length 2, defining the initial minimum and maximum
#'   resolution values to consider during optimization. Note that the actual optimization
#'   process may extend beyond this range.
#'
#' @return A list containing the optimized clustering resolution, the evaluation value,
#'   and the mapping of clusters to cell types.
#'
#' @export
optimize.loose <- function (obj, markers, itr = 5, resolution = c(0.001, 1)) {
  expr <- GetAssayData(obj)[markers$markers, ]
  expr <- (expr - min(expr)) / (max(expr) - min(expr))
  weights <- lapply(split(markers$markers, markers$type), function (X) {
    return(sum(expr[X, ]))
  }) %>% unlist()

  x.l <- resolution[1]
  x.r <- resolution[2]
  x.m1 <- x.l + (x.r - x.l) / 3
  x.m2 <- x.r - (x.r - x.l) / 3
  y.l <- eval.loose(obj, markers, weights, expr, x.l)
  y.r <- eval.loose(obj, markers, weights, expr, x.r)
  y.m1 <- eval.loose(obj, markers, weights, expr, x.m1)
  y.m2 <- eval.loose(obj, markers, weights, expr, x.m2)
  for (i in 1:(itr - 1)) {
    if (y.l$value < y.m1$value & y.m1$value < y.m2$value & y.m2$value < y.r$value) {
      x.r <- x.r + 2 * (x.r - x.l)
      y.r <- eval.loose(obj, markers, weights, expr, x.r)
    } else if (y.l$value > y.m1$value & y.m1$value > y.m2$value & y.m2$value > y.r$value) {
      x.l <- x.l - 2 * (x.r - x.l)
      y.l <- eval.loose(obj, markers, weights, expr, x.l)
    } else if (y.m1$value < y.m2$value) {
      x.l <- x.m1
      y.l <- eval.loose(obj, markers, weights, expr, x.l)
    } else {
      x.r <- x.m2
      y.r <- eval.loose(obj, markers, weights, expr, x.r)
    }
    x.m1 <- x.l + (x.r - x.l) / 3
    x.m2 <- x.r - (x.r - x.l) / 3
    y.m1 <- eval.loose(obj, markers, weights, expr, x.m1)
    y.m2 <- eval.loose(obj, markers, weights, expr, x.m2)
  }
  if (y.m1$value < y.m2$value) {
    return(y.m2)
  } else {
    return(y.m1)
  }
}

eval.loose <- function (obj, markers, weights, expr, resolution) {
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
  return(list(resolution = resolution, value = exp(mean(log(unlist(scores[2, ])))), mapping = unlist(scores[1, ])))
}
