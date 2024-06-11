#' Remove negative control probe and other technical features from Xenium In Situ data
#'
#' @param x a SpatialFeatureExperiment or SpatialExperiment object
#'
#' @return x, but only containing gene features
#' @export
removeNegativeControlFeatures <- function(x){
  x <- x[which(rowData(x)$Type=="Gene Expression"), ]
}
