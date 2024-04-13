#' Convert a SpatialFeatureExperiment object to a SpatialExperimentObject
#'
#' @param sfe a SpatialFeatureExperiment object
#'
#' @return a SpatialExperiment Object
#' @export
#' @import SingleCellExperiment
SFEtoSPE <- function(sfe){
  # First convert to an SCE
  sce <- SingleCellExperiment(assays=list(counts=counts(sfe)))
  colData(sce) <- colData(sfe)
  rowData(sce) <- rowData(sfe)
  rownames(sce) <- rownames(sfe)
  metadata(sce) <- metadata(sfe)
  reducedDims(sce) <- reducedDims(sfe)
  colnames(sce) <- colnames(sfe)

  # Get the x and y centroids for this region, and add them to the coldata
  sce$x_centroid <- spatialCoords(sfe)[,1]
  sce$y_centroid <- spatialCoords(sfe)[,2]

  # Convert to SPE
  spe <- toSpatialExperiment(sce, spatialCoordsNames=c("x_centroid", "y_centroid"))

  return(spe)
}
