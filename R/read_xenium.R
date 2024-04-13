#' Read in Xenium data and create an SFE object
#'
#' @param dir_name the name of the directory containing Xenium output files
#'
#' @return a SpatialFeatureExperiment object
#' @import SingleCellExperiment SpatialExperiment SpatialFeatureExperiment S4Vectors SummarizedExperiment
#' @export
readXenium <- function(dir_name){ # Note: this code for this function is based off of this resource: https://github.com/pachterlab/SFEData/blob/main/inst/scripts/make-data.R
  # but to my knowledge the code in the link is not available as an exported function, hence, my package.
  counts_path <- here::here(dir_name, "cell_feature_matrix.h5")
  cell_info_path <- here::here(dir_name, "cells.csv.gz")
  cell_poly_path <- here::here(dir_name, "cell_boundaries.csv.gz")
  nuc_poly_path <- here::here(dir_name, "nucleus_boundaries.csv.gz")

  # Read in the counts
  sce <- DropletUtils::read10xCounts(counts_path)
  counts(sce) <- methods::as(DelayedArray::realize(counts(sce)), "dgCMatrix")

  # read in cell info and cell and nucleus polygons
  cell_info <- vroom::vroom(cell_info_path)
  cell_poly <- vroom::vroom(cell_poly_path)
  nuc_poly <- vroom::vroom(nuc_poly_path)

  names(cell_poly)[1] <- "ID"
  names(nuc_poly)[1] <- "ID"

  # change df to sf object for cell/nuc images
  cells_sf <- df2sf(cell_poly, c("vertex_x", "vertex_y"), geometryType = "POLYGON")
  nuc_sf <- df2sf(nuc_poly, c("vertex_x", "vertex_y"), geometryType = "POLYGON")

  # QC check
  all(sf::st_is_valid(cells_sf))
  all(sf::st_is_valid(nuc_sf))

  # get rid if invalid cells/nucs
  ind_invalid <- !sf::st_is_valid(nuc_sf)
  nuc_sf[ind_invalid,] <- nngeo::st_remove_holes(sf::st_buffer(nuc_sf[ind_invalid,], 0))

  # add cell info
  colData(sce) <- cbind(colData(sce), cell_info)


  # make spatial objects
  spe <- toSpatialExperiment(sce, spatialCoordsNames = c("x_centroid", "y_centroid"))
  sfe <- toSpatialFeatureExperiment(spe)

  # add segmented cells/nuc to spatial object
  cellSeg(sfe, withDimnames = FALSE) <- cells_sf
  nucSeg(sfe, withDimnames = FALSE) <- nuc_sf

  # Add some QC Metrics
  colData(sfe)$nCounts <- Matrix::colSums(counts(sfe))
  colData(sfe)$nGenes <- Matrix::colSums(counts(sfe) > 0)

  is_blank <- stringr::str_detect(rownames(sfe), "^BLANK_")
  is_neg <- stringr::str_detect(rownames(sfe), "^NegControlProbe")
  is_neg2 <- stringr::str_detect(rownames(sfe), "^NegControlCodeword")
  is_anti <- stringr::str_detect(rownames(sfe), "^antisense")
  is_depr <- stringr::str_detect(rownames(sfe), "^DeprecatedCodeword")
  is_unassigned <- stringr::str_detect(rownames(sfe), "^Unassigned")

  is_any_neg <- is_blank | is_neg | is_neg2 | is_anti | is_depr | is_unassigned
  rowData(sfe)$is_neg <- is_any_neg

  n_panel <- nrow(sfe) - sum(is_any_neg)

  # normalize counts after QC
  colData(sfe)$nCounts_normed <- sfe$nCounts/n_panel
  colData(sfe)$nGenes_normed <- sfe$nGenes/n_panel
  colData(sfe)$prop_nuc <- sfe$nucleus_area / sfe$cell_area

  # add QC columns
  sfe <- scuttle::addPerCellQCMetrics(sfe, subsets = list(blank = is_blank,
                                                 negProbe = is_neg,
                                                 negCodeword = is_neg2,
                                                 anti = is_anti,
                                                 depr = is_depr,
                                                 unassigned = is_unassigned,
                                                 any_neg = is_any_neg))

  # Add cell ids and make gene names unique
  colnames(sfe) <- seq_len(ncol(sfe))
  rownames(sfe) <- scuttle::uniquifyFeatureNames(ID=rownames(sfe),  names=rowData(sfe)$Symbol)

  return(sfe)
}
