#' A wrapper around \code{scuttle::addPerCellQCMetrics} to add Xenium-specific QC metrics based on negative control probe counts.
#'
#' @param x a SpatialFeatureExperiment or SpatialExperiment object containing data from the Xenium In Situ platform
#'
#' @return x is returned with QC metrics added to \code{colData}
#' @export
addPerCellXeniumQCMetrics <- function(x){

  # Add some QC Metrics
  colData(x)$nCounts <- Matrix::colSums(counts(x))
  colData(x)$nGenes <- Matrix::colSums(counts(x) > 0)

  is_blank <- stringr::str_detect(rownames(x), "^BLANK_")
  is_neg <- stringr::str_detect(rownames(x), "^NegControlProbe")
  is_neg2 <- stringr::str_detect(rownames(x), "^NegControlCodeword")
  is_anti <- stringr::str_detect(rownames(x), "^antisense")
  is_depr <- stringr::str_detect(rownames(x), "^DeprecatedCodeword")
  is_unassigned <- stringr::str_detect(rownames(x), "^Unassigned")

  is_any_neg <- is_blank | is_neg | is_neg2 | is_anti | is_depr | is_unassigned
  rowData(x)$is_neg <- is_any_neg

  n_panel <- nrow(x) - sum(is_any_neg)

  # normalize counts after QC
  colData(x)$nCounts_normed <- x$nCounts/n_panel
  colData(x)$nGenes_normed <- x$nGenes/n_panel
  colData(x)$prop_nuc <- x$nucleus_area / x$cell_area

  # add QC columns
  x <- scuttle::addPerCellQCMetrics(x, subsets = list(blank = is_blank,
                                                          negProbe = is_neg,
                                                          negCodeword = is_neg2,
                                                          anti = is_anti,
                                                          depr = is_depr,
                                                          unassigned = is_unassigned,
                                                          any_neg = is_any_neg))

  return(x)
}
