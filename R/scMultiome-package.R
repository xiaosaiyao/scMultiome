#' @details
#' Single cell multiome data sets, paired and unpaired, were analyzed
#' with the \code{ArchR} package in order to obtain epiregulons.
#' \code{ArchR} projects were converted to \code{MultiAssayExperiment} objects.
#'
#' The creation of all datasets is described in detail in respective help files.
#' Run \code{listDatasets()} to view a list of available data sets or see \code{Datasets} below.
#' See \code{?<DATASET_NAME>} for details on particular data sets, e.g. \code{?prostateENZ}.
#'
#' @section Datasets:
#' ```{r child = system.file("scripts", "datasetList.Rmd", package = "scMultiome")}
#' ```
#'
# Bioconductor requires that AnnotationHub and ExperimentHub be fully imported
#' @import AnnotationHub
#' @import ExperimentHub
#'
#' @keywords internal
"_PACKAGE"
