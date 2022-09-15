#'
#' LNCaP Cells Treated with Enzalutamide
#'
#' Single-cell ATAC sequencing of parental LNCaP cells (DMSO treated, the control),
#' LNCaP cells treated with 10µM enzalutamide for 48 hours,
#' and LNCaP-derived enzalutamide-resistant RES-A and RES-B cells.
#'
#' @param experiments character vector of matrices to return; see \code{Format}
#' @param metadata logical flag specifying whether to return data or metadata only
#'
#' @return
#' \code{MultiAssayExperiment} made up of \code{SingleCellExperiment}s
#' with assays stored as \code{DelayedMatrix} objects.
#' If \code{metadata = TRUE}, an \code{ExperimentHub} object listing this data set's metadata.
#'
#' @format
#' \code{MultiAssayExperiment} obtained from an \code{ArchR} project. Annotated with the Hg38 genome build.
#' Contains the following experiments:
#' \itemize{
#'   \item{\strong{TileMatrix500}: SingleCellAccessibilityExperiment with 6062095 rows and 15522 columns}
#'   \item{\strong{GeneScoreMatrix}: SingleCellExperiment with 24919 rows and 15522 columns}
#'   \item{\strong{GeneIntegrationMatrix}: SingleCellExperiment with 23525 rows and 15522 columns}
#'   \item{\strong{PeakMatrix}: SingleCellExperiment with 80210 rows and 15522 columns}
#'   \item{\strong{MotifMatrix}: SingleCellExperiment with 870 rows and 15522 columns}
#' }
#'
#' @references
#' Single-cell ATAC and RNA sequencing reveal pre-existing and persistent cells
#' associated with prostate cancer relapse. \cr
#' Taavitsainen \emph{et al.}, \emph{Nature Communications} 2021 Sep 6;12(1):5307 \cr
#' \href{https://pubmed.ncbi.nlm.nih.gov/34489465/}{doi: 10.1038/s41467-021-25624-1}
#'
#' @section Data storage and access:
#' The \code{MultiAssayExperiments} is split into separate \code{SingleCellExperiment}
#' objects and they in turn are split into components, all of which are stored in a
#' single hdf5 file. Data and can be accessed with a special function that extracts
#' elements of the requested experiment(s), reassembles them, and builds an MAE.
#'
#' @section Data preparation:
#' ```{r child = system.file("scripts", "make-data-prostateENZ.Rmd", package = "scMultiome")}
#' ```
#'
#' @export
#'
prostateENZ <-
    function(metadata = FALSE,
             experiments = c("TileMatrix500",
                             "GeneScoreMatrix",
                             "GeneIntegrationMatrix",
                             "PeakMatrix",
                             "MotifMatrix")) {
        checkmate::assertFlag(metadata)
        experiments <- match.arg(experiments, several.ok = TRUE)

        retrieve(metadata, experiments, verbose = FALSE)
    }
