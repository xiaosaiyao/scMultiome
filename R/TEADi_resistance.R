#'
#' TEAD inhibitor resistance
#'
#' Single cell multiomics on cell lines resistant to TEAD inhibitor
#'
#' @inheritParams prostateENZ
#'
#' @inherit prostateENZ return
#'
#' @format
#' \code{MultiAssayExperiment} obtained from an \code{ArchR} project.
#' Annotated with the Hg38 genome build.
#' Contains the following experiments:
#' \itemize{
#'   \item{\strong{TileMatrix}: SingleCellExperiment with 6068436 rows and 4952 columns}
#'   \item{\strong{GeneScoreMatrix}: SingleCellExperiment with 57765 rows and 4952 columns}
#'   \item{\strong{GeneExpressionMatrix}: SingleCellExperiment with 36451 rows and 4952 columns}
#'   \item{\strong{PeakMatrix}: SingleCellExperiment with 103723 rows and 4952 column}
#'   \item{\strong{MotifMatrix}: SingleCellExperiment with 870 rows and 4952 columns}
#'   \item{\strong{TFPeaksDeviationsMatrix}: SingleCellExperiment with 1269 rows and 4952 columns}
#'   \item{\strong{TF_bindingMatrix}: SingleCellExperiment with 1504 rows and 4952 columns}
#' }
#'
#' @references
#' Manuscript under review...
#'
#' @inheritSection prostateENZ Data storage and access
#'
#' @section Data preparation:
#' ```{r child = system.file("scripts", "make-data-TEADi_resistance.Rmd", package = "scMultiome")}
#' ```
#'
#' @examples
#' TEADi_resistance()
#'
#' @export
#'
TEADi_resistance <-
    function(metadata = FALSE,
             experiments = c("TileMatrix",
                             "GeneScoreMatrix",
                             "GeneIntegrationMatrix",
                             "PeakMatrix",
                             "MotifMatrix",
                             "TFPeaksDeviationsMatrix",
                             "TF_bindingMatrix")) {
        checkmate::assertFlag(metadata)
        experiments <- match.arg(experiments, several.ok = TRUE)

        retrieve("TEADi_resistance", metadata, experiments, verbose = FALSE)
    }
