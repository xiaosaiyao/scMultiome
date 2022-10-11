#'
#' reprogramSeq
#'
#' scMultiome data of LNCaP infected with FOXA1, NKX2-1, GATA6
#'
#' @inheritParams prostateENZ
#'
#' @inherit prostateENZ return
#'
#' @format
#' \code{MultiAssayExperiment} obtained from an \code{ArchR} project.
#' Annotated with the hg38 genome build.
#' Contains the following experiments:
#' \itemize{
#'   \item{\strong{TileMatrix500}: SingleCellAccessibilityExperiment with 6062095 rows and 3903 columns}
#'   \item{\strong{GeneExpressionMatrix}: SingleCellExperiment with 36438 rows and 3903 columns}
#'   \item{\strong{GeneScoreMatrix}: SingleCellExperiment with 24919 rows and 3903 columns}
#'   \item{\strong{NEPCMatrix}: SingleCellExperiment with 2 rows and 3903 columns}
#'   \item{\strong{PeakMatrix}: SingleCellExperiment with 126602 rows and 3903 columns}
#'   \item{\strong{TF_bindingMatrix}: SingleCellExperiment with 1274 rows and 3903 columns}
#' }
#'
#' @references
#' Genentech dataset
#'
#' @inheritSection prostateENZ Data storage and access
#'
#' @section Data preparation:
#' ```{r child = system.file("scripts", "make-data-reprogramSeq.Rmd", package = "scMultiome")}
#' ```
#'
#' @examples
#' reprogramSeq()
#'
#' @export
#'
#'
reprogramSeq <-
    function(metadata = FALSE,
             experiments = c("TileMatrix500",
                             "GeneExpressionMatrix",
                             "GeneScoreMatrix",
                             "NEPCMatrix",
                             "PeakMatrix",
                             "TF_bindingMatrix")) {
        checkmate::assertFlag(metadata)
        experiments <- match.arg(experiments, several.ok = TRUE)

retrieve("reprogramSeq", metadata, experiments, verbose = FALSE)
    }

# place methods to convert your custom class to SingleCellExperiment and vice versa here
# skip if conversion can be done with `as`
# see ?convertSCE and the Adding Data Sets vignette (Custom Classes) for details
