#'
#' 10k PBMC data
#'
#' PBMC from a Healthy Donor. Granulocytes were removed by cell sorting.
#' Paired ATAC and Gene Expression libraries were generated from the isolated nuclei.
#' Targeted nuclei recovery was 10,000.
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
#'   \item{\strong{TileMatrix}: SingleCellExperiment with 6062095 rows and 9702 columns}
#'   \item{\strong{GeneScoreMatrix}: SingleCellExperiment with 24919 rows and 9702 columns}
#'   \item{\strong{GeneExpressionMatrix}: SingleCellExperiment with 36438 rows and 9702 columns}
#'   \item{\strong{PeakMatrix}: SingleCellExperiment with 159290 rows and 9702 columns}
#'   \item{\strong{MotifMatrix}: SingleCellExperiment with 870 rows and 9702 columns}
#' }
#'
#' @inheritSection prostateENZ Data storage and access
#'
#' @section Data preparation:
#' ```{r child = system.file("scripts", "make-data-PBMC_10x.Rmd", package = "scMultiome")}
#' ```
#'
#' @examples
#' PBMC_10x()
#'
#' @export
#'
PBMC_10x <-
    function(metadata = FALSE,
             experiments = c("TileMatrix",
                             "GeneScoreMatrix",
                             "GeneIntegrationMatrix",
                             "PeakMatrix",
                             "MotifMatrix")) {
        checkmate::assertFlag(metadata)
        experiments <- match.arg(experiments, several.ok = TRUE)

        retrieve("PBMC_10x", metadata, experiments, verbose = FALSE)
    }
