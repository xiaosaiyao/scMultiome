#'
#' scATAC-seq and unpaired scRNA-seq of hematopoetic cells
#'
#' Example scATAC-seq data of hematopoietic cells included in ArchR package was integrated with scRNAseq. ScATAC-seq data was obtained from GSE139369 and scRNA-seq obtained from https://jeffgranja.s3.amazonaws.com/ArchR/TestData/scRNA-Hematopoiesis-Granja-2019.rds
#'
#' @inheritParams prostateENZ
#'
#' @inherit prostateENZ return
#'
#' @format
#' \code{MultiAssayExperiment} obtained from an \code{ArchR} project.
#' Annotated with the hg19 genome build.
#' Contains the following experiments:
#' \itemize{
#'   \item{\strong{GeneIntegrationMatrix}: SingleCellExperiment with 17889 rows and 10250 columns}
#'   \item{\strong{GeneScoreMatrix}: SingleCellExperiment with 22217 rows and 10250 columns}
#'   \item{\strong{PeakMatrix}: SingleCellExperiment with 150046 rows and 10250 columns}
#'   \item{\strong{TileMatrix500}: SingleCellExperiment with 5762078 rows and 10250 columns}
#' }
#'
#' @references
#' Single-cell multiomic analysis identifies regulatory programs in mixed-phenotype acute leukemia
#' associated with prostate cancer relapse. \cr
#' Granja \emph{et al.}, \emph{Nature Biotechnology} 2019 Dec;37(12):1458-1465. \cr
#' \href{https://www.nature.com/articles/s41587-019-0332-7}{doi: 10.1038/s41587-019-0332-7}
#'
#' @inheritSection prostateENZ Data storage and access
#'
#' @section Data preparation:
#' ```{r child = system.file("scripts", "make-data-hematopoiesis.Rmd", package = "scMultiome")}
#' ```
#'
#' @examples
#' hematopoiesis()
#'
#' @export
#'
hematopoiesis <-
    function(metadata = FALSE,
             experiments = c("TileMatrix500",
                             "GeneScoreMatrix",
                             "GeneIntegrationMatrix",
                             "PeakMatrix")) {
        checkmate::assertFlag(metadata)
        experiments <- match.arg(experiments, several.ok = TRUE)

        retrieve("hematopoiesis", metadata, experiments, verbose = FALSE)
    }
