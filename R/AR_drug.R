#'
#' AR-dependent gene expression in prostate cancer cells
#'
#' Single cell gene expression and ATACseq from prostate cancer cell lines (LNCaP, VCaP, DU145, MDA-PCA-2B, 22Rv1 and NCI-H660).
#' after 24h of drug treatment or DMSO control.
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
#'   \item{\strong{TileMatrix}: SingleCellExperiment with 6068436 rows and 23118 columns}
#'   \item{\strong{GeneScoreMatrix}: SingleCellExperiment with 57765 rows and 23118 columns}
#'   \item{\strong{GeneExpressionMatrix}: SingleCellExperiment with 36451 rows and 23118 columns}
#'   \item{\strong{PeakMatrix}: SingleCellExperiment with 237856 rows and 23118 columns}
#'   \item{\strong{MotifMatrix}: SingleCellExperiment with 870 rows and 23118 columns}
#'   \item{\strong{TFPeaksDeviationsMatrix}: SingleCellExperiment with 1533 rows and 23118 columns}
#' }
#'
#' @references
#' Tomasz Włodarczyk, Aaron Lun, Diana Wu, Shreya Menon, Shushan Toneyan, Kerstin Seidel, Liang Wang,
#' Jenille Tan, Shang-Yang Chen, Timothy Keyes, Aleksander Chlebowski, Yu Guo, Ciara Metcalfe,
#' Marc Hafner, Christian W. Siebel, M. Ryan Corces, Robert Yauch, Shiqi Xie, Xiaosai Yao. 2023.
#' "Inference of single-cell transcription factor activity to dissect mechanisms of lineage
#' plasticity and drug response" bioRxiv 2023.11.27.568955; doi: https://doi.org/10.1101/2023.11.27.568955
#'
#' @inheritSection prostateENZ Data storage and access
#'
#' @section Data preparation:
#' ```{r child = system.file("scripts", "make-data-AR_drug.Rmd", package = "scMultiome")}
#' ```
#'
#' @examples
#' AR_drug()
#'
#' @export
#'
AR_drug <-
    function(metadata = FALSE,
             experiments = c("TileMatrix",
                             "GeneScoreMatrix",
                             "GeneIntegrationMatrix",
                             "PeakMatrix",
                             "MotifMatrix",
                             "TFPeaksDeviationsMatrix")) {
        checkmate::assertFlag(metadata)
        experiments <- match.arg(experiments, several.ok = TRUE)

        retrieve("AR_drug", metadata, experiments, verbose = FALSE)
    }
